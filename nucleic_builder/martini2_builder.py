"""Separate backend for the official Martini 2 DNA/RNA model."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

from .core import (
    ELASTIC_NETWORK_POLICIES,
    BuildResult,
    GeneratedStructureProvenance,
    ITPData,
    _CGAtom,
    _ElasticNetworkStats,
    _NAME_PATTERN,
    _chain_provenance,
    _elastic_network_stats,
    _parse_cg_pdb,
    _publish_output_pair,
    _remove_cross_chain_elastic_bonds,
    _validate_coordinate_units,
    _write_gro,
    parse_itp,
    validate_outputs,
)
from .dna_pdb_input import InputSummary as DNAInputSummary
from .dna_pdb_input import prepare_dna_input_pdb
from .errors import ConversionError, OutputExistsError, OutputValidationError
from .rna_pdb_input import InputSummary as RNAInputSummary
from .rna_pdb_input import prepare_input_pdb


MARTINI2_VENDOR_ROOT = Path(__file__).parent / "_vendor" / "martini_2_nucleic"
MARTINI2_CONVERTER = MARTINI2_VENDOR_ROOT / "martinize_nucleotide_py3.py"
MARTINI2_FORCE_FIELD = MARTINI2_VENDOR_ROOT / "martini_v2.1-dna.itp"
MARTINI2_POLARIZABLE_FORCE_FIELD = MARTINI2_VENDOR_ROOT / "martini_v2.1P-dna.itp"
MARTINI2_IONS = MARTINI2_VENDOR_ROOT / "martini_v2.0_ions.itp"
MARTINI2_FORCE_FIELD_FILES = (
    MARTINI2_FORCE_FIELD,
    MARTINI2_IONS,
)
MARTINI2_PACKAGED_PARAMETER_FILES = (
    MARTINI2_FORCE_FIELD,
    MARTINI2_POLARIZABLE_FORCE_FIELD,
    MARTINI2_IONS,
)

MARTINI2_SOURCE_PAGE = (
    "https://cgmartini.nl/docs/downloads/force-field-parameters/"
    "martini2/nucleic_acids.html"
)
MARTINI2_ARCHIVE = "na-tutorials_20170815.tar"
MARTINI2_ARCHIVE_SHA256 = (
    "15ba5bf45b9890603f0113d2021074f397a7f8d0264cb2093e970198f4b6c20b"
)
MARTINI2_CONVERTER_VERSION = "2.2 (2017 DNA/RNA release)"
MARTINI2_CONVERTER_SHA256 = (
    "e02a0ede1f444ccbd7fc9a7e2c0ee6910642887490210bf2f7c1076a2cce3edb"
)
MARTINI2_PORT_SHA256 = (
    "ee858476b4e09e0f13d0131ed9b2f617792ca87142af8364daa453913eb8e9fd"
)
MARTINI2_FORCE_FIELD_SHA256 = (
    "cc7c200dff400e97311213b93127697c6f8c21edb2350926072f0194eb90efe6"
)
MARTINI2_POLARIZABLE_FORCE_FIELD_SHA256 = (
    "b8dea4ffbef3a439db0baa465825528db919b026b1a821b374f1a6a605912ff0"
)
MARTINI2_IONS_SHA256 = (
    "c5b9b5b9541aa6d77b5b41a4b19dee62c1b8631c73e79cf76e1b281a144b4b4e"
)
MARTINI2_MODEL_WARNING = (
    "LEGACY MARTINI 2 DNA/RNA MODEL; USE ITS VERSION-MATCHED FORCE FIELD "
    "AND OBSERVE THE PUBLISHED LIMITATIONS"
)

_EXPECTED_BEADS = {
    "rna": {"A": 7, "C": 6, "G": 7, "U": 6},
    "dna": {"A": 7, "C": 6, "G": 7, "T": 6},
}


def _prepare_dna_names(path: Path) -> None:
    """Give the official Martini 2 converter the DA/DC/DG/DT names it requires."""

    rewritten: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line[0:6].strip() in {"ATOM", "HETATM"}:
            base = line[17:20].strip()
            if base in {"A", "C", "G", "T"}:
                line = line[:17] + f"{'D' + base:>3}" + line[20:]
        rewritten.append(line)
    path.write_text("\n".join(rewritten) + "\n", encoding="utf-8")


def _run_converter(
    work: Path,
    normalized_pdb: Path,
    cg_pdb: Path,
    *,
    name: str,
    topology_type: str,
    verbose: bool,
) -> Path:
    command = [
        sys.executable,
        str(MARTINI2_CONVERTER),
        "-type",
        topology_type,
        "-f",
        str(normalized_pdb),
        "-o",
        str(work / "converter-system.top"),
        "-x",
        str(cg_pdb),
        "-name",
        name,
    ]
    environment = os.environ.copy()
    environment.setdefault("PYTHONWARNINGS", "ignore::SyntaxWarning")
    try:
        completed = subprocess.run(
            command,
            cwd=work,
            env=environment,
            text=True,
            capture_output=True,
            check=False,
            timeout=240,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise ConversionError(f"Official Martini 2 converter failed: {exc}") from exc
    if verbose:
        if completed.stdout:
            print(completed.stdout, end="")
        if completed.stderr:
            print(completed.stderr, end="", file=sys.stderr)
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        if len(detail) > 4000:
            detail = detail[-4000:]
        raise ConversionError(
            "Official Martini 2 converter failed with exit code "
            f"{completed.returncode}: {detail or 'no diagnostic output'}"
        )
    generated = sorted(work.glob("*.itp"))
    if len(generated) != 1:
        names = ", ".join(path.name for path in generated) or "none"
        raise ConversionError(
            "Official Martini 2 converter did not produce exactly one merged molecule "
            f"ITP (found: {names})"
        )
    if not cg_pdb.is_file():
        raise ConversionError("Official Martini 2 converter generated no CG coordinates")
    return generated[0]


def _normalize_itp_atoms_and_name(
    path: Path,
    name: str,
    *,
    topology_type: str,
) -> None:
    """Normalize wrapper metadata without changing any Martini 2 bead type."""

    current_section: str | None = None
    molecule_written = False
    output: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if stripped.startswith("; Using the following options:"):
            output.append(
                "; Upstream converter options: "
                f"-type {topology_type} -name {name}; temporary paths managed by "
                "nucleic_builder"
            )
            continue
        if stripped.startswith("[") and stripped.endswith("]"):
            current_section = stripped[1:-1].strip().lower()
            output.append(raw_line)
            continue
        if not stripped or stripped.startswith((";", "#")):
            output.append(raw_line)
            continue
        data, separator, comment = raw_line.partition(";")
        fields = data.split()
        if current_section == "moleculetype" and not molecule_written:
            fields[0] = name
            molecule_written = True
            raw_line = "  ".join(fields)
            if separator:
                raw_line += f" ;{comment}"
        elif current_section == "atoms" and len(fields) >= 8:
            # Upstream omits each 5' phosphate coordinate but leaves a gap in
            # cgnr. Every bead already has its own charge group, so making cgnr
            # equal to the delivered atom index is semantics-preserving and
            # keeps the result inside the wrapper's strict validation boundary.
            fields[5] = fields[0]
            raw_line = " ".join(fields)
            if separator:
                raw_line += f" ;{comment}"
        output.append(raw_line)
    if not molecule_written:
        raise ConversionError("Official Martini 2 ITP has no [ moleculetype ] record")
    path.write_text("\n".join(output) + "\n", encoding="utf-8")


def _remove_all_elastic_bonds(path: Path) -> int:
    current_section: str | None = None
    output: list[str] = []
    removed = 0
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            current_section = stripped[1:-1].strip().lower()
        if current_section == "bonds" and stripped and not stripped.startswith((";", "#")):
            fields = raw_line.split(";", 1)[0].split()
            try:
                function = int(fields[2])
            except (ValueError, IndexError):
                pass
            else:
                if function == 6:
                    removed += 1
                    continue
        output.append(raw_line)
    path.write_text("\n".join(output) + "\n", encoding="utf-8")
    return removed


def _activate_rubber_bands(path: Path) -> None:
    """Activate this molecule's official network without leaking its macro."""

    lines = path.read_text(encoding="utf-8").splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if line.strip() == "#ifdef RUBBER_BANDS")
    except StopIteration as exc:
        raise ConversionError("Official Martini 2 ITP has no rubber-band block") from exc

    depth = 0
    end: int | None = None
    for index in range(start, len(lines)):
        directive = lines[index].strip().split(maxsplit=1)[0] if lines[index].strip() else ""
        if directive in {"#if", "#ifdef", "#ifndef"}:
            depth += 1
        elif directive == "#endif":
            depth -= 1
            if depth == 0:
                end = index
                break
    if end is None:
        raise ConversionError("Official Martini 2 rubber-band block is unterminated")

    before = [
        "#ifndef RUBBER_BANDS",
        "#define NUCLEIC_BUILDER_M2_RUBBER_BANDS",
        "#define RUBBER_BANDS",
        "#endif",
    ]
    after = [
        "#ifdef NUCLEIC_BUILDER_M2_RUBBER_BANDS",
        "#undef RUBBER_BANDS",
        "#undef NUCLEIC_BUILDER_M2_RUBBER_BANDS",
        "#endif",
    ]
    lines[start:start] = before
    end += len(before)
    lines[end + 1 : end + 1] = after
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _validate_mapping(
    cg_atoms: tuple[_CGAtom, ...],
    itp: ITPData,
    summary: RNAInputSummary | DNAInputSummary,
    *,
    polymer: str,
    elastic_network: str,
) -> None:
    expected_count = (
        sum(_EXPECTED_BEADS[polymer][res.canonical_name] for res in summary.residues)
        - summary.chain_count
    )
    if len(cg_atoms) != expected_count or len(itp.atoms) != expected_count:
        raise ConversionError(
            "Incomplete official Martini 2 atom-to-bead mapping: "
            f"coordinates={len(cg_atoms)}, topology={len(itp.atoms)}, "
            f"expected={expected_count}"
        )
    if [atom.index for atom in cg_atoms] != list(range(1, expected_count + 1)):
        raise ConversionError("Martini 2 CG PDB bead indices are not contiguous from 1")

    residue_by_number = {res.output_number: res for res in summary.residues}
    for cg_atom, itp_atom in zip(cg_atoms, itp.atoms):
        residue = residue_by_number.get(cg_atom.residue_number)
        if residue is None:
            raise ConversionError(
                f"Martini 2 CG PDB uses unexpected residue {cg_atom.residue_number}"
            )
        expected_name = (
            residue.canonical_name if polymer == "rna" else f"D{residue.canonical_name}"
        )
        if cg_atom.chain_id != residue.normalized_chain_id:
            raise ConversionError(
                f"Chain identity changed at Martini 2 bead {cg_atom.index}: "
                f"{cg_atom.chain_id!r} != {residue.normalized_chain_id!r}"
            )
        cg_identity = (cg_atom.residue_number, cg_atom.residue_name, cg_atom.atom_name)
        itp_identity = (
            itp_atom.residue_number,
            itp_atom.residue_name,
            itp_atom.atom_name,
        )
        expected_identity = (cg_atom.residue_number, expected_name, cg_atom.atom_name)
        if cg_identity != expected_identity or itp_identity != expected_identity:
            raise ConversionError(
                f"Martini 2 coordinate/topology identity mismatch at bead {cg_atom.index}: "
                f"CG={cg_identity!r}, ITP={itp_identity!r}, expected={expected_identity!r}"
            )

    chain_by_atom = {atom.index: atom.chain_id for atom in cg_atoms}
    for interaction in itp.interactions:
        chains = {chain_by_atom[index] for index in interaction.indices}
        if len(chains) <= 1:
            continue
        if (
            elastic_network == "legacy"
            and interaction.section == "bonds"
            and interaction.function == 6
        ):
            continue
        raise ConversionError(
            f"Martini 2 [{interaction.section}] interaction crosses a chain boundary "
            f"at ITP line {interaction.line_number}"
        )


def _sequence_provenance(
    generated: GeneratedStructureProvenance,
    *,
    polymer: str,
) -> list[str]:
    if polymer == "rna" and generated.strand_mode == "single":
        lines = [
            "; Input mode: sequence-derived idealized single-stranded RNA",
            "; Strand mode: single",
            f"; Strand (5'->3'): {generated.sequence_5to3}",
            "; Initial conformation: strand A retained from an ideal NAB A-form duplex",
            "; Structure prediction: none; folded solution structure was not predicted",
        ]
    else:
        form = "A-form dsRNA" if polymer == "rna" else "B-form dsDNA"
        lines = [
            f"; Input mode: sequence-derived ideal {form}",
            f"; Strand A (5'->3'): {generated.sequence_5to3}",
            f"; Strand B (5'->3'): {generated.complement_5to3}",
            f"; Strand B paired orientation (3'->5'): {generated.paired_complement_3to5}",
            "; Duplex orientation: strand A 5'->3'; strand B antiparallel 3'->5'",
        ]
    lines.extend(
        [
            f"; Structure generator: {generated.generator_name} {generated.generator_version}",
            f"; Structure generator repository: {generated.generator_repository}",
            f"; Structure generator commit: {generated.generator_commit}",
            f"; Structure generator settings: {generated.generator_settings}",
            f"; Intermediate PDB SHA-256: {generated.intermediate_pdb_sha256}",
        ]
    )
    return lines


def _add_provenance(
    path: Path,
    source_name: str,
    summary: RNAInputSummary | DNAInputSummary,
    *,
    polymer: str,
    topology_type: str,
    elastic_network: str,
    elastic_stats: _ElasticNetworkStats,
    generated_structure: GeneratedStructureProvenance | None,
) -> None:
    citation = (
        "Uusitalo et al., Biophys. J. 113, 246-256 (2017), "
        "DOI 10.1016/j.bpj.2017.05.043"
        if polymer == "rna"
        else "Uusitalo et al., JCTC 11, 3932-3945 (2015), "
        "DOI 10.1021/acs.jctc.5b00286"
    )
    limitation = (
        "not for RNA folding, hybridization, melting, hairpin formation, or "
        "intercalation; recommended maximum timestep 10 fs, sometimes smaller"
        if polymer == "rna"
        else "not for DNA hybridization, melting, hairpin formation, intercalation, "
        "or unconstrained dsDNA structural studies"
    )
    header = [
        "; Generated by nucleic_builder Martini 2 backend",
        "; Martini version: 2",
        f"; Polymer: canonical {polymer.upper()}",
        f"; Force-field source: {MARTINI2_SOURCE_PAGE}",
        f"; Force-field release: {MARTINI2_ARCHIVE} (README date 2017-06-05)",
        f"; Release archive SHA-256: {MARTINI2_ARCHIVE_SHA256}",
        "; Upstream converter status: emitted header labels the topology a "
        "development beta and says not to use it for production runs",
        f"; Converter version: {MARTINI2_CONVERTER_VERSION}",
        f"; Upstream converter SHA-256: {MARTINI2_CONVERTER_SHA256}",
        f"; Python 3 compatibility port SHA-256: {MARTINI2_PORT_SHA256}",
        "; Force-field parameters: martini_v2.1-dna.itp (standard-water model)",
        f"; Force-field parameter SHA-256: {MARTINI2_FORCE_FIELD_SHA256}",
        "; Packaged polarizable-water alternative: martini_v2.1P-dna.itp",
        "; Polarizable parameter SHA-256: "
        f"{MARTINI2_POLARIZABLE_FORCE_FIELD_SHA256}",
        "; Packaged ion definitions: martini_v2.0_ions.itp",
        f"; Ion parameter SHA-256: {MARTINI2_IONS_SHA256}",
        f"; Citation: {citation}",
        "; Upstream archive licensing: no standalone license or explicit README license found; "
        "see THIRD_PARTY_NOTICES.md",
        f"; Model limitation: {limitation}",
        f"; Official converter topology type: {topology_type}",
        f"; Elastic network policy: {elastic_network}",
        f"; Elastic bonds: {elastic_stats.total}; "
        f"cross-chain elastic bonds: {elastic_stats.cross_chain}",
        f"; Chain mapping: {_chain_provenance(summary)}",
    ]
    if generated_structure is None:
        header.extend(
            [
                f"; Input mode: all-atom {polymer.upper()} PDB",
                f"; Input: {source_name}",
                f"; Input SHA-256: {summary.input_sha256}",
            ]
        )
    else:
        header.extend(_sequence_provenance(generated_structure, polymer=polymer))
    body = path.read_text(encoding="utf-8")
    path.write_text("\n".join(header) + "\n" + body, encoding="utf-8")


def build_martini2(
    input_pdb: str | Path,
    name: str,
    output_dir: str | Path,
    *,
    polymer: str,
    force: bool = False,
    elastic_network: str = "legacy",
    verbose: bool = False,
    generated_structure: GeneratedStructureProvenance | None = None,
) -> BuildResult:
    """Build one molecule with the official Martini 2 DNA/RNA model."""

    if polymer not in {"rna", "dna"}:
        raise ConversionError("Martini 2 polymer must be 'rna' or 'dna'")
    if not _NAME_PATTERN.fullmatch(name):
        raise ConversionError(
            "Molecule name must start with a letter and contain only letters, digits, "
            "'.', '_', '+', or '-'"
        )
    if elastic_network not in ELASTIC_NETWORK_POLICIES:
        choices = ", ".join(sorted(ELASTIC_NETWORK_POLICIES))
        raise ConversionError(
            f"Unknown elastic-network policy {elastic_network!r}; choose one of: {choices}"
        )

    output_dir = Path(output_dir)
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise ConversionError(f"Could not create output directory {output_dir}: {exc}") from exc
    if not output_dir.is_dir():
        raise ConversionError(f"Output path is not a directory: {output_dir}")

    final_itp = output_dir / f"{name}.itp"
    final_gro = output_dir / f"{name}.gro"
    existing = [path for path in (final_itp, final_gro) if path.exists()]
    non_files = [path for path in existing if not path.is_file()]
    if non_files:
        raise OutputExistsError(
            "Output target exists but is not a regular file: "
            + ", ".join(str(path) for path in non_files)
        )
    if existing and not force:
        raise OutputExistsError(
            "Refusing to overwrite existing output(s): "
            + ", ".join(str(path) for path in existing)
            + ". Pass --force to replace them."
        )

    with tempfile.TemporaryDirectory(
        prefix=f".{name}.nucleic-builder-martini2-", dir=output_dir
    ) as temporary:
        work = Path(temporary)
        normalized_pdb = work / "normalized-input.pdb"
        cg_pdb = work / "coarse-grained.pdb"
        stage_itp = work / f"{name}.itp"
        stage_gro = work / f"{name}.gro"

        if polymer == "rna":
            summary: RNAInputSummary | DNAInputSummary = prepare_input_pdb(
                Path(input_pdb), normalized_pdb
            )
        else:
            summary = prepare_dna_input_pdb(Path(input_pdb), normalized_pdb)
            _prepare_dna_names(normalized_pdb)
        if (
            generated_structure is not None
            and summary.input_sha256 != generated_structure.intermediate_pdb_sha256
        ):
            raise ConversionError(
                "Private sequence-derived PDB checksum changed before Martini 2 conversion"
            )
        if summary.chain_count > 2:
            raise ConversionError(
                "Martini 2 supports one chain or exactly two chains in this pipeline; "
                f"the official merged-network converter is hard-coded for two strands, "
                f"but the input contains {summary.chain_count} chains"
            )

        if summary.chain_count == 1:
            topology_type = "ss" if elastic_network == "off" else "ss-stiff"
        else:
            topology_type = "ds-stiff"

        generated_itp = _run_converter(
            work,
            normalized_pdb,
            cg_pdb,
            name=name,
            topology_type=topology_type,
            verbose=verbose,
        )
        os.replace(generated_itp, stage_itp)
        _normalize_itp_atoms_and_name(
            stage_itp,
            name,
            topology_type=topology_type,
        )
        cg_atoms = _parse_cg_pdb(cg_pdb)

        if elastic_network == "off":
            _remove_all_elastic_bonds(stage_itp)
        elif elastic_network == "intrachain":
            _remove_cross_chain_elastic_bonds(stage_itp, cg_atoms)
            _activate_rubber_bands(stage_itp)
        else:
            _activate_rubber_bands(stage_itp)

        itp = parse_itp(stage_itp)
        elastic_stats = _elastic_network_stats(itp, cg_atoms)
        if elastic_network == "off" and elastic_stats.total:
            raise ConversionError("Martini 2 off policy left an elastic bond")
        if elastic_network == "intrachain" and elastic_stats.cross_chain:
            raise ConversionError("Martini 2 intrachain policy left a cross-chain bond")

        _validate_mapping(
            cg_atoms,
            itp,
            summary,
            polymer=polymer,
            elastic_network=elastic_network,
        )
        _add_provenance(
            stage_itp,
            Path(input_pdb).name,
            summary,
            polymer=polymer,
            topology_type=topology_type,
            elastic_network=elastic_network,
            elastic_stats=elastic_stats,
            generated_structure=generated_structure,
        )
        itp = parse_itp(stage_itp)
        translation = _write_gro(
            stage_gro,
            cg_atoms=cg_atoms,
            itp=itp,
            title=(
                f"{name}: Martini 2 {polymer.upper()}; rigidly translated, "
                "coordinates in nm"
            ),
            padding_nm=1.0,
        )
        _validate_coordinate_units(cg_atoms, stage_gro, translation)
        report = validate_outputs(
            stage_itp,
            stage_gro,
            expected_name=name,
            expected_charge=summary.expected_charge,
        )
        if report.residue_count != summary.residue_count:
            raise OutputValidationError(
                f"Martini 2 residue-count mismatch: generated {report.residue_count}, "
                f"input has {summary.residue_count}"
            )
        _publish_output_pair(stage_itp, stage_gro, final_itp, final_gro, work)

    return BuildResult(
        itp_path=final_itp,
        gro_path=final_gro,
        bead_count=report.bead_count,
        residue_count=report.residue_count,
        chain_count=summary.chain_count,
        total_charge=report.total_charge,
        ignored_water_atoms=summary.ignored_water_atoms,
        input_sha256=summary.input_sha256,
        elastic_network=elastic_network,
        elastic_bond_count=elastic_stats.total,
        cross_chain_elastic_bond_count=elastic_stats.cross_chain,
        input_mode="sequence" if generated_structure is not None else "pdb",
        sequence=(
            generated_structure.sequence_5to3 if generated_structure is not None else None
        ),
        complement=(
            generated_structure.complement_5to3 if generated_structure is not None else None
        ),
        intermediate_pdb_sha256=(
            generated_structure.intermediate_pdb_sha256
            if generated_structure is not None
            else None
        ),
        polymer_type=polymer,
        strand_mode=(
            generated_structure.strand_mode if generated_structure is not None else None
        ),
        martini_version=2,
        force_field_files=MARTINI2_FORCE_FIELD_FILES,
        backend_topology_type=topology_type,
    )
