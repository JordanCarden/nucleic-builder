"""Experimental wrapper around the separately pinned authors' DNA-alpha converter.

This module deliberately has its own converter loader and parameter path.  It
never calls the published Martini 3 RNA converter.
"""

from __future__ import annotations

import importlib.util
import logging
import sys
import tempfile
from pathlib import Path
from types import ModuleType

from .core import (
    ELASTIC_NETWORK_POLICIES,
    BuildResult,
    GeneratedStructureProvenance,
    ITPData,
    VENDOR_ROOT,
    _CGAtom,
    _ElasticNetworkStats,
    _NAME_PATTERN,
    _chain_provenance,
    _configure_upstream_logging,
    _elastic_network_stats,
    _parse_cg_pdb,
    _publish_output_pair,
    _remove_cross_chain_elastic_bonds,
    _validate_coordinate_units,
    _write_gro,
    parse_itp,
    validate_outputs,
)
from .dna_pdb_input import InputSummary, prepare_dna_input_pdb
from .errors import ConversionError, OutputExistsError, OutputValidationError


DNA_UPSTREAM_REPOSITORY = "https://github.com/DanYev/Martini-3-DNA-RNA"
DNA_UPSTREAM_COMMIT = "e761b7349fdf61dd485053c000dbb642f24ff9d8"
DNA_UPSTREAM_SCRIPT = VENDOR_ROOT / "martinize_dna_alpha.py"
DNA_UPSTREAM_PARAMETER_DIR = VENDOR_ROOT / "dna_alpha_itps"
DNA_MODEL_WARNING = "EXPERIMENTAL / UNPUBLISHED DNA-ALPHA MODEL"
_EXPECTED_DNA_BEADS = {"A": 8, "C": 7, "G": 9, "T": 8}


_dna_upstream_module: ModuleType | None = None


def _load_dna_upstream() -> ModuleType:
    """Load only the separately pinned DNA-alpha implementation."""

    global _dna_upstream_module
    if _dna_upstream_module is not None:
        return _dna_upstream_module
    spec = importlib.util.spec_from_file_location(
        "nucleic_builder._pinned_martinize_dna_alpha", DNA_UPSTREAM_SCRIPT
    )
    if spec is None or spec.loader is None:
        raise ConversionError(
            f"Could not load vendored DNA-alpha converter: {DNA_UPSTREAM_SCRIPT}"
        )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    root_logger = logging.getLogger()
    original_handlers = tuple(root_logger.handlers)
    original_level = root_logger.level
    try:
        spec.loader.exec_module(module)
    except Exception as exc:  # pragma: no cover - broken installation path
        raise ConversionError(
            f"Could not import vendored DNA-alpha converter: {exc}"
        ) from exc
    finally:
        for handler in tuple(root_logger.handlers):
            if handler not in original_handlers:
                root_logger.removeHandler(handler)
                handler.close()
        root_logger.setLevel(original_level)

    module.logger.handlers.clear()
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    module.logger.addHandler(handler)
    module.logger.propagate = False
    _dna_upstream_module = module
    return module


def _validate_dna_mapping(
    cg_atoms: tuple[_CGAtom, ...],
    itp: ITPData,
    summary: InputSummary,
    *,
    elastic_network: str,
) -> None:
    expected_count = (
        sum(_EXPECTED_DNA_BEADS[res.canonical_name] for res in summary.residues)
        - summary.chain_count
    )
    if len(cg_atoms) != expected_count:
        raise ConversionError(
            f"Incomplete DNA atom-to-bead mapping: generated {len(cg_atoms)} "
            f"coordinates, expected {expected_count}. Check the PDB for missing "
            "mapped atoms."
        )
    if len(itp.atoms) != expected_count:
        raise ConversionError(
            "Pinned DNA-alpha converter generated "
            f"{len(itp.atoms)} topology beads, expected {expected_count}"
        )
    if [atom.index for atom in cg_atoms] != list(range(1, expected_count + 1)):
        raise ConversionError("Temporary DNA-alpha CG PDB bead indices are not contiguous")

    residue_by_number = {res.output_number: res for res in summary.residues}
    for cg_atom, itp_atom in zip(cg_atoms, itp.atoms):
        expected_residue = residue_by_number.get(cg_atom.residue_number)
        if expected_residue is None:
            raise ConversionError(
                f"DNA-alpha CG PDB uses unexpected residue number {cg_atom.residue_number}"
            )
        if cg_atom.chain_id != expected_residue.normalized_chain_id:
            raise ConversionError(
                f"DNA chain identity changed at bead {cg_atom.index}: "
                f"{cg_atom.chain_id!r} != {expected_residue.normalized_chain_id!r}"
            )
        cg_identity = (cg_atom.residue_number, cg_atom.residue_name, cg_atom.atom_name)
        itp_identity = (
            itp_atom.residue_number,
            itp_atom.residue_name,
            itp_atom.atom_name,
        )
        if cg_identity != itp_identity:
            raise ConversionError(
                "Pinned DNA-alpha converter produced a coordinate/topology mismatch "
                f"at bead {cg_atom.index}: {cg_identity!r} != {itp_identity!r}"
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
            f"Non-elastic DNA [{interaction.section}] interaction crosses a chain "
            f"boundary at ITP line {interaction.line_number}"
        )


def _add_dna_itp_provenance(
    path: Path,
    source_name: str,
    summary: InputSummary,
    *,
    elastic_network: str,
    elastic_stats: _ElasticNetworkStats,
    generated_structure: GeneratedStructureProvenance | None,
) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    rewritten: list[str] = []
    inserted = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("; Generated on:"):
            continue
        if stripped == "; Created using the following options:":
            continue
        if stripped.startswith("; input_pdb="):
            continue
        rewritten.append(line)
        if not inserted and stripped.startswith("; MARTINI"):
            provenance_lines = [
                f"; !!! {DNA_MODEL_WARNING} !!!",
                "; This is not the published Martini 3 RNA model.",
                "; Generated by nucleic_builder Phase 2B",
                "; Polymer: canonical DNA",
            ]
            if generated_structure is None:
                provenance_lines.extend(
                    [
                        "; Input mode: all-atom DNA PDB",
                        f"; Input: {source_name}",
                        f"; Input SHA-256: {summary.input_sha256}",
                    ]
                )
            else:
                provenance_lines.extend(
                    [
                        "; Input mode: sequence-derived ideal B-form dsDNA",
                        f"; Strand A (5'->3'): {generated_structure.sequence_5to3}",
                        f"; Strand B (5'->3'): {generated_structure.complement_5to3}",
                        f"; Strand B paired orientation (3'->5'): "
                        f"{generated_structure.paired_complement_3to5}",
                        "; Duplex orientation: strand A 5'->3'; strand B antiparallel 3'->5'",
                        f"; Structure generator: {generated_structure.generator_name} "
                        f"{generated_structure.generator_version}",
                        f"; Structure generator repository: "
                        f"{generated_structure.generator_repository}",
                        f"; Structure generator commit: {generated_structure.generator_commit}",
                        f"; Structure generator settings: "
                        f"{generated_structure.generator_settings}",
                        f"; Intermediate PDB SHA-256: "
                        f"{generated_structure.intermediate_pdb_sha256}",
                    ]
                )
            provenance_lines.extend(
                [
                    f"; DNA-alpha upstream: {DNA_UPSTREAM_REPOSITORY}",
                    f"; DNA-alpha upstream commit: {DNA_UPSTREAM_COMMIT}",
                    "; DNA-alpha converter: martinize_dna_alpha.py",
                    "; DNA-alpha parameters: dna_alpha_itps/",
                    "; Model status: experimental alpha; DNA publication is pending upstream",
                    "; License: GPL-3.0-only; see distribution LICENSE and "
                    "nucleic_builder/_vendor/martini_3_dna_rna/NOTICE.md",
                    f"; Elastic network policy: {elastic_network}",
                    f"; Elastic bonds: {elastic_stats.total}; "
                    f"cross-chain elastic bonds: {elastic_stats.cross_chain}",
                    f"; Chain mapping: {_chain_provenance(summary)}",
                ]
            )
            rewritten.extend(provenance_lines)
            inserted = True
    if not inserted:
        raise ConversionError("Generated DNA ITP has no recognizable Martini header")
    path.write_text("\n".join(rewritten) + "\n", encoding="utf-8")


def build_dna(
    input_pdb: str | Path,
    name: str,
    output_dir: str | Path,
    *,
    force: bool = False,
    elastic_network: str = "legacy",
    verbose: bool = False,
    _generated_structure: GeneratedStructureProvenance | None = None,
) -> BuildResult:
    """Build experimental DNA-alpha ``NAME.itp`` and ``NAME.gro`` files."""

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
        raise ConversionError(
            f"Could not create output directory {output_dir}: {exc}"
        ) from exc
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
        prefix=f".{name}.nucleic-builder-dna-", dir=output_dir
    ) as tmp:
        work = Path(tmp)
        normalized_pdb = work / "normalized-dna-input.pdb"
        cg_pdb = work / "coarse-grained-dna-alpha.pdb"
        stage_itp = work / f"{name}.itp"
        stage_gro = work / f"{name}.gro"

        summary = prepare_dna_input_pdb(Path(input_pdb), normalized_pdb)
        if (
            _generated_structure is not None
            and summary.input_sha256 != _generated_structure.intermediate_pdb_sha256
        ):
            raise ConversionError(
                "Private sequence-derived DNA PDB checksum changed before conversion"
            )

        upstream = _load_dna_upstream()
        _configure_upstream_logging(upstream, verbose=verbose)
        try:
            upstream.martinize_dna(
                input_pdb=str(normalized_pdb),
                output_topology=str(stage_itp),
                output_structure=str(cg_pdb),
                molecule_name=name,
                elastic_network="no" if elastic_network == "off" else "yes",
                itp_dir=str(DNA_UPSTREAM_PARAMETER_DIR),
                debug=False,
            )
        except Exception as exc:
            raise ConversionError(
                f"Pinned upstream DNA-alpha conversion failed ({type(exc).__name__}): {exc}"
            ) from exc

        cg_atoms = _parse_cg_pdb(cg_pdb)
        itp = parse_itp(stage_itp)
        before_filter = _elastic_network_stats(itp, cg_atoms)
        if elastic_network == "intrachain":
            removed = _remove_cross_chain_elastic_bonds(stage_itp, cg_atoms)
            if removed != before_filter.cross_chain:
                raise ConversionError(
                    "Could not account for every cross-chain DNA elastic bond while "
                    "applying the intrachain policy"
                )
            itp = parse_itp(stage_itp)

        elastic_stats = _elastic_network_stats(itp, cg_atoms)
        if elastic_network == "off" and elastic_stats.total:
            raise ConversionError(
                "DNA-alpha converter generated elastic bonds while the network was disabled"
            )
        if elastic_network == "intrachain" and elastic_stats.cross_chain:
            raise ConversionError(
                "Intrachain DNA elastic-network policy left a cross-chain elastic bond"
            )

        _validate_dna_mapping(
            cg_atoms, itp, summary, elastic_network=elastic_network
        )
        _add_dna_itp_provenance(
            stage_itp,
            Path(input_pdb).name,
            summary,
            elastic_network=elastic_network,
            elastic_stats=elastic_stats,
            generated_structure=_generated_structure,
        )
        itp = parse_itp(stage_itp)
        translation = _write_gro(
            stage_gro,
            cg_atoms=cg_atoms,
            itp=itp,
            padding_nm=1.0,
            title=(
                f"{name}: {DNA_MODEL_WARNING}; rigidly translated, coordinates in nm"
            ),
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
                f"DNA residue-count mismatch: generated {report.residue_count}, "
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
        input_mode="sequence" if _generated_structure is not None else "pdb",
        sequence=(
            _generated_structure.sequence_5to3
            if _generated_structure is not None
            else None
        ),
        complement=(
            _generated_structure.complement_5to3
            if _generated_structure is not None
            else None
        ),
        intermediate_pdb_sha256=(
            _generated_structure.intermediate_pdb_sha256
            if _generated_structure is not None
            else None
        ),
        polymer_type="dna",
    )
