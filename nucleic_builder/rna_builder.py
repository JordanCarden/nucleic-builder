"""Orchestration around the pinned authors' Martini 3 RNA converter."""

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
from .errors import ConversionError, OutputExistsError, OutputValidationError
from .rna_pdb_input import InputSummary, prepare_input_pdb


UPSTREAM_REPOSITORY = "https://github.com/DanYev/Martini-3-DNA-RNA"
UPSTREAM_COMMIT = "e761b7349fdf61dd485053c000dbb642f24ff9d8"
UPSTREAM_SCRIPT = VENDOR_ROOT / "martinize_rna_v3.0.0.py"
UPSTREAM_PARAMETER_DIR = VENDOR_ROOT / "rna_v3.0.0_itps"
_EXPECTED_BEADS = {"A": 8, "C": 7, "G": 9, "U": 7}


_upstream_module: ModuleType | None = None


def _load_upstream() -> ModuleType:
    global _upstream_module
    if _upstream_module is not None:
        return _upstream_module
    spec = importlib.util.spec_from_file_location(
        "nucleic_builder._pinned_martinize_rna", UPSTREAM_SCRIPT
    )
    if spec is None or spec.loader is None:
        raise ConversionError(f"Could not load vendored converter: {UPSTREAM_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    root_logger = logging.getLogger()
    original_handlers = tuple(root_logger.handlers)
    original_level = root_logger.level
    try:
        spec.loader.exec_module(module)
    except Exception as exc:  # pragma: no cover - broken installation path
        raise ConversionError(f"Could not import vendored converter: {exc}") from exc
    finally:
        # The unmodified upstream module calls logging.basicConfig at import.
        # Undo only the handler(s) it added so importing this library does not
        # configure an application's root logger.
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
    _upstream_module = module
    return module


def _validate_upstream_mapping(
    cg_atoms: tuple[_CGAtom, ...],
    itp: ITPData,
    summary: InputSummary,
    *,
    elastic_network: str,
) -> None:
    expected_count = (
        sum(_EXPECTED_BEADS[res.canonical_name] for res in summary.residues)
        - summary.chain_count
    )
    if len(cg_atoms) != expected_count:
        raise ConversionError(
            f"Incomplete atom-to-bead mapping: generated {len(cg_atoms)} coordinates, "
            f"expected {expected_count}. Check the PDB for missing mapped atoms."
        )
    if len(itp.atoms) != expected_count:
        raise ConversionError(
            f"Pinned converter generated {len(itp.atoms)} topology beads, "
            f"expected {expected_count}"
        )
    if [atom.index for atom in cg_atoms] != list(range(1, expected_count + 1)):
        raise ConversionError("Temporary CG PDB bead indices are not contiguous from 1")

    residue_by_number = {res.output_number: res for res in summary.residues}
    for cg_atom, itp_atom in zip(cg_atoms, itp.atoms):
        expected_residue = residue_by_number.get(cg_atom.residue_number)
        if expected_residue is None:
            raise ConversionError(
                f"CG PDB uses unexpected residue number {cg_atom.residue_number}"
            )
        if cg_atom.chain_id != expected_residue.normalized_chain_id:
            raise ConversionError(
                f"Chain identity changed at bead {cg_atom.index}: "
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
                f"Pinned converter produced coordinate/topology mismatch at bead {cg_atom.index}: "
                f"{cg_identity!r} != {itp_identity!r}"
            )

    chain_by_atom = {atom.index: atom.chain_id for atom in cg_atoms}
    for interaction in itp.interactions:
        chains = {chain_by_atom[index] for index in interaction.indices}
        if len(chains) <= 1:
            continue
        # The authors' default elastic network intentionally includes
        # cross-chain distance restraints.  The legacy policy retains them;
        # no other policy or bonded term may bridge a chain break.
        if (
            elastic_network == "legacy"
            and interaction.section == "bonds"
            and interaction.function == 6
        ):
            continue
        raise ConversionError(
            f"Non-elastic [{interaction.section}] interaction crosses a chain boundary "
            f"at ITP line {interaction.line_number}"
        )


def _add_itp_provenance(
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
            if generated_structure is None:
                # Preserve the validated Phase 1.5.1 PDB-mode provenance fields.
                provenance_lines = [
                    "; Generated by nucleic_builder Phase 1.5.1",
                    f"; Input: {source_name}",
                    f"; Input SHA-256: {summary.input_sha256}",
                ]
            else:
                provenance_lines = [
                    "; Generated by nucleic_builder Phase 2A",
                    "; Input mode: sequence-derived ideal A-form dsRNA",
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
                    f"; Structure generator settings: {generated_structure.generator_settings}",
                    f"; Intermediate PDB SHA-256: "
                    f"{generated_structure.intermediate_pdb_sha256}",
                ]
            provenance_lines.extend(
                [
                    f"; Upstream: {UPSTREAM_REPOSITORY}",
                    f"; Upstream commit: {UPSTREAM_COMMIT}",
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
        raise ConversionError("Generated ITP has no recognizable Martini header")
    path.write_text("\n".join(rewritten) + "\n", encoding="utf-8")


def build_rna(
    input_pdb: str | Path,
    name: str,
    output_dir: str | Path,
    *,
    force: bool = False,
    elastic_network: str = "legacy",
    verbose: bool = False,
    _generated_structure: GeneratedStructureProvenance | None = None,
) -> BuildResult:
    """Build ``NAME.itp`` and ``NAME.gro`` for a canonical RNA PDB.

    ``elastic_network='legacy'`` preserves the Phase 1/upstream default: a
    0.3--1.2 nm network with a 200 kJ mol-1 nm-2 force constant, including
    cross-chain elastic bonds.  ``intrachain`` removes only cross-chain
    function-6 elastic bonds, and ``off`` disables the network.  Conditional
    position-restraint records remain unchanged for all policies.
    """

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
        prefix=f".{name}.nucleic-builder-rna-", dir=output_dir
    ) as tmp:
        work = Path(tmp)
        normalized_pdb = work / "normalized-input.pdb"
        cg_pdb = work / "coarse-grained.pdb"
        stage_itp = work / f"{name}.itp"
        stage_gro = work / f"{name}.gro"

        summary = prepare_input_pdb(Path(input_pdb), normalized_pdb)
        if (
            _generated_structure is not None
            and summary.input_sha256 != _generated_structure.intermediate_pdb_sha256
        ):
            raise ConversionError(
                "Private sequence-derived PDB checksum changed before Martini conversion"
            )
        upstream = _load_upstream()
        _configure_upstream_logging(upstream, verbose=verbose)
        try:
            upstream.martinize_rna(
                input_pdb=str(normalized_pdb),
                output_topology=str(stage_itp),
                output_structure=str(cg_pdb),
                molecule_name=name,
                elastic_network="no" if elastic_network == "off" else "yes",
                itp_dir=str(UPSTREAM_PARAMETER_DIR),
                debug=False,
            )
        except Exception as exc:
            raise ConversionError(
                f"Pinned upstream RNA conversion failed ({type(exc).__name__}): {exc}"
            ) from exc

        cg_atoms = _parse_cg_pdb(cg_pdb)
        itp = parse_itp(stage_itp)
        before_filter = _elastic_network_stats(itp, cg_atoms)
        if elastic_network == "intrachain":
            removed = _remove_cross_chain_elastic_bonds(stage_itp, cg_atoms)
            if removed != before_filter.cross_chain:
                raise ConversionError(
                    "Could not account for every cross-chain elastic bond while "
                    "applying the intrachain policy"
                )
            itp = parse_itp(stage_itp)

        elastic_stats = _elastic_network_stats(itp, cg_atoms)
        if elastic_network == "off" and elastic_stats.total:
            raise ConversionError(
                "Pinned converter generated elastic bonds even though the network was disabled"
            )
        if elastic_network == "intrachain" and elastic_stats.cross_chain:
            raise ConversionError(
                "Intrachain elastic-network policy left a cross-chain elastic bond"
            )

        _validate_upstream_mapping(
            cg_atoms, itp, summary, elastic_network=elastic_network
        )
        _add_itp_provenance(
            stage_itp,
            Path(input_pdb).name,
            summary,
            elastic_network=elastic_network,
            elastic_stats=elastic_stats,
            generated_structure=_generated_structure,
        )
        # Reparse after header rewriting to validate the exact delivered ITP.
        itp = parse_itp(stage_itp)
        translation = _write_gro(
            stage_gro,
            cg_atoms=cg_atoms,
            itp=itp,
            title=f"{name}: Martini 3 RNA; rigidly translated, coordinates in nm",
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
                f"Residue-count mismatch: generated {report.residue_count}, "
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
    )
