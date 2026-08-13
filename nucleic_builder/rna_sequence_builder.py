"""Pinned idealized RNA sequence construction with AmberClassic NAB.

This module only supplies atomic coordinates.  The generated PDB remains
private and is passed through the same strict PDB validator and pinned Martini
converter used for user-supplied structures.  Duplex mode retains the full
NAB A-RNA duplex.  Single mode retains only the entered strand from that
verified duplex geometry; it is an idealized starting conformation, not a
prediction of the strand's folded solution structure.
"""

from __future__ import annotations

import hashlib
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .core import (
    AMBERCLASSIC_COMMIT,
    AMBERCLASSIC_REPOSITORY,
    AMBERCLASSIC_TAG,
    BuildResult,
    GeneratedStructureProvenance,
    _run_backend,
    resolve_amberclassic_home,
)
from .errors import ConversionError, InputValidationError
from .rna_builder import build_rna
from .rna_pdb_input import prepare_input_pdb


NAB_SETTINGS = 'fd_helix("arna", sequence, "rna"); putpdb("-wwpdb")'
SINGLE_STRAND_NAB_SETTINGS = NAB_SETTINGS + "; retain generated strand A only"
RNA_STRAND_MODES = frozenset({"duplex", "single"})
_COMPLEMENT = str.maketrans("ACGU", "UGCA")


@dataclass(frozen=True)
class DuplexSequence:
    """Normalized strand sequences and their antiparallel interpretation."""

    strand_a_5to3: str
    strand_b_5to3: str
    strand_b_paired_3to5: str


def validate_duplex_sequence(sequence: str) -> DuplexSequence:
    """Validate one canonical RNA strand and derive its Watson-Crick partner."""

    if not isinstance(sequence, str) or not sequence:
        raise InputValidationError(
            "RNA sequence is empty; provide one strand using only A, C, G, and U"
        )
    normalized = sequence.upper()
    invalid = sorted(set(normalized) - set("ACGU"))
    if invalid:
        shown = ", ".join(repr(character) for character in invalid)
        raise InputValidationError(
            "RNA sequence contains unsupported character(s) "
            f"{shown}; Phase 2A accepts only A, C, G, and U (no DNA, "
            "modified nucleotides, whitespace, or FASTA syntax)"
        )
    paired = normalized.translate(_COMPLEMENT)
    return DuplexSequence(
        strand_a_5to3=normalized,
        strand_b_5to3=paired[::-1],
        strand_b_paired_3to5=paired,
    )


def _validate_strand_mode(strand_mode: str) -> None:
    if strand_mode not in RNA_STRAND_MODES:
        choices = ", ".join(sorted(RNA_STRAND_MODES))
        raise InputValidationError(
            f"Unknown RNA strand mode {strand_mode!r}; choose one of: {choices}"
        )


def _generate_atomic_duplex(
    sequence: DuplexSequence,
    work: Path,
    *,
    amberclassic_home: str | Path | None,
    verbose: bool,
) -> tuple[Path, str]:
    home = resolve_amberclassic_home(amberclassic_home)
    source = work / "ideal_arna_duplex.nab"
    executable = work / "ideal_arna_duplex"
    pdb = work / "ideal_arna_duplex.pdb"
    # Validation above makes interpolation safe and intentionally preserves the
    # exact user sequence in the pinned fd_helix invocation.
    source.write_text(
        "molecule duplex;\n\n"
        f'duplex = fd_helix( "arna", "{sequence.strand_a_5to3.lower()}", "rna" );\n'
        'putpdb( "ideal_arna_duplex.pdb", duplex, "-wwpdb" );\n',
        encoding="utf-8",
    )
    environment = os.environ.copy()
    environment["AMBERCLASSICHOME"] = str(home)
    environment["PATH"] = str(home / "bin") + os.pathsep + environment.get("PATH", "")
    _run_backend(
        [str(home / "bin" / "nab"), "-o", str(executable), str(source)],
        work=work,
        environment=environment,
        stage="compilation",
        verbose=verbose,
    )
    _run_backend(
        [str(executable)],
        work=work,
        environment=environment,
        stage="duplex generation",
        verbose=verbose,
    )
    if not pdb.is_file():
        raise ConversionError("AmberClassic NAB completed without writing its expected PDB")

    # Verify that fd_helix actually emitted the requested two antiparallel
    # strands before the private PDB enters the Martini converter.
    check_copy = work / "verified-duplex.pdb"
    summary = prepare_input_pdb(pdb, check_copy)
    chains = summary.residues_by_chain()
    observed = tuple("".join(res.canonical_name for res in chain) for chain in chains)
    expected = (sequence.strand_a_5to3, sequence.strand_b_5to3)
    if summary.chain_count != 2 or observed != expected:
        raise ConversionError(
            "AmberClassic NAB generated an unexpected duplex: "
            f"observed strand sequences {observed!r}, expected {expected!r}"
        )
    checksum = hashlib.sha256(pdb.read_bytes()).hexdigest()
    return pdb, checksum


def _generate_atomic_single_strand(
    sequence: DuplexSequence,
    work: Path,
    *,
    amberclassic_home: str | Path | None,
    verbose: bool,
) -> tuple[Path, str]:
    """Generate and verify one A-form-derived strand from the NAB duplex."""

    duplex_pdb, _ = _generate_atomic_duplex(
        sequence,
        work,
        amberclassic_home=amberclassic_home,
        verbose=verbose,
    )
    normalized_duplex = work / "normalized-duplex-for-single-strand.pdb"
    duplex_summary = prepare_input_pdb(duplex_pdb, normalized_duplex)
    strand_a_chain_id = duplex_summary.residues_by_chain()[0][0].normalized_chain_id

    single_pdb = work / "ideal_arna_single_strand.pdb"
    strand_lines = [
        line
        for line in normalized_duplex.read_text(encoding="utf-8").splitlines()
        if line[0:6].strip() in {"ATOM", "HETATM"}
        and line[21:22] == strand_a_chain_id
    ]
    if not strand_lines:
        raise ConversionError(
            "Could not extract strand A from the verified AmberClassic NAB duplex"
        )
    single_pdb.write_text(
        "\n".join([*strand_lines, "TER", "END"]) + "\n",
        encoding="utf-8",
    )

    verified_single = work / "verified-single-strand.pdb"
    single_summary = prepare_input_pdb(single_pdb, verified_single)
    observed = tuple(
        "".join(res.canonical_name for res in chain)
        for chain in single_summary.residues_by_chain()
    )
    expected = (sequence.strand_a_5to3,)
    if single_summary.chain_count != 1 or observed != expected:
        raise ConversionError(
            "AmberClassic NAB single-strand extraction was unexpected: "
            f"observed strand sequences {observed!r}, expected {expected!r}"
        )
    checksum = hashlib.sha256(single_pdb.read_bytes()).hexdigest()
    return single_pdb, checksum


def build_rna_from_sequence(
    sequence: str,
    name: str,
    output_dir: str | Path,
    *,
    force: bool = False,
    elastic_network: str,
    strand_mode: str = "duplex",
    verbose: bool = False,
    amberclassic_home: str | Path | None = None,
    martini_version: int | str = 3,
) -> BuildResult:
    """Build Martini files for an idealized single strand or A-form duplex.

    ``strand_mode='duplex'`` preserves the original behavior.  ``'single'``
    retains only the entered 5'->3' strand from the verified NAB A-form duplex
    geometry and does not predict its folded solution structure.
    """

    _validate_strand_mode(strand_mode)
    duplex = validate_duplex_sequence(sequence)
    with tempfile.TemporaryDirectory(prefix="nucleic-builder-rna-sequence-") as temporary:
        work = Path(temporary)
        if strand_mode == "single":
            pdb, checksum = _generate_atomic_single_strand(
                duplex,
                work,
                amberclassic_home=amberclassic_home,
                verbose=verbose,
            )
            complement_5to3 = None
            paired_complement_3to5 = None
            generator_settings = SINGLE_STRAND_NAB_SETTINGS
        else:
            pdb, checksum = _generate_atomic_duplex(
                duplex,
                work,
                amberclassic_home=amberclassic_home,
                verbose=verbose,
            )
            complement_5to3 = duplex.strand_b_5to3
            paired_complement_3to5 = duplex.strand_b_paired_3to5
            generator_settings = NAB_SETTINGS
        provenance = GeneratedStructureProvenance(
            sequence_5to3=duplex.strand_a_5to3,
            complement_5to3=complement_5to3,
            paired_complement_3to5=paired_complement_3to5,
            generator_name="AmberClassic NAB",
            generator_version=AMBERCLASSIC_TAG,
            generator_repository=AMBERCLASSIC_REPOSITORY,
            generator_commit=AMBERCLASSIC_COMMIT,
            generator_settings=generator_settings,
            intermediate_pdb_sha256=checksum,
            strand_mode=strand_mode,
        )
        return build_rna(
            pdb,
            name,
            output_dir,
            force=force,
            elastic_network=elastic_network,
            verbose=verbose,
            martini_version=martini_version,
            _generated_structure=provenance,
        )
