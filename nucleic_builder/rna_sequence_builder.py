"""Pinned ideal A-form dsRNA construction with AmberClassic NAB.

This module only supplies atomic coordinates.  The generated PDB remains
private and is passed through the same strict PDB validator and pinned Martini
converter used for user-supplied structures.
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


def build_rna_from_sequence(
    sequence: str,
    name: str,
    output_dir: str | Path,
    *,
    force: bool = False,
    elastic_network: str,
    verbose: bool = False,
    amberclassic_home: str | Path | None = None,
) -> BuildResult:
    """Build molecule-level Martini files for an ideal canonical A-form duplex."""

    duplex = validate_duplex_sequence(sequence)
    with tempfile.TemporaryDirectory(prefix="nucleic-builder-rna-sequence-") as temporary:
        work = Path(temporary)
        pdb, checksum = _generate_atomic_duplex(
            duplex,
            work,
            amberclassic_home=amberclassic_home,
            verbose=verbose,
        )
        provenance = GeneratedStructureProvenance(
            sequence_5to3=duplex.strand_a_5to3,
            complement_5to3=duplex.strand_b_5to3,
            paired_complement_3to5=duplex.strand_b_paired_3to5,
            generator_name="AmberClassic NAB",
            generator_version=AMBERCLASSIC_TAG,
            generator_repository=AMBERCLASSIC_REPOSITORY,
            generator_commit=AMBERCLASSIC_COMMIT,
            generator_settings=NAB_SETTINGS,
            intermediate_pdb_sha256=checksum,
        )
        return build_rna(
            pdb,
            name,
            output_dir,
            force=force,
            elastic_network=elastic_network,
            verbose=verbose,
            _generated_structure=provenance,
        )
