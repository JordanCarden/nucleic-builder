"""Pinned ideal B-form dsDNA construction with AmberClassic NAB."""

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
from .dna_builder import build_dna
from .dna_pdb_input import prepare_dna_input_pdb
from .errors import ConversionError, InputValidationError


DNA_NAB_SETTINGS = 'fd_helix("abdna", sequence, "dna"); putpdb("-wwpdb")'
_DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class DNADuplexSequence:
    """Normalized DNA strands and their antiparallel interpretation."""

    strand_a_5to3: str
    strand_b_5to3: str
    strand_b_paired_3to5: str


def validate_dna_duplex_sequence(sequence: str) -> DNADuplexSequence:
    """Accept only canonical A/C/G/T and derive the Watson-Crick partner."""

    if not isinstance(sequence, str) or not sequence:
        raise InputValidationError(
            "DNA sequence is empty; provide one strand using only A, C, G, and T"
        )
    normalized = sequence.upper()
    invalid = sorted(set(normalized) - set("ACGT"))
    if invalid:
        shown = ", ".join(repr(character) for character in invalid)
        raise InputValidationError(
            "DNA sequence contains unsupported character(s) "
            f"{shown}; experimental Phase 2B DNA accepts only A, C, G, and T "
            "(no RNA, modified nucleotides, whitespace, or FASTA syntax)"
        )
    paired = normalized.translate(_DNA_COMPLEMENT)
    return DNADuplexSequence(
        strand_a_5to3=normalized,
        strand_b_5to3=paired[::-1],
        strand_b_paired_3to5=paired,
    )


def _generate_atomic_dna_duplex(
    sequence: DNADuplexSequence,
    work: Path,
    *,
    amberclassic_home: str | Path | None,
    verbose: bool,
) -> tuple[Path, str]:
    home = resolve_amberclassic_home(amberclassic_home)
    source = work / "ideal_abdna_duplex.nab"
    executable = work / "ideal_abdna_duplex"
    pdb = work / "ideal_abdna_duplex.pdb"
    source.write_text(
        "molecule duplex;\n\n"
        f'duplex = fd_helix( "abdna", "{sequence.strand_a_5to3.lower()}", "dna" );\n'
        'putpdb( "ideal_abdna_duplex.pdb", duplex, "-wwpdb" );\n',
        encoding="utf-8",
    )
    environment = os.environ.copy()
    environment["AMBERCLASSICHOME"] = str(home)
    environment["PATH"] = str(home / "bin") + os.pathsep + environment.get("PATH", "")
    _run_backend(
        [str(home / "bin" / "nab"), "-o", str(executable), str(source)],
        work=work,
        environment=environment,
        stage="B-DNA compilation",
        verbose=verbose,
    )
    _run_backend(
        [str(executable)],
        work=work,
        environment=environment,
        stage="B-DNA duplex generation",
        verbose=verbose,
    )
    if not pdb.is_file():
        raise ConversionError(
            "AmberClassic NAB completed without writing its expected B-DNA PDB"
        )

    check_copy = work / "verified-dna-duplex.pdb"
    summary = prepare_dna_input_pdb(pdb, check_copy)
    chains = summary.residues_by_chain()
    observed = tuple("".join(res.canonical_name for res in chain) for chain in chains)
    expected = (sequence.strand_a_5to3, sequence.strand_b_5to3)
    if summary.chain_count != 2 or observed != expected:
        raise ConversionError(
            "AmberClassic NAB generated an unexpected B-DNA duplex: "
            f"observed strand sequences {observed!r}, expected {expected!r}"
        )
    checksum = hashlib.sha256(pdb.read_bytes()).hexdigest()
    return pdb, checksum


def build_dna_from_sequence(
    sequence: str,
    name: str,
    output_dir: str | Path,
    *,
    force: bool = False,
    elastic_network: str,
    verbose: bool = False,
    amberclassic_home: str | Path | None = None,
    martini_version: int | str = 3,
) -> BuildResult:
    """Build experimental DNA-alpha files for an ideal canonical B-form duplex."""

    duplex = validate_dna_duplex_sequence(sequence)
    with tempfile.TemporaryDirectory(prefix="nucleic-builder-dna-sequence-") as temporary:
        work = Path(temporary)
        pdb, checksum = _generate_atomic_dna_duplex(
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
            generator_settings=DNA_NAB_SETTINGS,
            intermediate_pdb_sha256=checksum,
        )
        return build_dna(
            pdb,
            name,
            output_dir,
            force=force,
            elastic_network=elastic_network,
            verbose=verbose,
            martini_version=martini_version,
            _generated_structure=provenance,
        )
