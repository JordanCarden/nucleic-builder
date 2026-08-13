from __future__ import annotations

import os
from pathlib import Path

import pytest

from nucleic_builder import build_rna_from_sequence
from nucleic_builder.core import (
    AMBERCLASSIC_COMMIT,
    AMBERCLASSIC_HOME_ENV,
    parse_itp,
    resolve_amberclassic_home,
    validate_outputs,
)
from nucleic_builder.errors import ConversionError, InputValidationError
from nucleic_builder.rna_sequence_builder import validate_duplex_sequence

from .conftest import BNT162B2_FRAGMENT_475, ROOT


def _backend_home() -> Path | None:
    value = os.environ.get(AMBERCLASSIC_HOME_ENV)
    if not value:
        return None
    path = Path(value)
    return path if path.exists() else None


def test_sequence_is_normalized_and_complemented_antiparallel() -> None:
    sequence = validate_duplex_sequence("acguu")
    assert sequence.strand_a_5to3 == "ACGUU"
    assert sequence.strand_b_paired_3to5 == "UGCAA"
    assert sequence.strand_b_5to3 == "AACGU"


@pytest.mark.parametrize(
    "sequence",
    ["", "ACGT", "ACGN", "AC GU", "AC-GU", ">strand\nACGU", "AΨGU"],
)
def test_sequence_rejects_noncanonical_or_sequence_file_syntax(sequence: str) -> None:
    with pytest.raises(InputValidationError, match="only A, C, G, and U"):
        validate_duplex_sequence(sequence)


def test_unknown_strand_mode_is_rejected_before_backend_lookup(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="Unknown RNA strand mode"):
        build_rna_from_sequence(
            "ACGU",
            "RNA",
            tmp_path,
            elastic_network="off",
            strand_mode="triplex",
        )


def test_missing_backend_has_actionable_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(AMBERCLASSIC_HOME_ENV, raising=False)
    monkeypatch.setenv("PATH", "")
    with pytest.raises(ConversionError, match="requires local AmberClassic NAB"):
        resolve_amberclassic_home()


def test_backend_revision_is_strictly_pinned(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "bin").mkdir()
    (tmp_path / "dat" / "leap").mkdir(parents=True)
    (tmp_path / "bin" / "nab").touch()
    (tmp_path / "bin" / "teLeap").touch()
    monkeypatch.setattr(
        "nucleic_builder.core._git_commit", lambda _: "0" * 40
    )
    with pytest.raises(ConversionError, match=AMBERCLASSIC_COMMIT):
        resolve_amberclassic_home(tmp_path)


@pytest.mark.skipif(
    _backend_home() is None,
    reason=(
        "sequence end-to-end test skipped: set NUCLEIC_BUILDER_AMBERCLASSIC_HOME "
        "to the pinned, built AmberClassic checkout"
    ),
)
def test_sequence_end_to_end_is_deterministic_and_private(tmp_path: Path) -> None:
    first = build_rna_from_sequence(
        "GCAUCG",
        "SEQRNA",
        tmp_path / "first",
        elastic_network="off",
    )
    second = build_rna_from_sequence(
        "gcaucg",
        "SEQRNA",
        tmp_path / "second",
        elastic_network="off",
        strand_mode="duplex",
    )

    assert first.itp_path.read_bytes() == second.itp_path.read_bytes()
    assert first.gro_path.read_bytes() == second.gro_path.read_bytes()
    checked_in = ROOT / "examples" / "sequence" / "generated"
    assert (checked_in / "SEQRNA.itp").read_bytes() == first.itp_path.read_bytes()
    assert (checked_in / "SEQRNA.gro").read_bytes() == first.gro_path.read_bytes()
    assert sorted(path.name for path in first.itp_path.parent.iterdir()) == [
        "SEQRNA.gro",
        "SEQRNA.itp",
    ]
    assert first.input_mode == "sequence"
    assert first.sequence == "GCAUCG"
    assert first.complement == "CGAUGC"
    assert first.intermediate_pdb_sha256 == (
        "8ce8a436b8d2da812de01dc23685fc7b9aee5bca0918891bd4854900bc6b3a43"
    )
    assert first.bead_count == 92
    assert first.residue_count == 12
    assert first.chain_count == 2
    assert first.total_charge == pytest.approx(-10)
    assert first.strand_mode == "duplex"
    validate_outputs(
        first.itp_path,
        first.gro_path,
        expected_name="SEQRNA",
        expected_charge=-10,
    )

    provenance = first.itp_path.read_text()
    assert "; Strand A (5'->3'): GCAUCG" in provenance
    assert "; Strand B (5'->3'): CGAUGC" in provenance
    assert "; Strand B paired orientation (3'->5'): CGUAGC" in provenance
    assert f"; Structure generator commit: {AMBERCLASSIC_COMMIT}" in provenance
    assert (
        '; Structure generator settings: fd_helix("arna", sequence, "rna")'
        in provenance
    )
    assert (
        f"; Intermediate PDB SHA-256: {first.intermediate_pdb_sha256}"
        in provenance
    )
    assert "ideal_arna_duplex.pdb" not in provenance
    assert "nucleic-builder-rna-sequence-" not in provenance


@pytest.mark.skipif(
    _backend_home() is None,
    reason=(
        "single-strand end-to-end test requires "
        "NUCLEIC_BUILDER_AMBERCLASSIC_HOME"
    ),
)
def test_single_strand_sequence_is_deterministic_and_private(tmp_path: Path) -> None:
    first = build_rna_from_sequence(
        "GCAUCG",
        "SEQSSRNA",
        tmp_path / "first",
        elastic_network="off",
        strand_mode="single",
    )
    second = build_rna_from_sequence(
        "gcaucg",
        "SEQSSRNA",
        tmp_path / "second",
        elastic_network="off",
        strand_mode="single",
    )

    assert first.itp_path.read_bytes() == second.itp_path.read_bytes()
    assert first.gro_path.read_bytes() == second.gro_path.read_bytes()
    assert sorted(path.name for path in first.itp_path.parent.iterdir()) == [
        "SEQSSRNA.gro",
        "SEQSSRNA.itp",
    ]
    assert first.input_mode == "sequence"
    assert first.sequence == "GCAUCG"
    assert first.complement is None
    assert first.strand_mode == "single"
    assert first.bead_count == 46
    assert first.residue_count == 6
    assert first.chain_count == 1
    assert first.total_charge == pytest.approx(-5)
    validate_outputs(
        first.itp_path,
        first.gro_path,
        expected_name="SEQSSRNA",
        expected_charge=-5,
    )

    provenance = first.itp_path.read_text()
    assert "; Input mode: sequence-derived idealized single-stranded RNA" in provenance
    assert "; Strand mode: single" in provenance
    assert "; Strand (5'->3'): GCAUCG" in provenance
    assert "; Strand B" not in provenance
    assert "; Structure prediction: none" in provenance
    assert "retain generated strand A only" in provenance
    assert "ideal_arna_single_strand.pdb" not in provenance
    assert "nucleic-builder-rna-sequence-" not in provenance


@pytest.mark.skipif(
    _backend_home() is None,
    reason=(
        "475-nt single-strand test requires NUCLEIC_BUILDER_AMBERCLASSIC_HOME"
    ),
)
def test_bnt162b2_475_single_mode_keeps_exactly_the_entered_chain(
    tmp_path: Path,
) -> None:
    result = build_rna_from_sequence(
        BNT162B2_FRAGMENT_475,
        "BNT162B2_475",
        tmp_path,
        elastic_network="off",
        strand_mode="single",
    )

    residues: dict[int, str] = {}
    for atom in parse_itp(result.itp_path).atoms:
        residues.setdefault(atom.residue_number, atom.residue_name)
    observed_sequence = "".join(residues.values())

    assert len(BNT162B2_FRAGMENT_475) == 475
    assert result.sequence == BNT162B2_FRAGMENT_475
    assert observed_sequence == BNT162B2_FRAGMENT_475
    assert result.complement is None
    assert result.strand_mode == "single"
    assert result.chain_count == 1
    assert result.residue_count == 475
    assert result.bead_count == 3651
    assert result.total_charge == pytest.approx(-474)
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "BNT162B2_475.gro",
        "BNT162B2_475.itp",
    ]
    provenance = result.itp_path.read_text()
    assert f"; Strand (5'->3'): {BNT162B2_FRAGMENT_475}" in provenance
    assert "; Strand B" not in provenance
