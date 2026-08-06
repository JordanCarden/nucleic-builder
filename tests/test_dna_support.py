from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pytest

from nucleic_builder import build_dna, build_dna_from_sequence
from nucleic_builder.core import (
    AMBERCLASSIC_COMMIT,
    AMBERCLASSIC_HOME_ENV,
    parse_gro,
    validate_outputs,
)
from nucleic_builder.dna_builder import (
    DNA_MODEL_WARNING,
    DNA_UPSTREAM_COMMIT,
    DNA_UPSTREAM_PARAMETER_DIR,
    DNA_UPSTREAM_SCRIPT,
)
from nucleic_builder.dna_pdb_input import prepare_dna_input_pdb
from nucleic_builder.dna_sequence_builder import (
    DNA_NAB_SETTINGS,
    validate_dna_duplex_sequence,
)
from nucleic_builder.errors import InputValidationError

from .conftest import ROOT, UPSTREAM_DSDNA


def _backend_home() -> Path | None:
    value = os.environ.get(AMBERCLASSIC_HOME_ENV)
    if not value:
        return None
    path = Path(value)
    return path if path.exists() else None


def test_dna_sequence_is_normalized_and_complemented_antiparallel() -> None:
    sequence = validate_dna_duplex_sequence("acgtt")
    assert sequence.strand_a_5to3 == "ACGTT"
    assert sequence.strand_b_paired_3to5 == "TGCAA"
    assert sequence.strand_b_5to3 == "AACGT"


@pytest.mark.parametrize(
    "sequence",
    ["", "ACGU", "ACGN", "AC GT", "AC-GT", ">strand\nACGT", "A5CGT"],
)
def test_dna_sequence_rejects_noncanonical_or_file_syntax(sequence: str) -> None:
    with pytest.raises(InputValidationError, match="only A, C, G, and T"):
        validate_dna_duplex_sequence(sequence)


def test_dna_pdb_validator_accepts_only_canonical_dna(tmp_path: Path) -> None:
    normalized = tmp_path / "normalized.pdb"
    summary = prepare_dna_input_pdb(UPSTREAM_DSDNA, normalized)
    assert summary.chain_count == 2
    assert summary.residue_count == 24
    assert summary.sequence == "ACGTACGTACGTACGTACGTACGT"
    names = {
        line[17:20].strip()
        for line in normalized.read_text().splitlines()
        if line.startswith("ATOM")
    }
    assert names == {"A", "C", "G", "T"}


def test_dna_pdb_validator_rejects_rna_uracil(tmp_path: Path) -> None:
    lines = UPSTREAM_DSDNA.read_text().splitlines()
    first_residue = [
        line[:17] + "  U" + line[20:]
        if line.startswith("ATOM") and line[21:22] == "A" and int(line[22:26]) == 1
        else line
        for line in lines
    ]
    source = tmp_path / "rna-in-dna.pdb"
    source.write_text("\n".join(first_residue) + "\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match=r"U at chain A residue 1"):
        prepare_dna_input_pdb(source, tmp_path / "unused.pdb")


def test_dna_pdb_validator_rejects_ribose_o2_prime(tmp_path: Path) -> None:
    lines = UPSTREAM_DSDNA.read_text().splitlines()
    rewritten: list[str] = []
    inserted = False
    for line in lines:
        rewritten.append(line)
        if (
            not inserted
            and line.startswith("ATOM")
            and line[21:22] == "A"
            and int(line[22:26]) == 1
        ):
            extra = line[:6] + f"{99999:5d}" + line[11:12] + " O2'" + line[16:]
            rewritten.append(extra)
            inserted = True
    source = tmp_path / "ribose-on-dna.pdb"
    source.write_text("\n".join(rewritten) + "\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="RNA ribose O2' atom"):
        prepare_dna_input_pdb(source, tmp_path / "unused.pdb")


def test_dna_pdb_validator_rejects_modified_nucleotide(tmp_path: Path) -> None:
    lines = UPSTREAM_DSDNA.read_text().splitlines()
    modified = [
        line[:17] + "5MC" + line[20:]
        if line.startswith("ATOM") and line[21:22] == "A" and int(line[22:26]) == 2
        else line
        for line in lines
    ]
    source = tmp_path / "modified-dna.pdb"
    source.write_text("\n".join(modified) + "\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="Modified DNA nucleotides"):
        prepare_dna_input_pdb(source, tmp_path / "unused.pdb")


def test_dna_alpha_pdb_conversion_is_isolated_and_consistent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "nucleic_builder.rna_builder._load_upstream",
        lambda: pytest.fail("DNA path attempted to load the RNA converter"),
    )
    result = build_dna(UPSTREAM_DSDNA, "DNA", tmp_path, elastic_network="off")
    assert result.polymer_type == "dna"
    assert result.bead_count == 190
    assert result.residue_count == 24
    assert result.chain_count == 2
    assert result.total_charge == pytest.approx(-22)
    assert result.input_sha256 == hashlib.sha256(UPSTREAM_DSDNA.read_bytes()).hexdigest()
    assert sorted(path.name for path in result.itp_path.parent.iterdir()) == [
        "DNA.gro",
        "DNA.itp",
    ]
    validate_outputs(
        result.itp_path,
        result.gro_path,
        expected_name="DNA",
        expected_charge=-22,
    )
    provenance = result.itp_path.read_text()
    assert f"; !!! {DNA_MODEL_WARNING} !!!" in provenance
    assert f"; DNA-alpha upstream commit: {DNA_UPSTREAM_COMMIT}" in provenance
    assert "; DNA-alpha converter: martinize_dna_alpha.py" in provenance
    assert "; Input SHA-256:" in provenance
    assert "normalized-dna-input.pdb" not in provenance
    assert DNA_MODEL_WARNING in parse_gro(result.gro_path).title


def test_dna_runtime_assets_are_separately_named_and_pinned() -> None:
    assert DNA_UPSTREAM_SCRIPT.name == "martinize_dna_alpha.py"
    assert DNA_UPSTREAM_PARAMETER_DIR.name == "dna_alpha_itps"
    assert (DNA_UPSTREAM_PARAMETER_DIR / "martini_alpha_dna.itp").is_file()
    assert DNA_UPSTREAM_COMMIT == "e761b7349fdf61dd485053c000dbb642f24ff9d8"
    assert DNA_UPSTREAM_SCRIPT.read_bytes() != (
        ROOT
        / "nucleic_builder"
        / "_vendor"
        / "martini_3_dna_rna"
        / "martinize_rna_v3.0.0.py"
    ).read_bytes()


@pytest.mark.skipif(
    _backend_home() is None,
    reason=(
        "DNA sequence end-to-end test skipped: set NUCLEIC_BUILDER_AMBERCLASSIC_HOME "
        "to the pinned, built AmberClassic checkout"
    ),
)
def test_dna_sequence_end_to_end_is_deterministic_private_and_labeled(
    tmp_path: Path,
) -> None:
    first = build_dna_from_sequence(
        "GCATCG", "SEQDNA", tmp_path / "first", elastic_network="off"
    )
    second = build_dna_from_sequence(
        "gcatcg", "SEQDNA", tmp_path / "second", elastic_network="off"
    )
    assert first.itp_path.read_bytes() == second.itp_path.read_bytes()
    assert first.gro_path.read_bytes() == second.gro_path.read_bytes()
    checked_in = ROOT / "examples" / "dna" / "generated"
    assert (checked_in / "SEQDNA.itp").read_bytes() == first.itp_path.read_bytes()
    assert (checked_in / "SEQDNA.gro").read_bytes() == first.gro_path.read_bytes()
    assert sorted(path.name for path in first.itp_path.parent.iterdir()) == [
        "SEQDNA.gro",
        "SEQDNA.itp",
    ]
    assert first.polymer_type == "dna"
    assert first.input_mode == "sequence"
    assert first.sequence == "GCATCG"
    assert first.complement == "CGATGC"
    assert first.intermediate_pdb_sha256 == (
        "fb960d1cdfe84ec74e22e69b3c4c08d8de852f640d641ed98f7b5d3b85e6f028"
    )
    assert first.bead_count == 94
    assert first.residue_count == 12
    assert first.chain_count == 2
    assert first.total_charge == pytest.approx(-10)

    provenance = first.itp_path.read_text()
    assert "; Strand A (5'->3'): GCATCG" in provenance
    assert "; Strand B (5'->3'): CGATGC" in provenance
    assert "; Strand B paired orientation (3'->5'): CGTAGC" in provenance
    assert f"; Structure generator commit: {AMBERCLASSIC_COMMIT}" in provenance
    assert f"; Structure generator settings: {DNA_NAB_SETTINGS}" in provenance
    assert f"; Intermediate PDB SHA-256: {first.intermediate_pdb_sha256}" in provenance
    assert "ideal_abdna_duplex.pdb" not in provenance
    assert "nucleic-builder-dna-sequence-" not in provenance
    assert DNA_MODEL_WARNING in first.gro_path.read_text().splitlines()[0]
