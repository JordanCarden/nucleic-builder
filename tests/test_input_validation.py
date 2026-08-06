from __future__ import annotations

from pathlib import Path

import pytest

from nucleic_builder.errors import InputValidationError
from nucleic_builder.rna_pdb_input import prepare_input_pdb


def atom_line(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain: str,
    residue_number: int,
    *,
    record: str = "ATOM",
    altloc: str = "",
) -> str:
    return (
        f"{record:<6s}{serial:5d} {atom_name:>4s}{altloc:1s}{residue_name:>3s} "
        f"{chain:1s}{residue_number:4d}    {serial:8.3f}{serial + 1:8.3f}"
        f"{serial + 2:8.3f}  1.00  0.00           C"
    )


BASE_HEAVY_ATOMS = {
    "A": {"N9", "C8", "N3", "C4", "N1", "C2", "N6", "C6", "N7", "C5"},
    "C": {"N1", "C5", "C6", "C2", "O2", "N3", "N4", "C4"},
    "G": {"C8", "N9", "C4", "N3", "C2", "N2", "N1", "C6", "O6", "C5", "N7"},
    "U": {"N1", "C5", "C6", "C2", "O2", "N3", "C4", "O4"},
}
SUGAR_HEAVY_ATOMS = {"C5'", "C4'", "O4'", "C3'", "C1'", "C2'", "O2'"}


def complete_residue_lines(
    first_serial: int,
    residue_name: str,
    canonical_name: str,
    chain: str,
    residue_number: int,
    *,
    phosphate: bool = False,
    o3: bool = False,
) -> list[str]:
    names = set(SUGAR_HEAVY_ATOMS | BASE_HEAVY_ATOMS[canonical_name])
    if phosphate:
        names.update({"P", "OP1", "OP2", "O5'"})
    if o3:
        names.add("O3'")
    return [
        atom_line(first_serial + offset, name, residue_name, chain, residue_number)
        for offset, name in enumerate(sorted(names))
    ]


def write_pdb(path: Path, lines: list[str]) -> Path:
    path.write_text("\n".join(lines + ["END"]) + "\n", encoding="utf-8")
    return path


def test_terminal_variants_are_position_checked_and_normalized(tmp_path: Path) -> None:
    source = write_pdb(
        tmp_path / "terminal.pdb",
        complete_residue_lines(1, "RA5", "A", "A", 7, o3=True)
        + complete_residue_lines(30, "U3", "U", "A", 8, phosphate=True, o3=True),
    )
    normalized = tmp_path / "normalized.pdb"
    summary = prepare_input_pdb(source, normalized)

    assert summary.sequence == "AU"
    assert [res.output_number for res in summary.residues] == [1, 2]
    records = [line for line in normalized.read_text().splitlines() if line.startswith("ATOM")]
    assert list(dict.fromkeys(line[17:20].strip() for line in records)) == ["A", "U"]
    assert list(dict.fromkeys(int(line[22:26]) for line in records)) == [1, 2]


def test_terminal_variant_in_the_middle_is_rejected(tmp_path: Path) -> None:
    source = write_pdb(
        tmp_path / "bad-terminal.pdb",
        [
            atom_line(1, "C1'", "A", "A", 1),
            atom_line(2, "C1'", "RC5", "A", 2),
            atom_line(3, "C1'", "G", "A", 3),
        ],
    )
    with pytest.raises(InputValidationError, match="5' terminal name RC5 is not first"):
        prepare_input_pdb(source, tmp_path / "unused.pdb")


@pytest.mark.parametrize("residue_name", ["PSU", "5MC", "1MG", "2MA"])
def test_modified_nucleotides_are_never_canonicalized(
    tmp_path: Path, residue_name: str
) -> None:
    source = write_pdb(
        tmp_path / f"{residue_name}.pdb",
        [atom_line(1, "C1'", residue_name, "A", 1, record="HETATM")],
    )
    with pytest.raises(InputValidationError, match="Modified nucleotides.*will not be mapped"):
        prepare_input_pdb(source, tmp_path / "unused.pdb")


def test_unknown_residue_is_rejected_with_location(tmp_path: Path) -> None:
    source = write_pdb(
        tmp_path / "unknown.pdb",
        [atom_line(1, "C1'", "DA", "B", 4)],
    )
    with pytest.raises(InputValidationError, match=r"DA at chain B residue 4"):
        prepare_input_pdb(source, tmp_path / "unused.pdb")


def test_crystallographic_water_is_explicitly_removed(tmp_path: Path) -> None:
    source = write_pdb(
        tmp_path / "water.pdb",
        complete_residue_lines(1, "A", "A", "A", 1)
        + [atom_line(30, "O", "HOH", "W", 2, record="HETATM")],
    )
    normalized = tmp_path / "normalized.pdb"
    summary = prepare_input_pdb(source, normalized)
    assert summary.ignored_water_atoms == 1
    assert "HOH" not in normalized.read_text()


def test_multiple_models_are_rejected(tmp_path: Path) -> None:
    source = write_pdb(
        tmp_path / "models.pdb",
        [
            "MODEL        1",
            atom_line(1, "C1'", "A", "A", 1),
            "ENDMDL",
            "MODEL        2",
            atom_line(2, "C1'", "A", "A", 1),
            "ENDMDL",
        ],
    )
    with pytest.raises(InputValidationError, match="exactly one PDB model"):
        prepare_input_pdb(source, tmp_path / "unused.pdb")


def test_alternate_locations_are_rejected(tmp_path: Path) -> None:
    source = write_pdb(
        tmp_path / "altloc.pdb",
        [atom_line(1, "C1'", "A", "A", 1, altloc="A")],
    )
    with pytest.raises(InputValidationError, match="Alternate location 'A'"):
        prepare_input_pdb(source, tmp_path / "unused.pdb")


def test_interleaved_chains_are_rejected(tmp_path: Path) -> None:
    source = write_pdb(
        tmp_path / "interleaved.pdb",
        [
            atom_line(1, "C1'", "A", "A", 1),
            atom_line(2, "C1'", "C", "B", 1),
            atom_line(3, "C1'", "G", "A", 2),
        ],
    )
    with pytest.raises(InputValidationError, match="multiple non-contiguous blocks"):
        prepare_input_pdb(source, tmp_path / "unused.pdb")


def test_ter_preserves_blank_reused_chain_ids_as_distinct_segments(tmp_path: Path) -> None:
    source = write_pdb(
        tmp_path / "blank-chains.pdb",
        complete_residue_lines(1, "A", "A", "", 1)
        + ["TER"]
        + complete_residue_lines(30, "U", "U", "", 1),
    )
    normalized = tmp_path / "normalized.pdb"
    summary = prepare_input_pdb(source, normalized)

    assert summary.chain_count == 2
    assert [res.chain_id for res in summary.residues] == ["", ""]
    assert [res.normalized_chain_id for res in summary.residues] == ["A", "B"]
    atom_chains = [
        line[21:22]
        for line in normalized.read_text().splitlines()
        if line.startswith("ATOM")
    ]
    assert list(dict.fromkeys(atom_chains)) == ["A", "B"]
