from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from nucleic_builder import build_rna
from nucleic_builder.core import parse_gro, parse_itp, validate_outputs
from nucleic_builder.errors import ConversionError, InputValidationError, OutputExistsError

from .conftest import HAIRPIN, INDEPENDENT_RNA, ROOT, SSRNA, UPSTREAM_DSRNA


def test_upstream_dsrna_counts_charge_and_identity(dsrna_build) -> None:
    assert dsrna_build.bead_count == 339
    assert dsrna_build.residue_count == 44
    assert dsrna_build.chain_count == 2
    assert dsrna_build.total_charge == pytest.approx(-42.0)
    assert dsrna_build.ignored_water_atoms == 0

    report = validate_outputs(
        dsrna_build.itp_path,
        dsrna_build.gro_path,
        expected_name="RNA",
        expected_charge=-42,
    )
    assert report.bead_count == 339


def test_independent_1rna_counts_charge_and_water_handling(independent_build) -> None:
    assert independent_build.bead_count == 208
    assert independent_build.residue_count == 28
    assert independent_build.chain_count == 2
    assert independent_build.total_charge == pytest.approx(-26.0)
    assert independent_build.ignored_water_atoms == 91
    assert independent_build.elastic_network == "off"
    assert independent_build.elastic_bond_count == 0
    validate_outputs(
        independent_build.itp_path,
        independent_build.gro_path,
        expected_name="RNA1",
        expected_charge=-26,
    )


def test_experimental_single_strand_counts_and_identity(ssrna_build) -> None:
    assert ssrna_build.bead_count == 48
    assert ssrna_build.residue_count == 6
    assert ssrna_build.chain_count == 1
    assert ssrna_build.total_charge == pytest.approx(-5.0)
    assert ssrna_build.elastic_network == "off"
    assert ssrna_build.elastic_bond_count == 0
    assert ssrna_build.cross_chain_elastic_bond_count == 0
    validate_outputs(
        ssrna_build.itp_path,
        ssrna_build.gro_path,
        expected_name="SSRNA6",
        expected_charge=-5,
    )


def test_experimental_hairpin_counts_and_identity(hairpin_build) -> None:
    assert hairpin_build.bead_count == 203
    assert hairpin_build.residue_count == 26
    assert hairpin_build.chain_count == 1
    assert hairpin_build.total_charge == pytest.approx(-25.0)
    assert hairpin_build.elastic_network == "intrachain"
    assert hairpin_build.elastic_bond_count == 230
    assert hairpin_build.cross_chain_elastic_bond_count == 0
    validate_outputs(
        hairpin_build.itp_path,
        hairpin_build.gro_path,
        expected_name="HAIRPIN26",
        expected_charge=-25,
    )


def test_checked_in_example_outputs_are_reproducible(independent_build) -> None:
    generated = ROOT / "examples" / "generated"
    assert (generated / "RNA1.itp").read_bytes() == independent_build.itp_path.read_bytes()
    assert (generated / "RNA1.gro").read_bytes() == independent_build.gro_path.read_bytes()
    assert sorted(path.suffix for path in generated.iterdir()) == [".gro", ".itp"]


@pytest.mark.parametrize(
    ("fixture_name", "relative_directory", "basename"),
    [
        ("ssrna_build", Path("ssrna/generated"), "SSRNA6"),
        ("hairpin_build", Path("hairpin/generated"), "HAIRPIN26"),
    ],
)
def test_new_checked_in_outputs_are_reproducible(
    request: pytest.FixtureRequest,
    fixture_name: str,
    relative_directory: Path,
    basename: str,
) -> None:
    result = request.getfixturevalue(fixture_name)
    generated = ROOT / "examples" / relative_directory
    assert (generated / f"{basename}.itp").read_bytes() == result.itp_path.read_bytes()
    assert (generated / f"{basename}.gro").read_bytes() == result.gro_path.read_bytes()
    assert sorted(path.name for path in generated.iterdir()) == [
        f"{basename}.gro",
        f"{basename}.itp",
    ]


@pytest.mark.parametrize("source", [UPSTREAM_DSRNA, INDEPENDENT_RNA, SSRNA, HAIRPIN])
def test_input_sha256_is_returned_and_recorded(source: Path, tmp_path: Path) -> None:
    expected = hashlib.sha256(source.read_bytes()).hexdigest()
    result = build_rna(source, "HashRNA", tmp_path / source.stem)
    assert result.input_sha256 == expected
    assert f"; Input SHA-256: {expected}" in result.itp_path.read_text()


def test_elastic_network_default_preserves_legacy_duplex_behavior(tmp_path: Path) -> None:
    default = build_rna(UPSTREAM_DSRNA, "RNA", tmp_path / "default")
    explicit = build_rna(
        UPSTREAM_DSRNA,
        "RNA",
        tmp_path / "explicit",
        elastic_network="legacy",
    )
    assert default.itp_path.read_bytes() == explicit.itp_path.read_bytes()
    assert default.gro_path.read_bytes() == explicit.gro_path.read_bytes()
    assert default.elastic_bond_count == 366
    assert default.cross_chain_elastic_bond_count == 90


def test_intrachain_policy_removes_only_cross_chain_elastic_bonds(tmp_path: Path) -> None:
    legacy = build_rna(UPSTREAM_DSRNA, "RNA", tmp_path / "legacy")
    intrachain = build_rna(
        UPSTREAM_DSRNA,
        "RNA",
        tmp_path / "intrachain",
        elastic_network="intrachain",
    )
    assert intrachain.elastic_bond_count == 276
    assert intrachain.cross_chain_elastic_bond_count == 0
    assert legacy.elastic_bond_count - intrachain.elastic_bond_count == 90
    assert legacy.gro_path.read_bytes() == intrachain.gro_path.read_bytes()


def test_off_policy_has_no_elastic_bonds(tmp_path: Path) -> None:
    result = build_rna(
        HAIRPIN,
        "HairpinOff",
        tmp_path,
        elastic_network="off",
    )
    assert result.elastic_bond_count == 0
    assert result.cross_chain_elastic_bond_count == 0
    assert "; Elastic network policy: off" in result.itp_path.read_text()


def test_unknown_elastic_network_policy_is_rejected_without_outputs(tmp_path: Path) -> None:
    output = tmp_path / "out"
    with pytest.raises(ConversionError, match="Unknown elastic-network policy"):
        build_rna(UPSTREAM_DSRNA, "RNA", output, elastic_network="automatic")
    assert not output.exists()


def test_output_has_only_requested_files_and_no_private_paths(tmp_path: Path) -> None:
    result = build_rna(UPSTREAM_DSRNA, "Clean", tmp_path)
    assert sorted(path.name for path in tmp_path.iterdir()) == ["Clean.gro", "Clean.itp"]
    text = result.itp_path.read_text()
    assert "normalized-input.pdb" not in text
    assert ".nucleic-builder-" not in text
    assert "; Generated on:" not in text


def test_name_is_both_filename_and_moleculetype(tmp_path: Path) -> None:
    result = build_rna(UPSTREAM_DSRNA, "Duplex_1", tmp_path)
    assert result.itp_path.name == "Duplex_1.itp"
    assert result.gro_path.name == "Duplex_1.gro"
    assert parse_itp(result.itp_path).molecule_name == "Duplex_1"


@pytest.mark.parametrize("name", ["../RNA", "RNA name", "1RNA", "RNA;bad", ""])
def test_unsafe_or_invalid_name_is_rejected(tmp_path: Path, name: str) -> None:
    with pytest.raises(ConversionError, match="Molecule name"):
        build_rna(UPSTREAM_DSRNA, name, tmp_path)


def test_existing_outputs_require_force(tmp_path: Path) -> None:
    first = build_rna(UPSTREAM_DSRNA, "RNA", tmp_path)
    original = first.itp_path.read_bytes()
    with pytest.raises(OutputExistsError, match="Refusing to overwrite"):
        build_rna(UPSTREAM_DSRNA, "RNA", tmp_path)
    assert first.itp_path.read_bytes() == original
    replacement = build_rna(UPSTREAM_DSRNA, "RNA", tmp_path, force=True)
    assert replacement.bead_count == 339


def test_chain_numbering_that_restarts_is_normalized_consistently(tmp_path: Path) -> None:
    lines = INDEPENDENT_RNA.read_text().splitlines()
    rewritten = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and line[21:22] == "B":
            old_number = int(line[22:26])
            line = line[:22] + f"{old_number - 14:4d}" + line[26:]
        rewritten.append(line)
    source = tmp_path / "restart.pdb"
    source.write_text("\n".join(rewritten) + "\n")

    result = build_rna(source, "Restart", tmp_path / "out")
    itp = parse_itp(result.itp_path)
    gro = parse_gro(result.gro_path)
    assert max(atom.residue_number for atom in itp.atoms) == 28
    assert max(atom.residue_number for atom in gro.atoms) == 28
    assert [atom.residue_number for atom in gro.atoms] == [
        atom.residue_number for atom in itp.atoms
    ]
    assert "B: input 1-14 => residues 15-28" in result.itp_path.read_text()


def test_ter_separates_strands_when_source_chain_id_is_reused(tmp_path: Path) -> None:
    rewritten = []
    for line in UPSTREAM_DSRNA.read_text().splitlines():
        if line.startswith("ATOM") and line[21:22] == "B":
            line = line[:21] + "A" + line[22:]
        rewritten.append(line)
    source = tmp_path / "reused-chain.pdb"
    source.write_text("\n".join(rewritten) + "\n")

    result = build_rna(source, "TERDuplex", tmp_path / "out")
    assert result.chain_count == 2
    assert result.bead_count == 339
    assert result.total_charge == pytest.approx(-42.0)
    provenance = result.itp_path.read_text()
    assert "A (segment 1): input 1-22 => residues 1-22" in provenance
    assert "A (segment 2): input 23-44 => residues 23-44" in provenance


def test_incomplete_mapping_fails_without_partial_outputs(tmp_path: Path) -> None:
    lines = []
    for line in UPSTREAM_DSRNA.read_text().splitlines():
        # Remove every atom contributing to BB2 in residue A:2.  The upstream
        # converter would silently omit the bead while still writing topology.
        if line.startswith("ATOM") and line[21:22] == "A" and int(line[22:26]) == 2:
            if line[12:16].strip() in {"C5'", "1H5'", "2H5'", "H5'", "H5''", "C4'", "H4'", "O4'", "C3'", "H3'"}:
                continue
        lines.append(line)
    source = tmp_path / "incomplete.pdb"
    source.write_text("\n".join(lines) + "\n")
    output = tmp_path / "out"
    with pytest.raises(InputValidationError, match="Incomplete all-atom mapping"):
        build_rna(source, "RNA", output)
    assert list(output.iterdir()) == []


def test_partial_mapping_group_is_rejected_before_conversion(tmp_path: Path) -> None:
    lines = [
        line
        for line in UPSTREAM_DSRNA.read_text().splitlines()
        if not (
            line.startswith("ATOM")
            and line[21:22] == "A"
            and int(line[22:26]) == 1
            and line[12:16].strip() == "O6"
        )
    ]
    source = tmp_path / "missing-o6.pdb"
    source.write_text("\n".join(lines) + "\n")
    output = tmp_path / "out"

    with pytest.raises(
        InputValidationError,
        match=r"G at chain A residue 1: missing required mapped heavy atom\(s\): O6",
    ):
        build_rna(source, "RNA", output)
    assert list(output.iterdir()) == []


def test_terminal_o3_required_to_complete_last_backbone_bead(tmp_path: Path) -> None:
    lines = [
        line
        for line in UPSTREAM_DSRNA.read_text().splitlines()
        if not (
            line.startswith("ATOM")
            and line[21:22] == "A"
            and int(line[22:26]) == 22
            and line[12:16].strip() == "O3'"
        )
    ]
    source = tmp_path / "missing-terminal-o3.pdb"
    source.write_text("\n".join(lines) + "\n")
    output = tmp_path / "out"

    with pytest.raises(
        InputValidationError,
        match=r"C at chain A residue 22: missing required mapped heavy atom\(s\): O3'",
    ):
        build_rna(source, "RNA", output)
    assert list(output.iterdir()) == []


def test_modified_input_creates_no_outputs(tmp_path: Path) -> None:
    lines = UPSTREAM_DSRNA.read_text().splitlines()
    rewritten = [
        line[:17] + "PSU" + line[20:]
        if line.startswith("ATOM") and line[21:22] == "A" and int(line[22:26]) == 1
        else line
        for line in lines
    ]
    source = tmp_path / "modified.pdb"
    source.write_text("\n".join(rewritten) + "\n")
    output = tmp_path / "out"
    with pytest.raises(InputValidationError, match="Modified nucleotides"):
        build_rna(source, "RNA", output)
    assert list(output.iterdir()) == []


def test_sequence_only_input_remains_rejected(tmp_path: Path) -> None:
    source = tmp_path / "sequence.txt"
    source.write_text("ACGUACGU\n", encoding="utf-8")
    output = tmp_path / "out"
    with pytest.raises(InputValidationError, match="no canonical RNA atom records"):
        build_rna(source, "RNA", output)
    assert list(output.iterdir()) == []
