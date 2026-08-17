from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pytest

from nucleic_builder import build_dna, build_rna, build_rna_from_sequence
from nucleic_builder.core import AMBERCLASSIC_HOME_ENV, parse_itp, validate_outputs
from nucleic_builder.martini2_builder import (
    MARTINI2_CONVERTER,
    MARTINI2_FORCE_FIELD,
    MARTINI2_IONS,
    MARTINI2_POLARIZABLE_FORCE_FIELD,
)

from .conftest import SSRNA, UPSTREAM_DSDNA, UPSTREAM_DSRNA


EXPECTED_CHECKSUMS = {
    MARTINI2_CONVERTER: "ee858476b4e09e0f13d0131ed9b2f617792ca87142af8364daa453913eb8e9fd",
    MARTINI2_FORCE_FIELD: "cc7c200dff400e97311213b93127697c6f8c21edb2350926072f0194eb90efe6",
    MARTINI2_POLARIZABLE_FORCE_FIELD: "b8dea4ffbef3a439db0baa465825528db919b026b1a821b374f1a6a605912ff0",
    MARTINI2_IONS: "c5b9b5b9541aa6d77b5b41a4b19dee62c1b8631c73e79cf76e1b281a144b4b4e",
}


def test_packaged_martini2_assets_match_the_pinned_release() -> None:
    for path, expected in EXPECTED_CHECKSUMS.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected


@pytest.mark.parametrize(
    (
        "polymer",
        "input_pdb",
        "name",
        "elastic_network",
        "expected_beads",
        "expected_residues",
        "expected_chains",
        "expected_charge",
        "expected_topology_type",
    ),
    [
        ("rna", SSRNA, "M2SSRNA", "off", 40, 6, 1, -5, "ss"),
        ("rna", UPSTREAM_DSRNA, "M2DSRNA", "legacy", 284, 44, 2, -42, "ds-stiff"),
        ("dna", UPSTREAM_DSDNA, "M2DSDNA", "legacy", 154, 24, 2, -22, "ds-stiff"),
    ],
)
def test_martini2_pdb_backends_produce_valid_versioned_outputs(
    tmp_path: Path,
    polymer: str,
    input_pdb: Path,
    name: str,
    elastic_network: str,
    expected_beads: int,
    expected_residues: int,
    expected_chains: int,
    expected_charge: int,
    expected_topology_type: str,
) -> None:
    builder = build_rna if polymer == "rna" else build_dna
    result = builder(
        input_pdb,
        name,
        tmp_path,
        elastic_network=elastic_network,
        martini_version=2,
    )

    assert result.martini_version == 2
    assert result.polymer_type == polymer
    assert result.bead_count == expected_beads
    assert result.residue_count == expected_residues
    assert result.chain_count == expected_chains
    assert result.total_charge == pytest.approx(expected_charge)
    assert result.backend_topology_type == expected_topology_type
    assert result.force_field_files == (MARTINI2_FORCE_FIELD, MARTINI2_IONS)
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        f"{name}.gro",
        f"{name}.itp",
    ]
    validate_outputs(
        result.itp_path,
        result.gro_path,
        expected_name=name,
        expected_charge=expected_charge,
    )
    provenance = result.itp_path.read_text(encoding="utf-8")
    assert "; Martini version: 2" in provenance
    assert "; Force-field parameters: martini_v2.1-dna.itp" in provenance


@pytest.mark.sequence_backend
@pytest.mark.skipif(
    not os.environ.get(AMBERCLASSIC_HOME_ENV),
    reason=f"Martini 2 sequence tests require {AMBERCLASSIC_HOME_ENV}",
)
@pytest.mark.parametrize(
    ("strand_mode", "expected_beads", "expected_residues", "expected_chains", "expected_charge"),
    [
        ("single", 38, 6, 1, -5),
        ("duplex", 76, 12, 2, -10),
    ],
)
def test_martini2_rna_sequence_modes_preserve_the_requested_strand_count(
    tmp_path: Path,
    strand_mode: str,
    expected_beads: int,
    expected_residues: int,
    expected_chains: int,
    expected_charge: int,
) -> None:
    result = build_rna_from_sequence(
        "GCAUCG",
        f"M2{strand_mode.upper()}",
        tmp_path,
        elastic_network="off",
        strand_mode=strand_mode,
        martini_version=2,
    )

    assert result.martini_version == 2
    assert result.strand_mode == strand_mode
    assert result.bead_count == expected_beads
    assert result.residue_count == expected_residues
    assert result.chain_count == expected_chains
    assert result.total_charge == pytest.approx(expected_charge)
    assert result.complement == (None if strand_mode == "single" else "CGAUGC")
    assert len(parse_itp(result.itp_path).atoms) == expected_beads
