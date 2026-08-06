from __future__ import annotations

from pathlib import Path

import pytest

from nucleic_builder import BuildResult, build_dna, build_rna


ROOT = Path(__file__).resolve().parents[1]
UPSTREAM_DSRNA = (
    ROOT
    / "nucleic_builder"
    / "_vendor"
    / "martini_3_dna_rna"
    / "tests"
    / "dsRNA.pdb"
)
UPSTREAM_DSDNA = (
    ROOT
    / "nucleic_builder"
    / "_vendor"
    / "martini_3_dna_rna"
    / "tests"
    / "dsDNA.pdb"
)
INDEPENDENT_RNA = ROOT / "examples" / "1RNA.pdb"
SSRNA = ROOT / "examples" / "ssrna" / "3G9Y-chain-C-alt-A.pdb"
HAIRPIN = ROOT / "examples" / "hairpin" / "6YMC-chain-A.pdb"


@pytest.fixture(scope="session")
def dsrna_build(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    return build_rna(
        UPSTREAM_DSRNA,
        "RNA",
        tmp_path_factory.mktemp("dsrna"),
        elastic_network="legacy",
    )


@pytest.fixture(scope="session")
def independent_build(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    return build_rna(
        INDEPENDENT_RNA,
        "RNA1",
        tmp_path_factory.mktemp("independent"),
        elastic_network="off",
    )


@pytest.fixture(scope="session")
def ssrna_build(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    return build_rna(
        SSRNA,
        "SSRNA6",
        tmp_path_factory.mktemp("ssrna"),
        elastic_network="off",
    )


@pytest.fixture(scope="session")
def hairpin_build(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    return build_rna(
        HAIRPIN,
        "HAIRPIN26",
        tmp_path_factory.mktemp("hairpin"),
        elastic_network="intrachain",
    )


@pytest.fixture(scope="session")
def dsdna_build(tmp_path_factory: pytest.TempPathFactory) -> BuildResult:
    return build_dna(
        UPSTREAM_DSDNA,
        "DNA",
        tmp_path_factory.mktemp("dsdna"),
        elastic_network="off",
    )
