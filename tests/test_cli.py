from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from .conftest import UPSTREAM_DSDNA, UPSTREAM_DSRNA


AMBERCLASSIC_HOME = os.environ.get("NUCLEIC_BUILDER_AMBERCLASSIC_HOME")


def test_module_cli_writes_only_itp_and_gro(tmp_path: Path) -> None:
    env = os.environ.copy()
    root = Path(__file__).resolve().parents[1]
    env["PYTHONPATH"] = str(root)
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "nucleic_builder",
            "--input",
            str(UPSTREAM_DSRNA),
            "--name",
            "CLI_RNA",
            "--output-dir",
            str(tmp_path),
            "--elastic-network",
            "legacy",
        ],
        cwd=root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "Validated 339 beads" in completed.stdout
    assert "Elastic network: legacy (366 bonds; 90 cross-chain)" in completed.stdout
    assert "INFO:" not in completed.stderr
    assert "Starting RNA Martinization" not in completed.stderr
    assert sorted(path.name for path in tmp_path.iterdir()) == ["CLI_RNA.gro", "CLI_RNA.itp"]


def test_cli_requires_explicit_elastic_network_policy(tmp_path: Path) -> None:
    env = os.environ.copy()
    root = Path(__file__).resolve().parents[1]
    env["PYTHONPATH"] = str(root)
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "nucleic_builder",
            "--input",
            str(UPSTREAM_DSRNA),
            "--name",
            "MissingPolicy",
            "--output-dir",
            str(tmp_path),
        ],
        cwd=root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 2
    assert "the following arguments are required: --elastic-network" in completed.stderr
    assert list(tmp_path.iterdir()) == []


def test_cli_verbose_enables_upstream_progress_messages(tmp_path: Path) -> None:
    env = os.environ.copy()
    root = Path(__file__).resolve().parents[1]
    env["PYTHONPATH"] = str(root)
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "nucleic_builder",
            "--input",
            str(UPSTREAM_DSRNA),
            "--name",
            "VerboseRNA",
            "--output-dir",
            str(tmp_path),
            "--elastic-network",
            "off",
            "--verbose",
        ],
        cwd=root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "INFO: === Starting RNA AA->CG Conversion ===" in completed.stderr
    assert "Elastic network: off (0 bonds; 0 cross-chain)" in completed.stdout


def test_cli_requires_exactly_one_input_mode(tmp_path: Path) -> None:
    env = os.environ.copy()
    root = Path(__file__).resolve().parents[1]
    env["PYTHONPATH"] = str(root)
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "nucleic_builder",
            "--input",
            str(UPSTREAM_DSRNA),
            "--sequence",
            "GCAUCG",
            "--name",
            "Both",
            "--output-dir",
            str(tmp_path),
            "--elastic-network",
            "off",
        ],
        cwd=root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 2
    assert "not allowed with argument" in completed.stderr
    assert list(tmp_path.iterdir()) == []


def test_cli_rejects_dna_sequence_before_backend_lookup(tmp_path: Path) -> None:
    env = os.environ.copy()
    env.pop("NUCLEIC_BUILDER_AMBERCLASSIC_HOME", None)
    root = Path(__file__).resolve().parents[1]
    env["PYTHONPATH"] = str(root)
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "nucleic_builder",
            "--sequence",
            "ACGT",
            "--name",
            "DNA",
            "--output-dir",
            str(tmp_path),
            "--elastic-network",
            "off",
        ],
        cwd=root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 2
    assert "Phase 2A accepts only A, C, G, and U" in completed.stderr
    assert not tmp_path.exists() or list(tmp_path.iterdir()) == []


def test_dna_pdb_cli_requires_explicit_selection_and_is_prominently_labeled(
    tmp_path: Path,
) -> None:
    env = os.environ.copy()
    root = Path(__file__).resolve().parents[1]
    env["PYTHONPATH"] = str(root)
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "nucleic_builder",
            "--polymer",
            "dna",
            "--input",
            str(UPSTREAM_DSDNA),
            "--name",
            "CLI_DNA",
            "--output-dir",
            str(tmp_path),
            "--elastic-network",
            "off",
        ],
        cwd=root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "Validated 190 beads, 24 residues, 2 chain(s), charge -22 e" in completed.stdout
    assert "WARNING: EXPERIMENTAL / UNPUBLISHED DNA-ALPHA MODEL" in completed.stderr
    assert "INFO:" not in completed.stderr
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "CLI_DNA.gro",
        "CLI_DNA.itp",
    ]


@pytest.mark.skipif(
    not AMBERCLASSIC_HOME,
    reason="sequence CLI test requires NUCLEIC_BUILDER_AMBERCLASSIC_HOME",
)
def test_sequence_cli_writes_only_itp_and_gro(tmp_path: Path) -> None:
    env = os.environ.copy()
    root = Path(__file__).resolve().parents[1]
    env["PYTHONPATH"] = str(root)
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "nucleic_builder",
            "--sequence",
            "GCAUCG",
            "--name",
            "SEQRNA",
            "--output-dir",
            str(tmp_path),
            "--elastic-network",
            "off",
        ],
        cwd=root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "Validated 92 beads, 12 residues, 2 chain(s), charge -10 e" in completed.stdout
    assert "Running:" not in completed.stdout
    assert "Running:" not in completed.stderr
    assert "INFO:" not in completed.stderr
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "SEQRNA.gro",
        "SEQRNA.itp",
    ]


@pytest.mark.skipif(
    not AMBERCLASSIC_HOME,
    reason="DNA sequence CLI test requires NUCLEIC_BUILDER_AMBERCLASSIC_HOME",
)
def test_dna_sequence_cli_is_explicit_labeled_and_writes_only_pair(
    tmp_path: Path,
) -> None:
    env = os.environ.copy()
    root = Path(__file__).resolve().parents[1]
    env["PYTHONPATH"] = str(root)
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "nucleic_builder",
            "--polymer",
            "dna",
            "--sequence",
            "GCATCG",
            "--name",
            "SEQDNA",
            "--output-dir",
            str(tmp_path),
            "--elastic-network",
            "off",
        ],
        cwd=root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "Validated 94 beads, 12 residues, 2 chain(s), charge -10 e" in completed.stdout
    assert "WARNING: EXPERIMENTAL / UNPUBLISHED DNA-ALPHA MODEL" in completed.stderr
    assert "Running:" not in completed.stdout + completed.stderr
    assert "INFO:" not in completed.stderr
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "SEQDNA.gro",
        "SEQDNA.itp",
    ]
