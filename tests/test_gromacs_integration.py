from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from nucleic_builder import (
    build_dna,
    build_dna_from_sequence,
    build_rna,
    build_rna_from_sequence,
)
from nucleic_builder.dna_builder import DNA_UPSTREAM_PARAMETER_DIR

from .conftest import (
    HAIRPIN,
    INDEPENDENT_RNA,
    ROOT,
    SSRNA,
    UPSTREAM_DSDNA,
    UPSTREAM_DSRNA,
)


GMX_COMMAND = os.environ.get("GMX", "gmx")
GMX = shutil.which(GMX_COMMAND)
pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        GMX is None,
        reason=(
            f"GROMACS integration test skipped: executable {GMX_COMMAND!r} was not found; "
            "set GMX to a working gmx path"
        ),
    ),
]

DATA = ROOT / "tests" / "data"
FORCEFIELD = DATA / "forcefield"


def run_gmx(
    workdir: Path,
    *arguments: str,
    stdin: str | None = None,
    timeout: int = 240,
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        [str(GMX), *arguments],
        cwd=workdir,
        input=stdin,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
        env={**os.environ, "GMX_MAXBACKUP": "-1"},
    )
    if completed.returncode != 0:
        pytest.fail(
            f"GROMACS command failed ({completed.returncode}): "
            f"{' '.join([str(GMX), *arguments])}\n"
            f"stdout:\n{completed.stdout[-8000:]}\n"
            f"stderr:\n{completed.stderr[-8000:]}"
        )
    return completed


def assert_grompp_has_no_warnings(completed: subprocess.CompletedProcess[str]) -> None:
    combined = completed.stdout + completed.stderr
    assert re.search(r"\bWARNING(?:\s+\d+)?\b", combined, re.IGNORECASE) is None, combined[-8000:]


def assert_mdrun_is_clean(
    completed: subprocess.CompletedProcess[str], log_text: str
) -> None:
    """Require a numerically clean stage, including generic warnings."""

    evidence = completed.stdout + completed.stderr + log_text
    assert re.search(r"\bWARNING(?:\s+\d+)?\b", evidence, re.IGNORECASE) is None, evidence[-8000:]
    assert re.search(r"LINCS\s+(?:WARNING|ERROR)", evidence, re.IGNORECASE) is None
    assert "fatal error" not in evidence.lower()
    assert re.search(r"\bnan\b", evidence, re.IGNORECASE) is None


def molecule_counts(topology: Path) -> dict[str, int]:
    section = None
    counts: dict[str, int] = {}
    for raw_line in topology.read_text().splitlines():
        line = raw_line.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip().lower()
            continue
        if section == "molecules" and not line.startswith("#"):
            name, count = line.split()[:2]
            counts[name] = int(count)
    return counts


@pytest.mark.parametrize(
    (
        "polymer",
        "input_pdb",
        "sequence",
        "name",
        "expected_charge",
        "expected_beads",
        "elastic_network",
    ),
    [
        ("rna", UPSTREAM_DSRNA, None, "RNA", -42, 339, "legacy"),
        ("rna", INDEPENDENT_RNA, None, "RNA1", -26, 208, "off"),
        ("rna", SSRNA, None, "SSRNA6", -5, 48, "off"),
        ("rna", HAIRPIN, None, "HAIRPIN26", -25, 203, "intrachain"),
        ("rna", None, "GCAUCG", "SEQRNA", -10, 92, "off"),
        ("dna", UPSTREAM_DSDNA, None, "DNA", -22, 190, "legacy"),
        ("dna", None, "GCATCG", "SEQDNA", -10, 94, "off"),
    ],
    ids=[
        "upstream-dsRNA-legacy",
        "independent-1RNA-off",
        "experimental-ssRNA-off",
        "experimental-hairpin-intrachain",
        "sequence-ideal-dsRNA-off",
        "experimental-upstream-dsDNA-legacy",
        "experimental-sequence-ideal-dsDNA-off",
    ],
)
def test_solvated_neutralized_em_restrained_and_unrestrained_md(
    tmp_path: Path,
    polymer: str,
    input_pdb: Path | None,
    sequence: str | None,
    name: str,
    expected_charge: int,
    expected_beads: int,
    elastic_network: str,
) -> None:
    molecule_dir = tmp_path / "molecule"
    if sequence is not None and polymer == "rna":
        if not os.environ.get("NUCLEIC_BUILDER_AMBERCLASSIC_HOME"):
            pytest.skip(
                "sequence GROMACS case requires NUCLEIC_BUILDER_AMBERCLASSIC_HOME"
            )
        result = build_rna_from_sequence(
            sequence,
            name,
            molecule_dir,
            elastic_network=elastic_network,
        )
    elif sequence is not None:
        if not os.environ.get("NUCLEIC_BUILDER_AMBERCLASSIC_HOME"):
            pytest.skip(
                "DNA sequence GROMACS case requires NUCLEIC_BUILDER_AMBERCLASSIC_HOME"
            )
        result = build_dna_from_sequence(
            sequence,
            name,
            molecule_dir,
            elastic_network=elastic_network,
        )
    else:
        assert input_pdb is not None
        if polymer == "rna":
            result = build_rna(
                input_pdb,
                name,
                molecule_dir,
                elastic_network=elastic_network,
            )
        else:
            result = build_dna(
                input_pdb,
                name,
                molecule_dir,
                elastic_network=elastic_network,
            )
    assert result.total_charge == expected_charge
    assert result.bead_count == expected_beads
    assert result.elastic_network == elastic_network

    topol_dir = tmp_path / "topol"
    topol_dir.mkdir()
    for filename in (
        "martini_v3.0.0.itp",
        "martini_v3.0.0_ions_v1.itp",
        "martini_v3.0.0_solvents_v1.itp",
    ):
        shutil.copy2(FORCEFIELD / filename, topol_dir / filename)
    if polymer == "rna":
        nucleic_ff_name = "martini_v3.0.0_rna.itp"
        nucleic_ff_source = FORCEFIELD / nucleic_ff_name
    else:
        nucleic_ff_name = "martini_alpha_dna.itp"
        nucleic_ff_source = DNA_UPSTREAM_PARAMETER_DIR / nucleic_ff_name
    shutil.copy2(nucleic_ff_source, topol_dir / nucleic_ff_name)
    shutil.copy2(result.itp_path, topol_dir / result.itp_path.name)
    shutil.copy2(DATA / "water.gro", tmp_path / "water.gro")

    topology = tmp_path / "system.top"
    topology.write_text(
        "\n".join(
            [
                '#include "topol/martini_v3.0.0.itp"',
                f'#include "topol/{nucleic_ff_name}"',
                '#include "topol/martini_v3.0.0_ions_v1.itp"',
                '#include "topol/martini_v3.0.0_solvents_v1.itp"',
                f'#include "topol/{name}.itp"',
                "",
                "[ system ]",
                f"Martini 3 {name} validation system",
                "",
                "[ molecules ]",
                f"{name}  1",
                "",
            ]
        ),
        encoding="utf-8",
    )

    run_gmx(
        tmp_path,
        "editconf",
        "-f",
        str(result.gro_path),
        "-o",
        "boxed.gro",
        "-bt",
        "dodecahedron",
        "-d",
        "1.2",
    )
    run_gmx(
        tmp_path,
        "solvate",
        "-cp",
        "boxed.gro",
        "-cs",
        "water.gro",
        "-radius",
        "0.21",
        "-o",
        "solvated.gro",
        "-p",
        "system.top",
    )
    ion_prep = run_gmx(
        tmp_path,
        "grompp",
        "-f",
        str(DATA / "em.mdp"),
        "-c",
        "solvated.gro",
        "-r",
        "solvated.gro",
        "-p",
        "system.top",
        "-o",
        "ions.tpr",
    )
    assert_grompp_has_no_warnings(ion_prep)
    run_gmx(
        tmp_path,
        "genion",
        "-s",
        "ions.tpr",
        "-o",
        "ionized.gro",
        "-p",
        "system.top",
        "-pname",
        "NA",
        "-nname",
        "CL",
        "-neutral",
        "-conc",
        "0.15",
        "-seed",
        "20260727",
        stdin="W\n",
    )
    counts = molecule_counts(topology)
    assert counts["NA"] - counts["CL"] == -expected_charge

    em_prep = run_gmx(
        tmp_path,
        "grompp",
        "-f",
        str(DATA / "em.mdp"),
        "-c",
        "ionized.gro",
        "-r",
        "ionized.gro",
        "-p",
        "system.top",
        "-o",
        "em.tpr",
    )
    assert_grompp_has_no_warnings(em_prep)
    assert "System has non-zero total charge" not in (em_prep.stdout + em_prep.stderr)
    em_run = run_gmx(
        tmp_path,
        "mdrun",
        "-deffnm",
        "em",
        "-ntmpi",
        "1",
        "-ntomp",
        "1",
        "-pin",
        "off",
    )
    assert (tmp_path / "em.gro").is_file()
    em_log = (tmp_path / "em.log").read_text()
    assert "converged to Fmax < 100" in em_log
    assert_mdrun_is_clean(em_run, em_log)

    restrained_prep = run_gmx(
        tmp_path,
        "grompp",
        "-f",
        str(DATA / "md.mdp"),
        "-c",
        "em.gro",
        "-r",
        "em.gro",
        "-p",
        "system.top",
        "-o",
        "restrained.tpr",
        "-pp",
        "restrained-processed.top",
    )
    assert_grompp_has_no_warnings(restrained_prep)
    restrained_topology = (tmp_path / "restrained-processed.top").read_text()
    assert "[ position_restraints ]" in restrained_topology
    restrained_run = run_gmx(
        tmp_path,
        "mdrun",
        "-deffnm",
        "restrained",
        "-ntmpi",
        "1",
        "-ntomp",
        "1",
        "-pin",
        "off",
    )
    assert (tmp_path / "restrained.gro").is_file()
    restrained_log = (tmp_path / "restrained.log").read_text()
    assert "Step           Time" in restrained_log
    assert_mdrun_is_clean(restrained_run, restrained_log)
    restrained_check = run_gmx(tmp_path, "check", "-f", "restrained.xtc")
    assert_grompp_has_no_warnings(restrained_check)
    restrained_check_text = restrained_check.stdout + restrained_check.stderr
    assert "Last frame" in restrained_check_text
    assert re.search(r"Last frame\s+10\s+time\s+15(?:\.0+)?", restrained_check_text)

    unrestrained_prep = run_gmx(
        tmp_path,
        "grompp",
        "-f",
        str(DATA / "md-unrestrained.mdp"),
        "-c",
        "restrained.gro",
        "-t",
        "restrained.cpt",
        "-p",
        "system.top",
        "-o",
        "unrestrained.tpr",
        "-pp",
        "unrestrained-processed.top",
    )
    assert_grompp_has_no_warnings(unrestrained_prep)
    unrestrained_topology = (tmp_path / "unrestrained-processed.top").read_text()
    assert "[ position_restraints ]" not in unrestrained_topology
    unrestrained_run = run_gmx(
        tmp_path,
        "mdrun",
        "-deffnm",
        "unrestrained",
        "-ntmpi",
        "1",
        "-ntomp",
        "1",
        "-pin",
        "off",
    )
    assert (tmp_path / "unrestrained.gro").is_file()
    unrestrained_log = (tmp_path / "unrestrained.log").read_text()
    assert "Step           Time" in unrestrained_log
    assert_mdrun_is_clean(unrestrained_run, unrestrained_log)
    unrestrained_check = run_gmx(tmp_path, "check", "-f", "unrestrained.xtc")
    assert_grompp_has_no_warnings(unrestrained_check)
    unrestrained_check_text = unrestrained_check.stdout + unrestrained_check.stderr
    assert "Last frame" in unrestrained_check_text
    assert re.search(
        r"Last frame\s+10\s+time\s+30(?:\.0+)?", unrestrained_check_text
    )
