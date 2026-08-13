"""Shared conversion, output, and validation mechanics."""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
from types import ModuleType
from typing import Any

from .errors import ConversionError, OutputValidationError


VENDOR_ROOT = Path(__file__).parent / "_vendor" / "martini_3_dna_rna"
ELASTIC_NETWORK_POLICIES = frozenset({"legacy", "intrachain", "off"})
MARTINI_VERSIONS = frozenset({2, 3})
AMBERCLASSIC_REPOSITORY = "https://github.com/Amber-MD/AmberClassic"
AMBERCLASSIC_TAG = "v2.0"
AMBERCLASSIC_COMMIT = "bdb3e0dee5b90f2be2950e26cfad1ae5a7440cae"
AMBERCLASSIC_HOME_ENV = "NUCLEIC_BUILDER_AMBERCLASSIC_HOME"
_NAME_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_.+-]{0,63}$")


@dataclass(frozen=True)
class BuildResult:
    """Paths and validated summary for one generated nucleic-acid molecule."""

    itp_path: Path
    gro_path: Path
    bead_count: int
    residue_count: int
    chain_count: int
    total_charge: float
    ignored_water_atoms: int
    input_sha256: str
    elastic_network: str
    elastic_bond_count: int
    cross_chain_elastic_bond_count: int
    input_mode: str = "pdb"
    sequence: str | None = None
    complement: str | None = None
    intermediate_pdb_sha256: str | None = None
    polymer_type: str = "rna"
    strand_mode: str | None = None
    martini_version: int = 3
    force_field_files: tuple[Path, ...] = ()
    backend_topology_type: str | None = None


@dataclass(frozen=True)
class GeneratedStructureProvenance:
    """Provenance for a private sequence-derived all-atom intermediate."""

    sequence_5to3: str
    complement_5to3: str | None
    paired_complement_3to5: str | None
    generator_name: str
    generator_version: str
    generator_repository: str
    generator_commit: str
    generator_settings: str
    intermediate_pdb_sha256: str
    strand_mode: str = "duplex"


@dataclass(frozen=True)
class _ElasticNetworkStats:
    total: int
    cross_chain: int


@dataclass(frozen=True)
class _CGAtom:
    index: int
    atom_name: str
    residue_name: str
    chain_id: str
    residue_number: int
    x_angstrom: float
    y_angstrom: float
    z_angstrom: float


@dataclass(frozen=True)
class ITPAtom:
    index: int
    atom_type: str
    residue_number: int
    residue_name: str
    atom_name: str
    charge_group: int
    charge: float
    mass: float


@dataclass(frozen=True)
class Interaction:
    section: str
    indices: tuple[int, ...]
    function: int | None
    line_number: int


@dataclass(frozen=True)
class ITPData:
    molecule_name: str
    atoms: tuple[ITPAtom, ...]
    interactions: tuple[Interaction, ...]

    @property
    def total_charge(self) -> float:
        return sum(atom.charge for atom in self.atoms)


@dataclass(frozen=True)
class GROAtom:
    residue_number: int
    residue_name: str
    atom_name: str
    index: int
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class GROData:
    title: str
    atoms: tuple[GROAtom, ...]
    box: tuple[float, ...]


@dataclass(frozen=True)
class ConsistencyReport:
    bead_count: int
    residue_count: int
    total_charge: float
    interaction_count: int


def normalize_martini_version(version: int | str) -> int:
    """Return a validated numeric Martini major version."""

    try:
        normalized = int(version)
    except (TypeError, ValueError) as exc:
        raise ConversionError(
            f"Unsupported Martini version {version!r}; choose 2 or 3"
        ) from exc
    if normalized not in MARTINI_VERSIONS or str(version).strip() not in {"2", "3"}:
        raise ConversionError(f"Unsupported Martini version {version!r}; choose 2 or 3")
    return normalized


_INDEX_COUNTS = {
    "bonds": 2,
    "angles": 3,
    "dihedrals": 4,
    "constraints": 2,
    "pairs": 2,
    "virtual_sites3": 4,
    "position_restraints": 1,
}


def parse_itp(path: Path) -> ITPData:
    """Parse the portions of a molecule ITP required for validation."""

    current_section: str | None = None
    molecule_name: str | None = None
    atoms: list[ITPAtom] = []
    interactions: list[Interaction] = []

    try:
        lines = Path(path).read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise OutputValidationError(f"Could not read generated ITP {path}: {exc}") from exc

    for line_number, raw_line in enumerate(lines, start=1):
        stripped = raw_line.strip()
        if not stripped or stripped.startswith(";") or stripped.startswith("#"):
            continue
        if stripped.startswith("[") and stripped.endswith("]"):
            current_section = stripped[1:-1].strip().lower()
            continue

        data = raw_line.split(";", 1)[0].split()
        if not data:
            continue
        try:
            if current_section == "moleculetype" and molecule_name is None:
                molecule_name = data[0]
            elif current_section == "atoms":
                if len(data) < 8:
                    raise ValueError("expected eight atom columns")
                atoms.append(
                    ITPAtom(
                        index=int(data[0]),
                        atom_type=data[1],
                        residue_number=int(data[2]),
                        residue_name=data[3],
                        atom_name=data[4],
                        charge_group=int(data[5]),
                        charge=float(data[6]),
                        mass=float(data[7]),
                    )
                )
            elif current_section == "exclusions":
                indices = tuple(int(token) for token in data)
                interactions.append(
                    Interaction(current_section, indices, None, line_number)
                )
            elif current_section in _INDEX_COUNTS:
                count = _INDEX_COUNTS[current_section]
                if len(data) < count:
                    raise ValueError(f"expected at least {count} index columns")
                indices = tuple(int(token) for token in data[:count])
                function = int(data[count]) if len(data) > count else None
                interactions.append(
                    Interaction(current_section, indices, function, line_number)
                )
        except ValueError as exc:
            raise OutputValidationError(
                f"Malformed [{current_section}] record in {path} at line {line_number}: {exc}"
            ) from exc

    if molecule_name is None:
        raise OutputValidationError(f"Generated ITP {path} has no [ moleculetype ] entry")
    if not atoms:
        raise OutputValidationError(f"Generated ITP {path} has no [ atoms ] records")
    return ITPData(molecule_name, tuple(atoms), tuple(interactions))


def parse_gro(path: Path) -> GROData:
    """Parse a coordinate-only GRO file using the standard fixed-width layout."""

    try:
        lines = Path(path).read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise OutputValidationError(f"Could not read generated GRO {path}: {exc}") from exc
    if len(lines) < 3:
        raise OutputValidationError(f"Generated GRO {path} is truncated")
    try:
        declared_count = int(lines[1].strip())
    except ValueError as exc:
        raise OutputValidationError(f"Invalid bead count in generated GRO {path}") from exc
    if len(lines) != declared_count + 3:
        raise OutputValidationError(
            f"GRO declares {declared_count} beads but has {len(lines) - 3} coordinate records"
        )

    atoms: list[GROAtom] = []
    for offset, line in enumerate(lines[2 : 2 + declared_count], start=3):
        if len(line) < 44:
            raise OutputValidationError(f"Malformed GRO coordinate record at line {offset}")
        try:
            atom = GROAtom(
                residue_number=int(line[0:5]),
                residue_name=line[5:10].strip(),
                atom_name=line[10:15].strip(),
                index=int(line[15:20]),
                x=float(line[20:28]),
                y=float(line[28:36]),
                z=float(line[36:44]),
            )
        except ValueError as exc:
            raise OutputValidationError(f"Malformed GRO field at line {offset}") from exc
        atoms.append(atom)
    try:
        box = tuple(float(value) for value in lines[-1].split())
    except ValueError as exc:
        raise OutputValidationError(f"Malformed box line in generated GRO {path}") from exc
    if len(box) not in {3, 9}:
        raise OutputValidationError("GRO box must contain three or nine values")
    return GROData(lines[0], tuple(atoms), box)


def validate_outputs(
    itp_path: Path,
    gro_path: Path,
    *,
    expected_name: str,
    expected_charge: float | None = None,
) -> ConsistencyReport:
    """Validate molecule identity, atom ordering, charge, and bonded indices."""

    itp = parse_itp(itp_path)
    gro = parse_gro(gro_path)
    if itp.molecule_name != expected_name:
        raise OutputValidationError(
            f"ITP molecule name {itp.molecule_name!r} does not match {expected_name!r}"
        )
    if len(itp.atoms) != len(gro.atoms):
        raise OutputValidationError(
            f"Bead-count mismatch: ITP has {len(itp.atoms)}, GRO has {len(gro.atoms)}"
        )
    expected_indices = list(range(1, len(itp.atoms) + 1))
    if [atom.index for atom in itp.atoms] != expected_indices:
        raise OutputValidationError("ITP atom indices are not contiguous from 1")
    if [atom.index for atom in gro.atoms] != expected_indices:
        raise OutputValidationError("GRO atom indices are not contiguous from 1")

    for itp_atom, gro_atom in zip(itp.atoms, gro.atoms):
        itp_identity = (
            itp_atom.residue_number,
            itp_atom.residue_name,
            itp_atom.atom_name,
        )
        gro_identity = (
            gro_atom.residue_number,
            gro_atom.residue_name,
            gro_atom.atom_name,
        )
        if itp_identity != gro_identity:
            raise OutputValidationError(
                f"ITP/GRO identity mismatch at bead {itp_atom.index}: "
                f"{itp_identity!r} != {gro_identity!r}"
            )
        if not all(isfinite(value) for value in (itp_atom.charge, itp_atom.mass)):
            raise OutputValidationError(f"Non-finite ITP value at bead {itp_atom.index}")
        if itp_atom.mass < 0:
            raise OutputValidationError(f"Negative mass at ITP bead {itp_atom.index}")
        expected_bead_charge = -1.0 if itp_atom.atom_name == "BB1" else 0.0
        if abs(itp_atom.charge - expected_bead_charge) > 1e-6:
            raise OutputValidationError(
                f"Unexpected charge {itp_atom.charge:g} e on {itp_atom.atom_name} "
                f"at ITP bead {itp_atom.index}; expected {expected_bead_charge:g} e"
            )
        if not all(isfinite(value) for value in (gro_atom.x, gro_atom.y, gro_atom.z)):
            raise OutputValidationError(f"Non-finite GRO coordinate at bead {gro_atom.index}")
        if itp_atom.charge_group < 1 or itp_atom.charge_group > len(itp.atoms):
            raise OutputValidationError(
                f"Charge-group index out of range at ITP bead {itp_atom.index}"
            )

    residue_numbers = [atom.residue_number for atom in itp.atoms]
    unique_residue_numbers = list(dict.fromkeys(residue_numbers))
    if unique_residue_numbers != list(range(1, max(unique_residue_numbers) + 1)):
        raise OutputValidationError("ITP residue ordering is not contiguous from 1")

    for interaction in itp.interactions:
        for index in interaction.indices:
            if index < 1 or index > len(itp.atoms):
                raise OutputValidationError(
                    f"[{interaction.section}] index {index} is outside 1..{len(itp.atoms)} "
                    f"at ITP line {interaction.line_number}"
                )

    if not all(isfinite(value) and value > 0 for value in gro.box[:3]):
        raise OutputValidationError("GRO box has non-positive or non-finite principal lengths")
    total_charge = itp.total_charge
    if expected_charge is not None and abs(total_charge - expected_charge) > 1e-6:
        raise OutputValidationError(
            f"Unexpected total charge: generated {total_charge:g} e, "
            f"expected {expected_charge:g} e"
        )

    return ConsistencyReport(
        bead_count=len(itp.atoms),
        residue_count=len(unique_residue_numbers),
        total_charge=total_charge,
        interaction_count=len(itp.interactions),
    )


def _git_commit(path: Path) -> str:
    try:
        completed = subprocess.run(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise ConversionError(
            f"Could not verify the AmberClassic source pin at {path}: {exc}"
        ) from exc
    if completed.returncode != 0:
        detail = completed.stderr.strip() or "git metadata is unavailable"
        raise ConversionError(
            f"Could not verify the AmberClassic source pin at {path}: {detail}"
        )
    return completed.stdout.strip().lower()


def resolve_amberclassic_home(home: str | Path | None = None) -> Path:
    """Locate and strictly verify the supported local AmberClassic checkout."""

    candidate: Path | None = Path(home) if home is not None else None
    if candidate is None and os.environ.get(AMBERCLASSIC_HOME_ENV):
        candidate = Path(os.environ[AMBERCLASSIC_HOME_ENV])
    if candidate is None:
        nab = shutil.which("nab")
        if nab:
            candidate = Path(nab).resolve().parent.parent
    if candidate is None:
        raise ConversionError(
            "Sequence input requires local AmberClassic NAB. Clone/build the pinned "
            f"{AMBERCLASSIC_TAG} source and set {AMBERCLASSIC_HOME_ENV}; "
            "the registration-gated 3DNA download and remote 3DNA API are not used"
        )

    candidate = candidate.expanduser().resolve()
    required = (
        candidate / "bin" / "nab",
        candidate / "bin" / "teLeap",
        candidate / "dat" / "leap",
    )
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise ConversionError(
            "AmberClassic NAB installation is incomplete; missing: " + ", ".join(missing)
        )
    actual_commit = _git_commit(candidate)
    if actual_commit != AMBERCLASSIC_COMMIT:
        raise ConversionError(
            "Unsupported AmberClassic source revision: "
            f"found {actual_commit or '<unknown>'}, require {AMBERCLASSIC_COMMIT} "
            f"({AMBERCLASSIC_TAG})"
        )
    return candidate


def _run_backend(
    command: list[str],
    *,
    work: Path,
    environment: dict[str, str],
    stage: str,
    verbose: bool,
) -> None:
    try:
        completed = subprocess.run(
            command,
            cwd=work,
            env=environment,
            text=True,
            capture_output=True,
            check=False,
            timeout=240,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise ConversionError(f"AmberClassic NAB {stage} failed: {exc}") from exc
    if verbose:
        if completed.stdout:
            print(completed.stdout, end="")
        if completed.stderr:
            print(completed.stderr, end="", file=sys.stderr)
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        if len(detail) > 4000:
            detail = detail[-4000:]
        raise ConversionError(
            f"AmberClassic NAB {stage} failed with exit code "
            f"{completed.returncode}: {detail or 'no diagnostic output'}"
        )


def _configure_upstream_logging(upstream: ModuleType, *, verbose: bool) -> None:
    upstream.logger.setLevel(logging.INFO if verbose else logging.WARNING)


def _parse_cg_pdb(path: Path) -> tuple[_CGAtom, ...]:
    atoms: list[_CGAtom] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ConversionError(f"Could not read temporary upstream CG PDB: {exc}") from exc
    for line_number, line in enumerate(lines, start=1):
        if line[0:6].strip() not in {"ATOM", "HETATM"}:
            continue
        try:
            atoms.append(
                _CGAtom(
                    index=int(line[6:11]),
                    atom_name=line[12:16].strip(),
                    residue_name=line[17:20].strip(),
                    chain_id=line[21:22].strip(),
                    residue_number=int(line[22:26]),
                    x_angstrom=float(line[30:38]),
                    y_angstrom=float(line[38:46]),
                    z_angstrom=float(line[46:54]),
                )
            )
        except (ValueError, IndexError) as exc:
            raise ConversionError(
                f"Malformed temporary CG PDB record at line {line_number}"
            ) from exc
    if not atoms:
        raise ConversionError("Pinned upstream converter generated no CG coordinates")
    return tuple(atoms)


def _elastic_network_stats(
    itp: ITPData, cg_atoms: tuple[_CGAtom, ...]
) -> _ElasticNetworkStats:
    chain_by_atom = {atom.index: atom.chain_id for atom in cg_atoms}
    elastic = [
        interaction
        for interaction in itp.interactions
        if interaction.section == "bonds" and interaction.function == 6
    ]
    cross_chain = sum(
        len({chain_by_atom[index] for index in interaction.indices}) > 1
        for interaction in elastic
    )
    return _ElasticNetworkStats(total=len(elastic), cross_chain=cross_chain)


def _remove_cross_chain_elastic_bonds(
    path: Path, cg_atoms: tuple[_CGAtom, ...]
) -> int:
    """Remove only function-6 elastic bonds that join different chains."""

    chain_by_atom = {atom.index: atom.chain_id for atom in cg_atoms}
    current_section: str | None = None
    output: list[str] = []
    removed = 0
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            current_section = stripped[1:-1].strip().lower()
        if current_section == "bonds" and stripped and not stripped.startswith((";", "#")):
            fields = raw_line.split(";", 1)[0].split()
            try:
                atom_a, atom_b, function = map(int, fields[:3])
            except (ValueError, IndexError):
                pass
            else:
                if function == 6 and chain_by_atom[atom_a] != chain_by_atom[atom_b]:
                    removed += 1
                    continue
        output.append(raw_line)
    path.write_text("\n".join(output) + "\n", encoding="utf-8")
    return removed


def _write_gro(
    path: Path,
    *,
    cg_atoms: tuple[_CGAtom, ...],
    itp: ITPData,
    title: str,
    padding_nm: float = 1.0,
) -> tuple[float, float, float]:
    if len(cg_atoms) > 99999:
        raise OutputValidationError("GRO output supports at most 99,999 beads")
    raw_nm = [
        (atom.x_angstrom / 10.0, atom.y_angstrom / 10.0, atom.z_angstrom / 10.0)
        for atom in cg_atoms
    ]
    minima = tuple(min(point[axis] for point in raw_nm) for axis in range(3))
    maxima = tuple(max(point[axis] for point in raw_nm) for axis in range(3))
    translation = tuple(padding_nm - minimum for minimum in minima)
    box = tuple(
        maxima[axis] - minima[axis] + 2.0 * padding_nm for axis in range(3)
    )

    lines = [title, f"{len(cg_atoms):5d}"]
    for cg_atom, itp_atom, point in zip(cg_atoms, itp.atoms, raw_nm):
        x, y, z = tuple(point[axis] + translation[axis] for axis in range(3))
        if max(abs(x), abs(y), abs(z)) >= 1000:
            raise OutputValidationError(
                "Coordinate magnitude is too large for standard GRO format"
            )
        lines.append(
            f"{itp_atom.residue_number:5d}{itp_atom.residue_name:<5.5s}"
            f"{itp_atom.atom_name:>5.5s}{itp_atom.index:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    lines.append("".join(f"{length:10.5f}" for length in box))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return translation


def _validate_coordinate_units(
    cg_atoms: tuple[_CGAtom, ...], gro_path: Path, translation: tuple[float, float, float]
) -> None:
    gro = parse_gro(gro_path)
    for cg_atom, gro_atom in zip(cg_atoms, gro.atoms):
        expected = (
            cg_atom.x_angstrom / 10.0 + translation[0],
            cg_atom.y_angstrom / 10.0 + translation[1],
            cg_atom.z_angstrom / 10.0 + translation[2],
        )
        actual = (gro_atom.x, gro_atom.y, gro_atom.z)
        if any(abs(left - right) > 0.00051 for left, right in zip(expected, actual)):
            raise OutputValidationError(
                f"PDB Å-to-GRO nm conversion mismatch at bead {cg_atom.index}"
            )


def _chain_provenance(summary: Any) -> str:
    descriptions = []
    chains = summary.residues_by_chain()
    source_id_counts = Counter(chain[0].chain_id for chain in chains)
    for chain in chains:
        chain_id = chain[0].chain_id or "<blank>"
        if not chain[0].chain_id or source_id_counts[chain[0].chain_id] > 1:
            chain_id += f" (segment {chain[0].segment_index + 1})"
        descriptions.append(
            f"{chain_id}: input {chain[0].source_label}-{chain[-1].source_label} "
            f"=> residues {chain[0].output_number}-{chain[-1].output_number}"
        )
    return "; ".join(descriptions)


def _publish_output_pair(
    stage_itp: Path,
    stage_gro: Path,
    final_itp: Path,
    final_gro: Path,
    work: Path,
) -> None:
    """Publish both files with rollback if either atomic rename fails."""

    backups: list[tuple[Path, Path]] = []
    published: list[tuple[Path, Path]] = []
    try:
        for final in (final_itp, final_gro):
            if final.exists():
                backup = work / f"previous-{final.name}"
                os.replace(final, backup)
                backups.append((final, backup))
        for staged, final in ((stage_itp, final_itp), (stage_gro, final_gro)):
            os.replace(staged, final)
            published.append((staged, final))
    except OSError as exc:
        rollback_errors: list[str] = []
        for staged, final in reversed(published):
            try:
                if final.exists():
                    os.replace(final, staged)
            except OSError as rollback_exc:  # pragma: no cover - exceptional filesystem failure
                rollback_errors.append(f"{final}: {rollback_exc}")
        for final, backup in reversed(backups):
            try:
                if backup.exists():
                    os.replace(backup, final)
            except OSError as rollback_exc:  # pragma: no cover - exceptional filesystem failure
                rollback_errors.append(f"{final}: {rollback_exc}")
        detail = ""
        if rollback_errors:
            detail = "; rollback also failed for " + "; ".join(rollback_errors)
        raise ConversionError(f"Could not publish the ITP/GRO output pair: {exc}{detail}") from exc
