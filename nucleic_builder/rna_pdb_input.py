"""Strict Phase 1 PDB validation and normalization.

The vendored converter intentionally has permissive residue aliases and a few
PDB parsing behaviours that are unsuitable for a validation boundary.  This
module establishes that boundary before any upstream code is called.
"""

from __future__ import annotations

import hashlib
from collections import OrderedDict
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
from typing import Iterable

from .errors import InputValidationError


CANONICAL_RESIDUES = frozenset({"A", "C", "G", "U"})

# These are the terminal aliases verified end-to-end against the pinned
# upstream implementation.  They are normalized in a temporary copy so the
# user input itself is never changed.
TERMINAL_ALIASES = {
    **{f"R{base}5": (base, "5") for base in CANONICAL_RESIDUES},
    **{f"R{base}3": (base, "3") for base in CANONICAL_RESIDUES},
    **{f"{base}5": (base, "5") for base in CANONICAL_RESIDUES},
    **{f"{base}3": (base, "3") for base in CANONICAL_RESIDUES},
}

ACCEPTED_RESIDUES = {
    **{base: (base, None) for base in CANONICAL_RESIDUES},
    **TERMINAL_ALIASES,
}

# The upstream mapping aliases these modified nucleotides to canonical bases.
# Phase 1 rejects them explicitly: aliasing is not parameterization.
KNOWN_MODIFIED_RESIDUES = frozenset(
    {
        "2MA",
        "6MA",
        "RAP",
        "DMA",
        "DHA",
        "SPA",
        "5MC",
        "3MP",
        "MRC",
        "NMC",
        "1MG",
        "2MG",
        "7MG",
        "MRG",
        "4SU",
        "DHU",
        "PSU",
        "5MU",
        "3MU",
        "MRU",
    }
)

# Ordinary crystallographic waters are not part of the requested RNA molecule.
# They are removed from the private normalized copy and reported to the caller.
WATER_RESIDUES = frozenset({"HOH", "WAT"})

# The upstream mapper creates a bead when *any* listed atom is present.  That
# permissive behaviour is useful during model development but unsafe at this
# wrapper boundary: a missing heavy atom would silently move the bead centre.
# Hydrogens remain optional because crystallographic PDBs commonly omit them.
_SUGAR_HEAVY_ATOMS = frozenset({"C5'", "C4'", "O4'", "C3'", "C1'", "C2'", "O2'"})
_BASE_HEAVY_ATOMS = {
    "A": frozenset({"N9", "C8", "N3", "C4", "N1", "C2", "N6", "C6", "N7", "C5"}),
    "C": frozenset({"N1", "C5", "C6", "C2", "O2", "N3", "N4", "C4"}),
    "G": frozenset(
        {"C8", "N9", "C4", "N3", "C2", "N2", "N1", "C6", "O6", "C5", "N7"}
    ),
    "U": frozenset({"N1", "C5", "C6", "C2", "O2", "N3", "C4", "O4"}),
}

# PDB chain identifiers have one column.  The normalized copy assigns a unique
# identifier to every chain segment, including TER-separated segments whose
# source chain identifiers are blank or reused.
_NORMALIZED_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


@dataclass(frozen=True)
class ResidueInfo:
    """One input residue and its globally renumbered topology position."""

    chain_id: str
    segment_index: int
    normalized_chain_id: str
    source_number: int
    insertion_code: str
    source_name: str
    canonical_name: str
    output_number: int

    @property
    def source_label(self) -> str:
        suffix = self.insertion_code or ""
        return f"{self.source_number}{suffix}"


@dataclass(frozen=True)
class InputSummary:
    """Validated structural metadata retained across conversion."""

    source: Path
    input_sha256: str
    residues: tuple[ResidueInfo, ...]
    ignored_water_atoms: int

    @property
    def chain_ids(self) -> tuple[str, ...]:
        """Unique internal chain IDs used by the normalized PDB."""

        return tuple(dict.fromkeys(residue.normalized_chain_id for residue in self.residues))

    @property
    def chain_count(self) -> int:
        return len(self.chain_ids)

    @property
    def residue_count(self) -> int:
        return len(self.residues)

    @property
    def sequence(self) -> str:
        return "".join(residue.canonical_name for residue in self.residues)

    @property
    def expected_charge(self) -> float:
        # The Martini RNA mapping omits the 5' phosphate bead on each chain.
        return float(-(self.residue_count - self.chain_count))

    def residues_by_chain(self) -> tuple[tuple[ResidueInfo, ...], ...]:
        groups: list[list[ResidueInfo]] = []
        segment_indices = tuple(dict.fromkeys(res.segment_index for res in self.residues))
        for segment_index in segment_indices:
            groups.append([res for res in self.residues if res.segment_index == segment_index])
        return tuple(tuple(group) for group in groups)


@dataclass(frozen=True)
class _CoordinateRecord:
    line_number: int
    line: str
    segment_index: int
    chain_id: str
    residue_number: int
    insertion_code: str
    residue_name: str
    atom_name: str

    @property
    def residue_key(self) -> tuple[int, str, int, str]:
        return (
            self.segment_index,
            self.chain_id,
            self.residue_number,
            self.insertion_code,
        )


def _chain_label(chain_id: str) -> str:
    return chain_id if chain_id else "<blank>"


def _format_locations(records: Iterable[_CoordinateRecord]) -> str:
    locations = []
    seen = set()
    for record in records:
        key = (record.residue_name, record.residue_key)
        if key in seen:
            continue
        seen.add(key)
        locations.append(
            f"{record.residue_name} at chain {_chain_label(record.chain_id)} "
            f"residue {record.residue_number}{record.insertion_code}"
        )
    return ", ".join(locations)


def _validate_required_heavy_atoms(
    summary: InputSummary,
    residue_atom_names: dict[tuple[int, str, int, str], set[str]],
) -> None:
    """Reject partial Martini mapping groups before invoking upstream code."""

    for chain_residues in summary.residues_by_chain():
        for position, residue in enumerate(chain_residues):
            key = (
                residue.segment_index,
                residue.chain_id,
                residue.source_number,
                residue.insertion_code,
            )
            present = residue_atom_names[key]
            required = set(_SUGAR_HEAVY_ATOMS | _BASE_HEAVY_ATOMS[residue.canonical_name])

            # move_o3() in the authors' converter transfers each residue's O3'
            # into the following residue's BB1.  On a multi-residue chain it
            # also needs the terminal O3' as the trigger that transfers the
            # penultimate O3'; the terminal coordinate itself is discarded.
            if len(chain_residues) > 1:
                required.add("O3'")
            if position > 0:
                required.update({"P", "O5'"})

            missing = sorted(required - present)
            if position > 0 and not ({"OP1", "O1P"} & present):
                missing.append("OP1/O1P")
            if position > 0 and not ({"OP2", "O2P"} & present):
                missing.append("OP2/O2P")
            if missing:
                raise InputValidationError(
                    "Incomplete all-atom mapping for "
                    f"{residue.source_name} at chain {_chain_label(residue.chain_id)} "
                    f"residue {residue.source_label}: missing required mapped heavy atom(s): "
                    + ", ".join(missing)
                )


def prepare_input_pdb(source: Path, destination: Path) -> InputSummary:
    """Validate *source* and write the normalized RNA-only PDB to *destination*.

    Residue and chain record order are retained.  Residues are numbered globally
    in a private copy because the upstream topology always uses global 1..N
    residue numbers, while many duplex PDBs restart numbering in each chain.
    """

    source = Path(source)
    if not source.is_file():
        raise InputValidationError(f"Input PDB does not exist or is not a file: {source}")

    try:
        source_bytes = source.read_bytes()
    except OSError as exc:
        raise InputValidationError(f"Could not read input PDB {source}: {exc}") from exc
    try:
        lines = source_bytes.decode("utf-8").splitlines()
    except UnicodeError as exc:
        raise InputValidationError(f"Input PDB is not valid UTF-8 text: {source}") from exc
    input_sha256 = hashlib.sha256(source_bytes).hexdigest()

    explicit_model_count = 0
    in_explicit_model = False
    saw_explicit_model = False
    coordinate_records: list[_CoordinateRecord] = []
    unsupported_records: list[_CoordinateRecord] = []
    modified_records: list[_CoordinateRecord] = []
    ignored_water_atoms = 0

    residues: OrderedDict[tuple[int, str, int, str], str] = OrderedDict()
    residue_atom_names: dict[tuple[int, str, int, str], set[str]] = {}
    finished_residue_keys: set[tuple[int, str, int, str]] = set()
    previous_residue_key: tuple[int, str, int, str] | None = None
    current_segment_index: int | None = None
    current_source_chain: str | None = None
    seen_source_chain_blocks: set[str] = set()
    ter_break_pending = False

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.rstrip("\r\n")
        record_name = line[0:6].strip().upper() if line else ""

        if record_name == "MODEL":
            explicit_model_count += 1
            saw_explicit_model = True
            in_explicit_model = True
            if explicit_model_count > 1:
                raise InputValidationError(
                    "Phase 1 accepts exactly one PDB model; multiple MODEL records were found"
                )
            continue
        if record_name == "ENDMDL":
            in_explicit_model = False
            continue
        if record_name == "TER":
            ter_break_pending = True
            continue
        if record_name not in {"ATOM", "HETATM"}:
            continue
        if saw_explicit_model and not in_explicit_model:
            raise InputValidationError(
                f"Coordinate record outside the single MODEL block at line {line_number}"
            )
        if len(line) < 54:
            raise InputValidationError(
                f"Malformed PDB coordinate record at line {line_number}: expected at least 54 columns"
            )

        try:
            int(line[6:11])
            residue_number = int(line[22:26])
            xyz = tuple(float(line[start:end]) for start, end in ((30, 38), (38, 46), (46, 54)))
        except ValueError as exc:
            raise InputValidationError(
                f"Malformed numeric field in PDB coordinate record at line {line_number}"
            ) from exc
        if not all(isfinite(value) for value in xyz):
            raise InputValidationError(f"Non-finite coordinate at PDB line {line_number}")

        alt_loc = line[16:17].strip()
        if alt_loc:
            raise InputValidationError(
                f"Alternate location {alt_loc!r} at PDB line {line_number} is unsupported; "
                "provide a PDB with one selected conformer"
            )

        residue_name = line[17:20].strip().upper()
        atom_name = line[12:16].strip()
        chain_id = line[21:22].strip()
        insertion_code = line[26:27].strip()
        provisional_segment = current_segment_index if current_segment_index is not None else 0
        record = _CoordinateRecord(
            line_number=line_number,
            line=line,
            segment_index=provisional_segment,
            chain_id=chain_id,
            residue_number=residue_number,
            insertion_code=insertion_code,
            residue_name=residue_name,
            atom_name=atom_name,
        )

        if record_name == "HETATM" and residue_name in WATER_RESIDUES:
            ignored_water_atoms += 1
            continue
        if residue_name in KNOWN_MODIFIED_RESIDUES:
            modified_records.append(record)
            continue
        if residue_name not in ACCEPTED_RESIDUES:
            unsupported_records.append(record)
            continue

        starts_segment = (
            current_segment_index is None
            or ter_break_pending
            or chain_id != current_source_chain
        )
        if starts_segment:
            if (
                current_segment_index is not None
                and not ter_break_pending
                and chain_id in seen_source_chain_blocks
            ):
                raise InputValidationError(
                    f"Chain {_chain_label(chain_id)} appears in multiple non-contiguous blocks"
                )
            if previous_residue_key is not None:
                finished_residue_keys.add(previous_residue_key)
            if current_segment_index is None:
                current_segment_index = 0
            else:
                current_segment_index += 1
            if current_segment_index >= len(_NORMALIZED_CHAIN_IDS):
                raise InputValidationError(
                    f"Phase 1 supports at most {len(_NORMALIZED_CHAIN_IDS)} PDB chain segments"
                )
            current_source_chain = chain_id
            seen_source_chain_blocks.add(chain_id)
            previous_residue_key = None
            ter_break_pending = False

        assert current_segment_index is not None
        residue_key = (current_segment_index, chain_id, residue_number, insertion_code)
        record = _CoordinateRecord(
            line_number=line_number,
            line=line,
            segment_index=current_segment_index,
            chain_id=chain_id,
            residue_number=residue_number,
            insertion_code=insertion_code,
            residue_name=residue_name,
            atom_name=atom_name,
        )

        if residue_key != previous_residue_key:
            if residue_key in finished_residue_keys:
                raise InputValidationError(
                    f"Residue records are not contiguous for chain {_chain_label(chain_id)} "
                    f"residue {residue_number}{insertion_code}"
                )
            if previous_residue_key is not None:
                finished_residue_keys.add(previous_residue_key)
            previous_residue_key = residue_key

        old_name = residues.get(residue_key)
        if old_name is not None and old_name != residue_name:
            raise InputValidationError(
                f"Conflicting names {old_name!r} and {residue_name!r} for chain "
                f"{_chain_label(chain_id)} residue {residue_number}{insertion_code}"
            )
        residues.setdefault(residue_key, residue_name)

        seen_names = residue_atom_names.setdefault(residue_key, set())
        if atom_name in seen_names:
            raise InputValidationError(
                f"Duplicate atom name {atom_name!r} in chain {_chain_label(chain_id)} "
                f"residue {residue_number}{insertion_code}; alternate conformers must be resolved"
            )
        seen_names.add(atom_name)
        coordinate_records.append(record)

    if modified_records:
        raise InputValidationError(
            "Modified nucleotides are not supported in Phase 1 and will not be mapped to "
            f"canonical bases: {_format_locations(modified_records)}"
        )
    if unsupported_records:
        raise InputValidationError(
            "Unsupported non-canonical residue(s) in input PDB: "
            f"{_format_locations(unsupported_records)}. Accepted RNA residue names are "
            "A, C, G, U and verified 3'/5' terminal variants."
        )
    if not coordinate_records:
        raise InputValidationError("Input PDB contains no canonical RNA atom records")

    residue_infos: list[ResidueInfo] = []
    for output_number, (
        (segment_index, chain_id, source_number, icode),
        source_name,
    ) in enumerate(
        residues.items(), start=1
    ):
        canonical_name, terminus = ACCEPTED_RESIDUES[source_name]
        residue_infos.append(
            ResidueInfo(
                chain_id=chain_id,
                segment_index=segment_index,
                normalized_chain_id=_NORMALIZED_CHAIN_IDS[segment_index],
                source_number=source_number,
                insertion_code=icode,
                source_name=source_name,
                canonical_name=canonical_name,
                output_number=output_number,
            )
        )

    summary = InputSummary(
        source=source,
        input_sha256=input_sha256,
        residues=tuple(residue_infos),
        ignored_water_atoms=ignored_water_atoms,
    )

    for chain_residues in summary.residues_by_chain():
        for position, residue in enumerate(chain_residues):
            _, terminus = ACCEPTED_RESIDUES[residue.source_name]
            if terminus == "5" and position != 0:
                raise InputValidationError(
                    f"5' terminal name {residue.source_name} is not first in chain "
                    f"{_chain_label(residue.chain_id)}"
                )
            if terminus == "3" and position != len(chain_residues) - 1:
                raise InputValidationError(
                    f"3' terminal name {residue.source_name} is not last in chain "
                    f"{_chain_label(residue.chain_id)}"
                )

    _validate_required_heavy_atoms(summary, residue_atom_names)

    if summary.residue_count > 9999:
        raise InputValidationError("Phase 1 PDB/GRO output supports at most 9,999 residues")

    residue_by_key = {
        (res.segment_index, res.chain_id, res.source_number, res.insertion_code): res
        for res in residue_infos
    }
    normalized_lines: list[str] = []
    previous_segment: int | None = None
    for record in coordinate_records:
        residue = residue_by_key[record.residue_key]
        if previous_segment is not None and residue.segment_index != previous_segment:
            normalized_lines.append("TER")
        previous_segment = residue.segment_index

        padded = record.line.ljust(80)
        normalized = (
            padded[:16]
            + " "
            + f"{residue.canonical_name:>3s}"
            + padded[20:21]
            + residue.normalized_chain_id
            + f"{residue.output_number:4d}"
            + padded[26:]
        )
        normalized_lines.append(normalized.rstrip())
    normalized_lines.extend(("TER", "END"))

    try:
        destination.write_text("\n".join(normalized_lines) + "\n", encoding="utf-8")
    except OSError as exc:
        raise InputValidationError(f"Could not write temporary normalized PDB: {exc}") from exc

    return summary
