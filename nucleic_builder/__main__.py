"""Command-line interface for ``python -m nucleic_builder``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .core import ELASTIC_NETWORK_POLICIES
from .dna_builder import DNA_MODEL_WARNING, build_dna
from .dna_sequence_builder import build_dna_from_sequence
from .errors import BuilderError
from .martini2_builder import MARTINI2_MODEL_WARNING
from .rna_builder import build_rna
from .rna_sequence_builder import build_rna_from_sequence


def _program_name() -> str:
    if Path(sys.argv[0]).stem == "nucleic-builder":
        return "nucleic-builder"
    return "python -m nucleic_builder"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=_program_name(),
        description=(
            "Convert canonical RNA or DNA into molecule-level ITP and GRO files. "
            "Martini 3 is the default; its DNA backend is experimental DNA-alpha."
        ),
    )
    parser.add_argument(
        "--martini-version",
        choices=(2, 3),
        default=3,
        type=int,
        help="Martini major version; defaults to 3",
    )
    parser.add_argument(
        "--polymer",
        choices=("rna", "dna"),
        default="rna",
        help=(
            "polymer selection; defaults to rna. DNA must be selected explicitly; "
            "Martini 3 DNA is experimental/unpublished, while Martini 2 DNA is published"
        ),
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input", type=Path, help="existing canonical all-atom PDB")
    source.add_argument(
        "--sequence",
        help=(
            "one canonical strand, 5' to 3' (RNA: A/C/G/U, DNA: A/C/G/T); "
            "RNA strand count is selected with --strand-mode, while DNA always "
            "builds its antiparallel complement"
        ),
    )
    parser.add_argument(
        "--strand-mode",
        choices=("duplex", "single"),
        help=(
            "RNA sequence input only; 'single' retains exactly the entered strand "
            "in an idealized A-form-derived starting conformation, while 'duplex' "
            "builds the current antiparallel duplex (default: duplex)"
        ),
    )
    parser.add_argument("--name", required=True, help="GROMACS molecule and output basename")
    parser.add_argument("--output-dir", required=True, type=Path, help="output directory")
    parser.add_argument(
        "--force", action="store_true", help="replace existing NAME.itp and NAME.gro"
    )
    parser.add_argument(
        "--elastic-network",
        choices=sorted(ELASTIC_NETWORK_POLICIES),
        required=True,
        help=(
            "elastic-network policy: legacy retains the selected backend's network "
            "including cross-chain bonds; intrachain removes cross-chain elastic "
            "bonds; off disables the network"
        ),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="show detailed messages from the structure generator and converter",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    if args.strand_mode is not None and (
        args.polymer != "rna" or args.sequence is None
    ):
        parser.error("--strand-mode is valid only with RNA --sequence input")
    if args.martini_version == 2:
        print(f"WARNING: {MARTINI2_MODEL_WARNING}", file=sys.stderr)
    elif args.polymer == "dna":
        print(f"WARNING: {DNA_MODEL_WARNING}", file=sys.stderr)
    try:
        if args.polymer == "dna" and args.sequence is not None:
            result = build_dna_from_sequence(
                args.sequence,
                args.name,
                args.output_dir,
                force=args.force,
                elastic_network=args.elastic_network,
                verbose=args.verbose,
                martini_version=args.martini_version,
            )
        elif args.polymer == "dna":
            result = build_dna(
                args.input,
                args.name,
                args.output_dir,
                force=args.force,
                elastic_network=args.elastic_network,
                verbose=args.verbose,
                martini_version=args.martini_version,
            )
        elif args.sequence is not None:
            result = build_rna_from_sequence(
                args.sequence,
                args.name,
                args.output_dir,
                force=args.force,
                elastic_network=args.elastic_network,
                strand_mode=args.strand_mode or "duplex",
                verbose=args.verbose,
                martini_version=args.martini_version,
            )
        else:
            result = build_rna(
                args.input,
                args.name,
                args.output_dir,
                force=args.force,
                elastic_network=args.elastic_network,
                verbose=args.verbose,
                martini_version=args.martini_version,
            )
    except BuilderError as exc:
        print(f"nucleic-builder: error: {exc}", file=sys.stderr)
        return 2

    if result.ignored_water_atoms:
        if result.polymer_type == "rna":
            print(
                f"Ignored {result.ignored_water_atoms} crystallographic water atom(s); "
                "only the RNA molecule was converted.",
                file=sys.stderr,
            )
        else:
            print(
                f"Ignored {result.ignored_water_atoms} crystallographic water atom(s); "
                "only the DNA molecule was converted.",
                file=sys.stderr,
            )
    print(f"Wrote {result.itp_path}")
    print(f"Wrote {result.gro_path}")
    print(
        f"Validated {result.bead_count} beads, {result.residue_count} residues, "
        f"{result.chain_count} chain(s), charge {result.total_charge:g} e"
    )
    print(
        f"Elastic network: {result.elastic_network} "
        f"({result.elastic_bond_count} bonds; "
        f"{result.cross_chain_elastic_bond_count} cross-chain)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
