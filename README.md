# Nucleic Builder — Martini 3 RNA/DNA

**Validated Phase 2B deliverable.** `nucleic_builder` converts a canonical RNA or
DNA sequence—or an existing canonical all-atom PDB—into exactly two
molecule-level GROMACS files:

- `NAME.itp` — molecule topology
- `NAME.gro` — coarse-grained coordinates in nm

## Current deliverable

| Input | Structure used | Martini model | Output |
|---|---|---|---|
| RNA sequence, `A/C/G/U` | Ideal A-form duplex | Published Martini 3 RNA | `NAME.itp`, `NAME.gro` |
| DNA sequence, `A/C/G/T` | Ideal B-form duplex | **Experimental DNA-alpha** | `NAME.itp`, `NAME.gro` |
| Canonical RNA PDB | Supplied geometry | Published Martini 3 RNA | `NAME.itp`, `NAME.gro` |
| Canonical DNA PDB | Supplied geometry | **Experimental DNA-alpha** | `NAME.itp`, `NAME.gro` |

Sequence input currently constructs a duplex; it does not predict whether a
sequence folds, hybridizes, or remains stable in solution. Existing-PDB input
also supports conversion tests for a single strand or known hairpin, but no
single-strand or hairpin structure is predicted from sequence. Modified
nucleotides, RNA/DNA hybrids, user-facing solvation, and structure prediction
are deliberately out of scope.

The program does not recreate or fit force-field parameters. RNA remains a
strict wrapper around the implementation published by Danis Yangaliev and
S. Banu Ozkan. Its behavior and pin remain unchanged at
[`e761b7349fdf61dd485053c000dbb642f24ff9d8`](https://github.com/DanYev/Martini-3-DNA-RNA/commit/e761b7349fdf61dd485053c000dbb642f24ff9d8).
DNA is an isolated wrapper around that snapshot's separately named
`martinize_dna_alpha.py` and `dna_alpha_itps/`, independently recorded at the
same exact commit. DNA never enters the RNA converter or RNA parameter set.

The importable Python package is named `nucleic_builder`; its public API covers
both model paths.

## Experimental DNA-alpha warning

> **DNA output is EXPERIMENTAL / UNPUBLISHED DNA-ALPHA.** The authors' upstream
> README says the DNA force field is alpha and that a publication is being
> prepared. The published Martini 3 paper cited below is an RNA paper, not DNA
> validation. Successful conversion and a short GROMACS smoke test establish
> file compatibility and numerical execution only—not scientific accuracy,
> production readiness, duplex stability, melting behavior, or transferability.

DNA must be selected explicitly with `--polymer dna`. The warning is written
into the ITP provenance, the GRO title, and every DNA CLI run.

## Scientific scope warning

**Successful conversion and a short numerically stable trajectory do not show
that the model is suitable for every RNA question.** The model paper identifies
the lack of directional hydrogen bonding as a limitation for studying
**hybridization, melting, and intercalation**. It also says the model is
optimized for dsRNA and needs further refinement for **ssRNA and complex
secondary structures**, including tRNA-like folds. Do not use these smoke-test
results as validation of those phenomena. See the
[published paper](https://doi.org/10.1016/j.bpj.2025.07.034) and
[public full preprint](https://doi.org/10.1101/2025.04.13.648640).
An ideal A-form starting structure does not remove any of these force-field
limitations and is not evidence that a sequence will hybridize or remain
duplexed in solution.

## Install

Python 3.9 or newer and NumPy are required. Existing-PDB conversion has no
structure-generator dependency. For development and tests:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[test]'
```

The installed console command is `nucleic-builder`; `python -m nucleic_builder`
is the equivalent module form used in the examples below.

GROMACS is not required to generate molecule files. A working `gmx`
executable is required for the integration tests.

### Sequence backend

For RNA, 3DNA v2.3 `fiber` was preferred, but
[3DNA's own v2.3 release note](https://x3dna.org/highlights/3dna-c-source-code-is-available)
says users must register on its forum to download it. It therefore cannot be
installed non-interactively or redistributed with this project. Phase 2A uses
the requested open-source fallback:
AmberClassic NAB `fd_helix("arna", sequence, "rna")`. DNA uses the same pinned
NAB installation with its documented right-handed Arnott B-DNA setting,
`fd_helix("abdna", sequence, "dna")`. Neither path calls a remote 3DNA API.

The backend is pinned to AmberClassic tag `v2.0`, commit
`bdb3e0dee5b90f2be2950e26cfad1ae5a7440cae`. Install it from the official
source with GNU C/C++/Fortran compilers, Make, Bison, and Flex available:

```bash
scripts/install_amberclassic_nab.sh "$PWD/.tools/AmberClassic-v2.0"
export NUCLEIC_BUILDER_AMBERCLASSIC_HOME="$PWD/.tools/AmberClassic-v2.0"
```

That installer and the supported NAB build path are for Linux or WSL2. Native
Windows is not a supported AmberClassic/NAB backend for this project.

The installer refuses to overwrite a target, checks out the exact commit,
uses the upstream build, and verifies `nab`, `teLeap`, and the Git revision.
The Python runtime repeats those checks and rejects a different or incomplete
backend. The inspected AmberClassic binary package was not used because it
omitted the `nab` translator even though it contained NAB runtime data.

## Quick start

Every command requires an explicit elastic-network choice. For an ordinary
sequence-derived duplex, use `off` unless a specific validated protocol calls
for a different policy.

### RNA sequence to ideal A-form duplex

```bash
python -m nucleic_builder \
  --sequence GCAUCG \
  --name SEQRNA \
  --output-dir output \
  --elastic-network off
```

This means strand A is `5′-GCAUCG-3′`. The generated strand B is
`5′-CGAUGC-3′`, paired antiparallel as `3′-CGUAGC-5′`.

### DNA sequence to ideal B-form duplex

```bash
python -m nucleic_builder \
  --polymer dna \
  --sequence GCATCG \
  --name SEQDNA \
  --output-dir output \
  --elastic-network off
```

This means strand A is `5′-GCATCG-3′`; strand B is `5′-CGATGC-3′`,
paired as `3′-CGTAGC-5′`. DNA must be selected explicitly and every run
and generated file is marked experimental/unpublished.

### Existing canonical all-atom PDB

```bash
python -m nucleic_builder \
  --input structure.pdb \
  --name RNA \
  --output-dir output \
  --elastic-network off
```

RNA is the default polymer. Use `--polymer dna --input structure.pdb` for
canonical DNA. `--input` and
`--sequence` are mutually exclusive. Sequence matching is case-insensitive,
but whitespace, FASTA text, ambiguity codes, the other polymer's base spelling,
and modified-base notation are rejected.

Default logging is quiet: only the result summary and actionable warnings or
errors are shown. Add `--verbose` to show detailed messages from the pinned
structure generator and upstream converter.

The command refuses to replace an existing `NAME.itp` or `NAME.gro`; pass
`--force` when replacement is intentional. It creates only those two molecule
files. Temporary atomic/CG PDBs are used internally as needed and deleted
before exit.

The molecule name must start with a letter and may contain letters, digits,
`.`, `_`, `+`, and `-`. It is used exactly as the output basename and the ITP
`[ moleculetype ]` name.

### Elastic-network policy

The CLI requires an explicit choice. It does not infer a structural class or
silently select the upstream `legacy` behavior:

```bash
python -m nucleic_builder \
  --input structure.pdb \
  --name RNA \
  --output-dir output \
  --elastic-network intrachain
```

| Value | Behavior |
|---|---|
| `legacy` | Authors' 0.3–1.2 nm, 200 kJ mol⁻¹ nm⁻² network on BB1/BB3 beads, including possible cross-chain bonds |
| `intrachain` | Same network, with only cross-chain function-6 elastic bonds removed |
| `off` | No elastic network |

For the published RNA model, recommended starting policies are `off` for an ordinary duplex or flexible
single strand, and `intrachain` for a known hairpin whose experimental fold
must be retained. Use `off` instead for a hairpin folding/unfolding question.
Use `legacy` only when Phase 1 reproducibility or deliberate cross-strand
locking is required. These recommendations are explicit rather than
automatically inferred because a PDB does not reliably label its structural
class.

No scientifically validated default is claimed for the unpublished DNA-alpha
model. The DNA integration matrix exercises authors-compatible `legacy` and a
sequence-derived `off` case only to validate both code paths. The CLI still
requires the user to choose explicitly.

### Accepted input

The default RNA PDB route accepts a single model containing unmodified RNA residues:

```text
A C G U
RA5 RA3 A5 A3   RC5 RC3 C5 C3
RG5 RG3 G5 G3   RU5 RU3 U5 U3
```

Terminal aliases are checked for the correct chain end and normalized to
`A/C/G/U` in a private temporary copy. Residue and chain order are retained;
residues are numbered globally in the output so duplexes whose chains both
start at residue 1 remain consistent with GROMACS topology numbering.

The sequence route accepts one non-empty strand containing only `A/C/G/U`.
It is duplex construction, not a sequence-file parser: no FASTA header,
whitespace cleanup, DNA spelling, single-strand output, or folding is implied.

`HETATM` records named `HOH` or `WAT` are removed from the private RNA-only
copy and reported by the CLI. Other non-RNA residues are rejected. Known
modified bases—including `PSU`, `5MC`, `1MG`, and `2MA`—are rejected and are
never mapped to canonical bases. DNA residue names, sequence-only text,
multiple models, unresolved alternate locations, interleaved chain blocks,
and incomplete heavy-atom mapping groups fail with a clear diagnostic.
Hydrogens are optional.

With `--polymer dna`, the PDB route instead accepts canonical DNA residue names
`DA/DC/DG/DT`, plain `A/C/G/T`, and ordinary `DA5/DA3` or `A5/A3` terminal
variants (and corresponding C/G/T names). Terminal positions are checked and a
private copy is normalized to `A/C/G/T`. An RNA ribose `O2'` atom is rejected,
so an A/C/G-only RNA cannot silently pass as DNA. `U`, RNA terminal names, proteins,
ligands, ions, and modified DNA nucleotides are rejected. The DNA sequence
route accepts exactly one non-empty `A/C/G/T` strand and constructs only its
ideal antiparallel B-form duplex.

### Coordinates, chains, and provenance

The authors' atom-to-bead arithmetic averages are preserved. Å coordinates
from the temporary CG PDB are converted to nm, then the molecule is moved by a
single rigid translation into an orthorhombic box with 1 nm padding. A rigid
translation does not change molecular geometry.

GRO has no chain-ID field. Source chain changes and explicit PDB `TER` records
define chain segments, preserving separate strands even when source chain IDs
are blank or reused. Chain order, breaks, source IDs, selected elastic policy,
elastic and cross-chain bond counts, the exact upstream pin, and the SHA-256
of the original input bytes are recorded in the ITP provenance header. No
ordinary backbone or side-chain bonded term crosses a chain break. The ITP
also retains the authors' conditional `POSRES` block.

For sequence input, provenance instead records both strands in their 5′→3′
forms, the paired 3′→5′ orientation, AmberClassic repository/tag/commit and
NAB settings, the private all-atom PDB SHA-256, and the unchanged Martini
upstream provenance. The temporary NAB source, executable, atomic PDB, and CG
PDB are deleted; only `NAME.itp` and `NAME.gro` are published.

Every DNA ITP additionally records the separate DNA-alpha converter, parameter
directory, independent upstream pin, and experimental/unpublished status. The
GRO title carries the same warning. RNA and DNA provenance identify their
polymer and cannot silently share a converter.

## Use the molecule files in GROMACS

The generated ITP contains one molecule type, not a complete system topology.
An RNA simulation topology must include Martini 3 particle definitions, the
authors' RNA nonbonded parameters, ions/solvent as needed, and the generated
file in this order:

```gromacs
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_rna.itp"
#include "martini_v3.0.0_ions_v1.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "RNA.itp"

[ system ]
RNA test system

[ molecules ]
RNA  1
```

Experimental DNA uses `martini_alpha_dna.itp` in place of
`martini_v3.0.0_rna.itp` and includes the DNA molecule ITP afterward. Never
include the RNA nonbonded file as a substitute for the DNA-alpha file. The
DNA-alpha nonbonded file ships with the runtime because GROMACS requires it,
but it is not copied to the two-file user output directory. Common test-only
force-field files are under `tests/data/forcefield` and are not installed.

## Validate

Run all automated checks:

```bash
python -m pytest -v
```

Run only the GROMACS matrix, optionally selecting a non-default executable:

```bash
NUCLEIC_BUILDER_AMBERCLASSIC_HOME=/path/to/AmberClassic-v2.0 \
GMX=/path/to/gmx python -m pytest -m integration -vv
```

The integration test clearly skips when `gmx` is unavailable. It exercises the
upstream dsRNA, independent 1RNA duplex, experimental 3G9Y single strand,
experimental 6YMC hairpin, sequence-derived ideal A-RNA duplex, the authors'
unchanged dsDNA fixture, and a sequence-derived ideal B-DNA duplex. Each case is
solvated, neutralized with 0.15 M NaCl,
preprocessed without `-maxwarn`, minimized with position restraints, run for
1,000 restrained MD steps, and continued for 1,000 unrestrained steps with no
`-DPOSRES`. The preprocessed topologies are checked to prove that the first MD
stage contains `[ position_restraints ]` and the second does not. Both stages
must finish without GROMACS warnings, LINCS errors, NaNs, or fatal errors. The
matrix is dsRNA/`legacy`, 1RNA/`off`, ssRNA/`off`,
hairpin/`intrachain`, sequence-derived dsRNA/`off`, dsDNA/`legacy`, and
sequence-derived dsDNA/`off`. Sequence tests skip
clearly unless `NUCLEIC_BUILDER_AMBERCLASSIC_HOME` points to the pinned backend.
The single-rank setting avoids domain-decomposition
limits from long elastic bonds; it is not warning suppression.

The unchanged upstream RNA regression and DNA-alpha smoke test can be run directly:

```bash
python -m pytest \
  nucleic_builder/_vendor/martini_3_dna_rna/tests/test_1_dna.py \
  nucleic_builder/_vendor/martini_3_dna_rna/tests/test_2_rna.py -v
```

The RNA suite compares the file-based converter against the authors' hard-coded
reference implementation, coordinates, and every emitted ITP section. The DNA
smoke test invokes the authors' separate alpha converter and requires non-empty
PDB/ITP output with atom and topology records.

Build and inspect the wheel with:

```bash
python -m pip install build
python -m build --wheel
unzip -l dist/nucleic_builder-*.whl
```

The wheel contains the two distinctly named runtime converters and exactly one
required copy of each runtime RNA and DNA-alpha map/parameter set. Vendored
regression tests, their reference converter, test PDBs, symlinked parameter
directory, and duplicate reference parameter files remain in the source
checkout but are excluded from installation.

## Examples

| Structural class | Input | Generated output | Validation policy |
|---|---|---|---|
| Duplex | [`examples/1RNA.pdb`](examples/1RNA.pdb) | [`examples/generated`](examples/generated) | `off` |
| Single strand | [`examples/ssrna/3G9Y-chain-C-alt-A.pdb`](examples/ssrna/3G9Y-chain-C-alt-A.pdb) | [`examples/ssrna/generated`](examples/ssrna/generated) | `off` |
| Hairpin/stem-loop | [`examples/hairpin/6YMC-chain-A.pdb`](examples/hairpin/6YMC-chain-A.pdb) | [`examples/hairpin/generated`](examples/hairpin/generated) | `intrachain` |
| Ideal A-form duplex | `--sequence GCAUCG` | [`examples/sequence/generated`](examples/sequence/generated) | `off` |
| **Experimental ideal B-DNA duplex** | `--polymer dna --sequence GCATCG` | [`examples/dna/generated`](examples/dna/generated) | `off` (execution test only) |

The 3G9Y RNA was crystallized while bound to a ZRANB2 zinc-finger protein; the
checked-in input contains its experimentally determined six-nucleotide AGGUAA
chain C with alternate conformers resolved to conformer A. The checked-in 6YMC
input contains the 26-mer stem-loop from a crystal with two RNA copies, barium,
and water. These inputs validate selection, conversion, topology
consistency, and short-run numerical execution. **They do not validate the
free-solution behavior of an isolated ssRNA or hairpin.** Both checked-in
inputs retain the selected experimental coordinates; no structure was
predicted or repaired. Full source, selection, and checksum provenance is in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Attribution and license

The RNA conversion implementation, experimental DNA-alpha implementation, and
their separate parameter directories come from
[DanYev/Martini-3-DNA-RNA](https://github.com/DanYev/Martini-3-DNA-RNA) and are
used under the upstream README's GNU GPL v3 declaration. The complete GPLv3
text is in [LICENSE](LICENSE), and the combined project is distributed as
`GPL-3.0-only`.

Please cite:

> D. Yangaliev and S. B. Ozkan, “Coarse-grained RNA model for the Martini 3
> force field,” *Biophysical Journal* 125(2), 445–456 (2026),
> https://doi.org/10.1016/j.bpj.2025.07.034.

That publication validates RNA. It must not be cited as a published validation
of the DNA-alpha model; upstream says the DNA publication is still pending.

Copied files, historical copyright, test-data licensing, water-box source, and
PDB provenance are documented in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Known limitations

- RNA sequence construction supports only an ideal canonical A-form duplex; experimental DNA sequence construction supports only an ideal canonical B-form duplex. Neither predicts folds or represents a free-solution ensemble.
- The sequence backend requires the exact local AmberClassic checkout; the registration-gated local 3DNA distribution is not used and no remote API is called.
- RNA accepts only A/C/G/U. Experimental DNA accepts only A/C/G/T and uses the separately pinned unpublished alpha model. RNA/DNA hybrids and all modified nucleotides are rejected.
- Only one PDB model is accepted; alternate conformers must be resolved first.
- Nucleic-acid chain segments are represented in one GROMACS molecule type, matching the authors' converters.
- GRO cannot store chain labels; chain breaks remain in topology connectivity and provenance comments.
- Elastic-network numerical parameters are fixed to the authors' values; only network scope can be selected.
- An elastic network preserves selected starting distances and can bias dynamics; a smoke test is not scientific validation of a policy for every RNA or DNA system.
- The paper specifically limits use for hybridization, melting, and intercalation and calls for further refinement for ssRNA and complex secondary structures.
- The extracted 3G9Y and 6YMC coordinates came from protein-bound and barium-containing crystals, respectively, not validated free-solution ensembles.
- Solvation and ionization remain validation-only, not user-facing output features.
- No automated structure repair, secondary-structure detection, or class inference is provided.
- The DNA-alpha force field is unpublished and experimental; the smoke tests establish GROMACS execution, not physical accuracy or production suitability.
