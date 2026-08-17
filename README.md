# Nucleic Builder — Martini 2 and Martini 3 RNA/DNA

`nucleic_builder` converts a canonical RNA or
DNA sequence—or an existing canonical all-atom PDB—into exactly two
molecule-level GROMACS files:

- `NAME.itp` — molecule topology
- `NAME.gro` — coarse-grained coordinates in nm

## Current deliverable

`--martini-version {2,3}` selects the backend. Martini 3 remains the default,
so every existing command that omits the option retains its current Martini 3
behavior and output.

| Input | Structure used | Martini model | Output |
|---|---|---|---|
| RNA sequence, `A/C/G/U`, duplex or single | Ideal A-form duplex or one A-form-derived strand | Martini 2 RNA or published Martini 3 RNA | `NAME.itp`, `NAME.gro` |
| DNA sequence, `A/C/G/T` | Ideal B-form duplex | Martini 2 DNA or **experimental Martini 3 DNA-alpha** | `NAME.itp`, `NAME.gro` |
| Canonical RNA PDB, one or two chains | Supplied geometry | Martini 2 RNA or published Martini 3 RNA | `NAME.itp`, `NAME.gro` |
| Canonical DNA PDB, one or two chains | Supplied geometry | Martini 2 DNA or **experimental Martini 3 DNA-alpha** | `NAME.itp`, `NAME.gro` |

RNA sequence input supports an explicit strand mode. `duplex` preserves the
original ideal A-form duplex behavior and remains the default; `single` keeps
only the entered 5′→3′ strand in an idealized conformation derived from that
A-form geometry. Neither mode predicts whether a sequence folds, hybridizes,
or remains stable in solution. In particular, single-strand mode does **not**
predict the RNA's folded solution structure. Existing-PDB input retains the
supplied single-strand, duplex, or hairpin geometry. Modified nucleotides,
RNA/DNA hybrids, user-facing solvation, and structure prediction are
deliberately out of scope.

The program does not recreate or fit force-field parameters. The default
Martini 3 RNA route remains a
strict wrapper around the implementation published by Danis Yangaliev and
S. Banu Ozkan. Its behavior and pin remain unchanged at
[`e761b7349fdf61dd485053c000dbb642f24ff9d8`](https://github.com/DanYev/Martini-3-DNA-RNA/commit/e761b7349fdf61dd485053c000dbb642f24ff9d8).
DNA is an isolated wrapper around that snapshot's separately named
`martinize_dna_alpha.py` and `dna_alpha_itps/`, independently recorded at the
same exact commit. DNA never enters the RNA converter or RNA parameter set.

The importable Python package is named `nucleic_builder`; its public API covers
both versioned model paths.

## Martini 2 backend

Martini 2 uses the official Martini Force Field Initiative DNA/RNA release,
`na-tutorials_20170815.tar`. The archive is pinned at SHA-256
`15ba5bf45b9890603f0113d2021074f397a7f8d0264cb2093e970198f4b6c20b`.
Its converter identifies itself as version 2.2; the standard-water force-field
file retains its upstream name `martini_v2.1-dna.itp`. That file is the combined
Martini 2 DNA/RNA particle definition file—its name is not evidence that RNA
beads were translated to DNA types. The backend preserves the official native
types, including `Q0`, `SN0`, `SNda`, and the `T*` base beads.

The official converter is Python 2 code. The distribution therefore contains
both the byte-exact upstream source and a Python 3 compatibility port. The port
changes Python iteration, slicing, and printing mechanics only; its mappings,
bead names/types, bonded parameters, and topology-mode settings are unchanged.
The upstream source checksum is
`e02a0ede1f444ccbd7fc9a7e2c0ee6910642887490210bf2f7c1076a2cce3edb`;
the compatibility port checksum is
`ee858476b4e09e0f13d0131ed9b2f617792ca87142af8364daa453913eb8e9fd`.
The converter's own emitted topology header calls the model a development beta
and says not to use it for production runs. That upstream warning is preserved
verbatim in every Martini 2 ITP and summarized in its provenance header.

Supported Martini 2 structures are:

- RNA sequence input in `single` or `duplex` mode.
- DNA sequence input as the pipeline's existing B-form duplex.
- Canonical RNA or DNA PDB input containing one chain or exactly two chains.
- Single-chain hairpin geometry supplied by PDB; no folding is predicted.

PDBs with more than two nucleic-acid chains are rejected. The released
converter says its merged elastic-network implementation is hard-coded for two
chains, so silently generalizing it would not be the official backend.

The existing elastic policy has version-specific meaning for Martini 2:

| Policy | One-chain Martini 2 mode | Two-chain Martini 2 mode |
|---|---|---|
| `legacy` | Official `ss-stiff` | Official `ds-stiff`, including cross-chain network bonds |
| `intrachain` | Official `ss-stiff` | Official `ds-stiff`, with only cross-chain function-6 bonds removed |
| `off` | Official `ss` | Official `ds-stiff` bonded model, with all function-6 network bonds removed |

The stiff official network uses a 1.0 nm upper cutoff and
500 kJ mol⁻¹ nm⁻² force constant. The official soft duplex mode is not exposed
by this three-policy interface. A duplex with `off` is mechanically generatable
but is not scientifically recommended: the DNA and RNA papers state that the
unrestrained prehybridized duplex is not reliably stable.

The exact installed Martini 2 files are under
`nucleic_builder/_vendor/martini_2_nucleic/`:

- `martini_v2.1-dna.itp`, SHA-256
  `cc7c200dff400e97311213b93127697c6f8c21edb2350926072f0194eb90efe6`
- `martini_v2.1P-dna.itp`, SHA-256
  `b8dea4ffbef3a439db0baa465825528db919b026b1a821b374f1a6a605912ff0`
- `martini_v2.0_ions.itp`, SHA-256
  `c5b9b5b9541aa6d77b5b41a4b19dee62c1b8631c73e79cf76e1b281a144b4b4e`

Use the standard file unless intentionally following the polarizable-water
protocol. Do not include the standard and polarizable files together.

Please cite the model matching the polymer:

- DNA: J. J. Uusitalo et al., *J. Chem. Theory Comput.* 11, 3932–3945
  (2015), https://doi.org/10.1021/acs.jctc.5b00286.
- RNA: J. J. Uusitalo et al., *Biophys. J.* 113, 246–256 (2017),
  https://doi.org/10.1016/j.bpj.2017.05.043.

The Martini 2 papers explicitly rule out folding/hairpin formation,
hybridization, melting, and intercalation because base pairing lacks directional
hydrogen bonding and duplex/tertiary structures rely on elastic networks. The
RNA paper recommends at most a 10 fs timestep and notes that smaller steps may
be needed. Martini 2 is now legacy material on the official site. These limits,
the exact source, version, citations, and checksums are written into every
Martini 2 ITP provenance header.

## Martini 3 experimental DNA-alpha warning

> **Martini 3 DNA output is EXPERIMENTAL / UNPUBLISHED DNA-ALPHA.** The authors' upstream
> README says the DNA force field is alpha and that a publication is being
> prepared. The published Martini 3 paper cited below is an RNA paper, not DNA
> validation. Successful conversion and a short GROMACS smoke test establish
> file compatibility and numerical execution only—not scientific accuracy,
> production readiness, duplex stability, melting behavior, or transferability.

DNA must be selected explicitly with `--polymer dna`. For the default Martini 3
backend, the warning is written into the ITP provenance, GRO title, and every
DNA CLI run. Martini 2 DNA is the separately published Uusitalo model described
above, not this DNA-alpha backend.

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
duplexed in solution. Likewise, extracting one strand from that geometry is
only a reproducible idealized starting conformation; it is not a predicted
free-solution fold or ensemble.

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
AmberClassic NAB `fd_helix("arna", sequence, "rna")`. Duplex mode retains the
two verified NAB strands. Single mode verifies the same duplex and retains
only its entered strand before passing it through the normal strict RNA PDB
validation and Martini conversion path. DNA uses the same pinned NAB
installation with its documented right-handed Arnott B-DNA setting,
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

Every command requires an explicit elastic-network choice. The examples in the
next four subsections omit `--martini-version` and therefore run Martini 3.
For an ordinary Martini 3 sequence-derived duplex, use `off` unless a specific
validated protocol calls for a different policy.

### RNA sequence to ideal A-form duplex

```bash
python -m nucleic_builder \
  --sequence GCAUCG \
  --strand-mode duplex \
  --name SEQRNA \
  --output-dir output \
  --elastic-network off
```

This means strand A is `5′-GCAUCG-3′`. The generated strand B is
`5′-CGAUGC-3′`, paired antiparallel as `3′-CGUAGC-5′`.
Omitting `--strand-mode` preserves this duplex behavior for compatibility.

### RNA sequence to one idealized single strand

```bash
python -m nucleic_builder \
  --sequence GCAUCG \
  --strand-mode single \
  --name SEQSSRNA \
  --output-dir output \
  --elastic-network off
```

This produces exactly one `5′-GCAUCG-3′` chain. Its coordinates are obtained
by retaining strand A from the verified ideal NAB A-form duplex. This gives a
deterministic, extended helical starting conformation; it does not predict the
RNA's secondary structure, tertiary fold, or free-solution ensemble.

For the 475-nt BNT162b2 fragment used by Grzetic et al.:

```bash
python -m nucleic_builder \
  --sequence GAGAAUAAACUAGUAUUCUUCUGGUCCCCACAGACUCAGAGAGAACCCGCCACCAUGUUCGUGUUCCUGGUGCUGCUGCCUCUGGUGUCCAGCCAGUGUGUGAACCUGACCACCAGAACACAGCUGCCUCCAGCCUACACCAACAGCUUUACCAGAGGCGUGUACUACCCCGACAAGGUGUUCAGAUCCAGCGUGCUGCACUCUACCCAGGACCUGUUCCUGCCUUUCUUCAGCAACGUGACCUGGUUCCACGCCAUCCACGUGUCCGGCACCAAUGGCACCAAGAGAUUCGACAACCCCGUGCUGCCCUUCAACGACGGGGUGUACUUUGCCAGCACCGAGAAGUCCAACAUCAUCAGAGGCUGGAUCUUCGGCACCACACUGGACAGCAAGACCCAGAGCCUGCUGAUCGUGAACAACGCCACCAACGUGGUCAUCAAAGUGUGCGAGUUCCAGUUCUGCAACGACCCCUUCC \
  --strand-mode single \
  --name BNT162B2_475 \
  --output-dir output \
  --elastic-network off
```

This command creates a 475-residue, one-chain RNA molecule. It does not add a
complementary strand or any lipid-nanoparticle components.

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
paired as `3′-CGTAGC-5′`. DNA must be selected explicitly. With the default
Martini 3 backend, every DNA run and generated file is marked
experimental/unpublished.

### Explicit Martini versions

The equivalent explicit Martini 3 RNA command is:

```bash
python -m nucleic_builder \
  --martini-version 3 \
  --sequence GCAUCG \
  --strand-mode duplex \
  --name M3RNA \
  --output-dir output \
  --elastic-network off
```

Use the published Martini 2 RNA backend and its official stiff duplex network
with:

```bash
python -m nucleic_builder \
  --martini-version 2 \
  --sequence GCAUCG \
  --strand-mode duplex \
  --name M2RNA \
  --output-dir output \
  --elastic-network legacy
```

For Martini 2 DNA, add both the version and polymer selections:

```bash
python -m nucleic_builder \
  --martini-version 2 \
  --polymer dna \
  --sequence GCATCG \
  --name M2DNA \
  --output-dir output \
  --elastic-network legacy
```

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

| Value | Martini 3 behavior |
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

The RNA sequence route accepts one non-empty strand containing only `A/C/G/U`.
`--strand-mode duplex` constructs its antiparallel complement;
`--strand-mode single` retains exactly the entered strand. It is not a
sequence-file parser: no FASTA header, whitespace cleanup, DNA spelling,
modified bases, or folding is implied. `--strand-mode` is rejected for PDB and
DNA input; DNA sequence behavior remains an ideal B-form duplex.

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

For duplex sequence input, provenance records both strands in their 5′→3′
forms, the paired 3′→5′ orientation, and the existing ideal-dsRNA input mode.
For single sequence input, it explicitly records `Strand mode: single`, the
one retained sequence, the A-form-derived extraction, and that no folded
solution structure was predicted. Both modes record the AmberClassic
repository/tag/commit and NAB settings, private all-atom PDB SHA-256, and
unchanged Martini upstream provenance. Temporary NAB source, executable,
atomic PDB, and CG PDB files are deleted; only `NAME.itp` and `NAME.gro` are
published.

Every DNA ITP additionally records the separate DNA-alpha converter, parameter
directory, independent upstream pin, and experimental/unpublished status. The
GRO title carries the same warning. RNA and DNA provenance identify their
polymer and cannot silently share a converter.

## Use the molecule files in GROMACS

The generated ITP contains one molecule type, not a complete system topology.

### Martini 3 topology

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

### Martini 2 topology

Use the exact standard-water file packaged with this project, followed by ions
if the system contains them, and then the generated molecule:

```gromacs
#include "martini_v2.1-dna.itp"
#include "martini_v2.0_ions.itp"
#include "M2RNA.itp"

[ system ]
Martini 2 RNA system

[ molecules ]
M2RNA  1
```

Despite its historical filename, `martini_v2.1-dna.itp` is the exact combined
DNA/RNA file shipped by the official 2017 release. It must not be replaced with
a Martini 3 file or a file whose bead types have been renamed. The optional
`martini_v2.1P-dna.itp` is a complete alternative for the polarizable-water
protocol; do not include it together with the standard file.

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
unchanged dsDNA fixture, a sequence-derived idealized RNA single strand, and a
sequence-derived ideal B-DNA duplex. Each case is
solvated, neutralized with 0.15 M NaCl,
preprocessed without `-maxwarn`, minimized with position restraints, run for
1,000 restrained MD steps, and continued for 1,000 unrestrained steps with no
`-DPOSRES`. The preprocessed topologies are checked to prove that the first MD
stage contains `[ position_restraints ]` and the second does not. Both stages
must finish without GROMACS warnings, LINCS errors, NaNs, or fatal errors. The
matrix is dsRNA/`legacy`, 1RNA/`off`, ssRNA/`off`,
hairpin/`intrachain`, sequence-derived dsRNA/`off`, sequence-derived
ssRNA/`off`, dsDNA/`legacy`, and
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

The wheel contains the two distinctly named Martini 3 runtime converters, the
Martini 2 compatibility converter and byte-exact upstream source, the three
official Martini 2 parameter files, and exactly one required copy of each
runtime Martini 3 RNA and DNA-alpha map/parameter set. Vendored
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
| Idealized A-form-derived single strand | `--sequence GCAUCG --strand-mode single` | `NAME.itp`, `NAME.gro` | `off` |
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

The Martini 2 backend and parameter files come from the official Martini Force
Field Initiative `na-tutorials_20170815.tar` download. That archive contains
copyright notices in its force-field files but no standalone license and no
explicit license statement in its README or converter; this absence is
recorded in the notices rather than being replaced with an invented license.

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

- RNA sequence construction supports an ideal canonical A-form duplex and an idealized single strand extracted from that geometry; experimental DNA sequence construction supports only an ideal canonical B-form duplex. None predicts folds or represents a free-solution ensemble.
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
- Martini 2 input is limited to one or two chains because the official converter's merged-network implementation is hard-coded for two strands.
- Martini 2 cannot model folding/hairpin formation, hybridization, melting, or intercalation; its duplex and tertiary-structure use relies on elastic networks, and RNA generally requires a timestep of 10 fs or less.
