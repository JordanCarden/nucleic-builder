# Third-party notices and provenance

## Official Martini 2 DNA/RNA release

- Project: Martini Force Field Initiative DNA/RNA model
- Source page: https://cgmartini.nl/docs/downloads/force-field-parameters/martini2/nucleic_acids.html
- Release archive: `na-tutorials_20170815.tar`
- Archive URL: https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini2/nucleic_acids/rna/na-tutorials_20170815.tar
- Archive SHA-256: `15ba5bf45b9890603f0113d2021074f397a7f8d0264cb2093e970198f4b6c20b`
- Release README date: 2017-06-05
- Retrieved: 2026-08-12
- Converter-reported version: 2.2

Vendored byte-for-byte under
`nucleic_builder/_vendor/martini_2_nucleic/`:

- `martinize-nucleotide.py.upstream`, SHA-256
  `e02a0ede1f444ccbd7fc9a7e2c0ee6910642887490210bf2f7c1076a2cce3edb`
- `martini_v2.1-dna.itp`, SHA-256
  `cc7c200dff400e97311213b93127697c6f8c21edb2350926072f0194eb90efe6`
- `martini_v2.1P-dna.itp`, SHA-256
  `b8dea4ffbef3a439db0baa465825528db919b026b1a821b374f1a6a605912ff0`
- `martini_v2.0_ions.itp`, SHA-256
  `c5b9b5b9541aa6d77b5b41a4b19dee62c1b8631c73e79cf76e1b281a144b4b4e`
- `README.upstream`, SHA-256
  `be27f9805bece1d02eb599797b8b0b0a32255688e79b93b631c47ad9382a3002`

The standard-water `martini_v2.1-dna.itp` is the release's combined DNA/RNA
particle definition file. Its historical upstream filename is retained. The
backend does not translate or rename its native nucleic-acid bead types.

The official converter uses Python 2 syntax. `martinize_nucleotide_py3.py` is a
compatibility port made with `2to3` plus the missed Python 3 iterator, slicing,
and `zip` fixes. Its SHA-256 is
`ee858476b4e09e0f13d0131ed9b2f617792ca87142af8364daa453913eb8e9fd`.
The force-field classes, mappings, bead names/types, bonded parameters, and
topology-mode settings are unchanged. Wrapper post-processing changes only the
requested molecule name, normalizes the converter's gapped one-bead-per-charge-
group numbering, replaces ephemeral private paths in the converter-options
comment, scopes its `RUBBER_BANDS` preprocessor symbol, selects network scope,
and adds provenance.

The downloaded archive has copyright lines in the force-field files but no
standalone `LICENSE` or `COPYING` file and no explicit license statement in its
README or converter. This notice records that fact and does not infer a license
for those upstream files. The surrounding project remains GPL-3.0-only.

Scientific citations:

> J. J. Uusitalo, H. I. Ingolfsson, P. Akhshi, D. P. Tieleman, and
> S. J. Marrink, “Martini Coarse-Grained Force Field: Extension to DNA,”
> *J. Chem. Theory Comput.* 11, 3932–3945 (2015), DOI
> 10.1021/acs.jctc.5b00286.

> J. J. Uusitalo, H. I. Ingolfsson, S. J. Marrink, and I. Faustino,
> “Martini Coarse-Grained Force Field: Extension to RNA,” *Biophys. J.* 113,
> 246–256 (2017), DOI 10.1016/j.bpj.2017.05.043.

The publications state that the models do not support folding/hairpin
formation, hybridization, melting, or intercalation because base pairing lacks
directional hydrogen bonding and duplex/tertiary structures depend on elastic
networks. The RNA publication recommends a maximum 10 fs timestep and notes
that smaller steps may be necessary. The official website now categorizes the
Martini 2 tutorial as legacy material. The official converter additionally
emits a topology header that calls the model a development beta and says not
to use it for production runs; generated ITPs preserve that warning.

## Martini 3 RNA converter and RNA parameters

- Project: `DanYev/Martini-3-DNA-RNA`
- Source: https://github.com/DanYev/Martini-3-DNA-RNA
- Exact snapshot: `e761b7349fdf61dd485053c000dbb642f24ff9d8`
- Commit date: 2026-07-21
- Retrieved: 2026-07-27
- Historical upstream notice: `Copyright (c) 2025, DY`
- License declared by the upstream README: GNU General Public License v3.0

The inspected snapshot contains no `LICENSE` or `COPYING` file and no source
license header. This project conservatively treats the declaration as
`GPL-3.0-only`, includes the complete GPLv3 text at the repository root, and
retains the upstream source, file names, docstring citation, regression test,
and reference implementation.

Vendored runtime files under `nucleic_builder/_vendor/martini_3_dna_rna`:

- `martinize_rna_v3.0.0.py`
- `rna_v3.0.0_itps/rna_bb.{itp,map}`
- `rna_v3.0.0_itps/rna_{A,C,G,U}.{itp,map}`
- `rna_v3.0.0_itps/martini_v3.0.0_rna.itp`
- the upstream README

Vendored regression assets, unchanged:

- `tests/dsRNA.pdb`
- `tests/martinize_rna_ref.py`
- `tests/test_2_rna.py`
- `tests/martinize_ref_itps/` and its upstream symlink

The regression assets and their duplicate reference parameters are retained in
the source tree to preserve the authors' test unchanged. They are explicitly
excluded from the installed wheel; the wheel contains only the converter and
the single runtime parameter/map directory.

No vendored upstream source or RNA parameter file was locally modified. New
code in `nucleic_builder` wraps it with stricter validation, temporary residue-name
normalization, deterministic provenance, GRO output, and consistency checks.

Scientific citation:

> D. Yangaliev and S. B. Ozkan, “Coarse-grained RNA model for the Martini 3
> force field,” *Biophysical Journal* 125(2), 445–456 (2026), DOI
> 10.1016/j.bpj.2025.07.034.

The current upstream README accidentally says “DNA model” in its citation and
links a Zenodo record for a different project. This project uses the title and
DOI from the journal record and does not attribute the pinned snapshot to that
Zenodo record.

## Experimental DNA-alpha converter and parameters

- Project: `DanYev/Martini-3-DNA-RNA`
- Source: https://github.com/DanYev/Martini-3-DNA-RNA
- Separately recorded DNA-alpha snapshot: `e761b7349fdf61dd485053c000dbb642f24ff9d8`
- Commit date: 2026-07-21
- DNA-alpha files vendored: 2026-07-29
- License declared by the upstream README: GNU General Public License v3.0
- Upstream model status: **alpha; publication being prepared**

The DNA-alpha converter is separately named `martinize_dna_alpha.py`; its maps,
bonded parameters, and nonbonded definitions are under `dna_alpha_itps/`. They
are copied byte-for-byte and loaded only by `nucleic_builder.dna_builder`. They are
never substituted for, or routed through, the published RNA converter and RNA
parameter directory.

Vendored runtime DNA-alpha files:

- `martinize_dna_alpha.py`
- `dna_alpha_itps/dna_bb.{itp,map}`
- `dna_alpha_itps/dna_{A,C,G,T}.{itp,map}`
- `dna_alpha_itps/martini_alpha_dna.itp`

Vendored regression assets, unchanged:

- `tests/dsDNA.pdb`
- `tests/test_1_dna.py`

The regression files stay in the source distribution and are excluded from the
installed wheel. The wheel includes the one required runtime DNA-alpha
parameter directory. The full GPLv3 text and historical-license treatment are
the same as described for the exact upstream snapshot above.

The published Yangaliev/Ozkan paper cited for RNA is not presented as a DNA
publication. This project labels every DNA ITP, GRO title, and CLI run
`EXPERIMENTAL / UNPUBLISHED DNA-ALPHA` and treats GROMACS smoke tests as file
and execution validation only.

## AmberClassic NAB sequence-structure backend

- Project: `Amber-MD/AmberClassic`
- Source: https://github.com/Amber-MD/AmberClassic
- Exact tag: `v2.0`
- Exact commit: `bdb3e0dee5b90f2be2950e26cfad1ae5a7440cae`
- License: GNU GPL v2 or (at the user's option) any later version, with
  compatible component licenses described by that project
- Phase 2A call: `fd_helix("arna", sequence, "rna")`
- Phase 2B DNA call: `fd_helix("abdna", sequence, "dna")` (right-handed
  B-DNA, Arnott fiber-diffraction setting)
- PDB writer option: `-wwpdb`

AmberClassic is an external local dependency and is not copied into this
repository or wheel. `scripts/install_amberclassic_nab.sh` clones and builds
the exact official revision; generated provenance records the pin and settings.
The GPL-2.0-or-later terms are compatible with this combined GPLv3-only
distribution.

Local 3DNA v2.3 `fiber` was evaluated first, but its download requires 3DNA
registration and was not suitable for non-interactive reproducible installation
or redistribution. No 3DNA code, templates, credentials, or generated files are
included, and the remote 3DNA API is never used.

## Martini 3 validation force-field files

The following test-only assets were copied byte-for-byte from the pinned
GPLv3 upstream snapshot's `tests/cgmd_params` directory:

- `martini_v3.0.0.itp`
- `martini_v3.0.0_ions_v1.itp`
- `martini_v3.0.0_solvents_v1.itp`

Only `martini_v3.0.0.itp` was also verified byte-for-byte against the core
Martini force-field project at
https://github.com/marrink-lab/martini-forcefields, which is distributed under
the Apache License 2.0. A copy of that license is retained as
`tests/data/forcefield/LICENSE.Apache-2.0`. The copied ion and solvent files
differ from the current official-repository versions, so this project does not
represent them as unchanged Apache upstream files. The RNA-specific nonbonded
file is from the GPLv3 snapshot described above.

The 400-bead Martini water box at `tests/data/water.gro` was downloaded from
the official Martini tutorial URL:

https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/example_applications/solvent_systems/water.gro

SHA-256: `a14626e1b935673f42391c7aa4d51a2ea70fed98a7abfc1b677940f04a005f24`.

## Independent PDB example

`examples/1RNA.pdb` is the unmodified RCSB archive file for PDB ID 1RNA:

- Record: https://www.rcsb.org/structure/1RNA
- Structure DOI: https://doi.org/10.2210/pdb1RNA/pdb
- Download: https://files.rcsb.org/download/1RNA.pdb
- Retrieved: 2026-07-27
- SHA-256: `a6bbc4f9768643cb01a1b52cefdb3bc2138910c3eef8ae0c27e8f9fd1b4ce906`

PDB archive coordinate files are made available under the CC0 1.0 Universal
Public Domain Dedication. Attribution to the structure authors and RCSB PDB is
nevertheless retained here.

## Experimental single-strand PDB example

`examples/ssrna/3G9Y-chain-C-alt-A.pdb` is a deterministic RNA-only extract
from the RCSB archive file for PDB ID 3G9Y:

- Record: https://www.rcsb.org/structure/3G9Y
- Structure DOI: https://doi.org/10.2210/pdb3G9Y/pdb
- Download: https://files.rcsb.org/download/3G9Y.pdb
- Method/resolution: X-ray diffraction, 1.40 Å
- Retrieved: 2026-07-29
- Archive SHA-256: `f0843de3388e95f655587d4c2a544017c4889c84acf509a1a5c6f323eab861e0`
- Derived input SHA-256: `20c3c7ff17b6076f0f9ec6e9905d441ba05e49ffbb6080b050ab9755098fde59`

The selected experimental component is the complete six-nucleotide AGGUAA
ssRNA chain C. It was crystallized bound to the second zinc finger of
ZRANB2/ZNF265; it is not an experimentally determined free-solution ssRNA
ensemble. Residues 5 and 6 have alternate conformers in the archive. The
derivation retains blank-location atoms plus conformer A, clears the selected
altloc labels, and excludes the bound protein, zinc ions, and solvent. This is
an explicit example-data decision; the user-facing converter still rejects
unresolved alternate locations.

## Experimental hairpin PDB example

`examples/hairpin/6YMC-chain-A.pdb` is a deterministic RNA-only extract from
the RCSB archive file for PDB ID 6YMC:

- Record: https://www.rcsb.org/structure/6YMC
- Structure DOI: https://doi.org/10.2210/pdb6YMC/pdb
- Download: https://files.rcsb.org/download/6YMC.pdb
- Method/resolution: X-ray diffraction, 2.00 Å
- Retrieved: 2026-07-29
- Archive SHA-256: `f0d051d0fe2981c7f021c7e14491fd21ee5af1b6ba27fe6a5a23aba2a3587bba`
- Derived input SHA-256: `c4c1181ca5d5056f5b543a672c645ffe52471f6c4866c6ce10aabb2076f550b5`

The archive describes a 26-mer stem-loop and its crystal contains two
crystallographic RNA copies plus barium and water. The selected chain is not an
experimentally validated isolated free-solution ensemble. The derivation
retains canonical RNA chain A only, in its deposited order and coordinates,
and excludes all other records. No alternate conformer selection or structure
modeling was required.

The 3G9Y and 6YMC archive coordinate files are available under the same wwPDB
CC0 terms noted for 1RNA. Attribution and exact source checksums are retained
even though CC0 does not require them.
