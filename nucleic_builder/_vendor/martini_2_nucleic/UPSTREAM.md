# Official Martini 2 DNA/RNA release provenance

- Source: Martini Force Field Initiative, `na-tutorials_20170815.tar`
- Source page: https://cgmartini.nl/docs/downloads/force-field-parameters/martini2/nucleic_acids.html
- Archive URL: https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini2/nucleic_acids/rna/na-tutorials_20170815.tar
- Archive SHA-256: `15ba5bf45b9890603f0113d2021074f397a7f8d0264cb2093e970198f4b6c20b`
- Release README date: 2017-06-05
- Retrieved: 2026-08-12

Vendored byte-for-byte from the archive:

- `martinize-nucleotide.py.upstream`: `e02a0ede1f444ccbd7fc9a7e2c0ee6910642887490210bf2f7c1076a2cce3edb`
- `martini_v2.1-dna.itp`: `cc7c200dff400e97311213b93127697c6f8c21edb2350926072f0194eb90efe6`
- `martini_v2.1P-dna.itp`: `b8dea4ffbef3a439db0baa465825528db919b026b1a821b374f1a6a605912ff0`
- `martini_v2.0_ions.itp`: `c5b9b5b9541aa6d77b5b41a4b19dee62c1b8631c73e79cf76e1b281a144b4b4e`
- `README.upstream`: `be27f9805bece1d02eb599797b8b0b0a32255688e79b93b631c47ad9382a3002`

`martinize_nucleotide_py3.py` is a Python 3 compatibility port of the exact
upstream converter. It was produced with Python's `2to3`, then received only
the compatibility fixes that `2to3` missed: Python 3 slice dispatch,
subscriptable `zip` results, and iterator `next()` calls. Its force-field
classes, mappings, bead names/types, bonded parameters, and topology-mode
settings are unchanged. The port SHA-256 is
`ee858476b4e09e0f13d0131ed9b2f617792ca87142af8364daa453913eb8e9fd`.

The downloaded archive contains copyright lines in the force-field files but
no standalone `LICENSE` or `COPYING` file and no explicit license statement in
its README or converter. This repository records that absence rather than
inventing a license for those upstream files. The surrounding
`nucleic_builder` project remains GPL-3.0-only.

Scientific citations:

- J. J. Uusitalo, H. I. Ingolfsson, P. Akhshi, D. P. Tieleman, and
  S. J. Marrink, “Martini Coarse-Grained Force Field: Extension to DNA,”
  *J. Chem. Theory Comput.* 11, 3932–3945 (2015),
  https://doi.org/10.1021/acs.jctc.5b00286.
- J. J. Uusitalo, H. I. Ingolfsson, S. J. Marrink, and I. Faustino,
  “Martini Coarse-Grained Force Field: Extension to RNA,”
  *Biophys. J.* 113, 246–256 (2017),
  https://doi.org/10.1016/j.bpj.2017.05.043.
