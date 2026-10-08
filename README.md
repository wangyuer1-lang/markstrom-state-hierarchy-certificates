# State-hierarchy certificates

Finite computational supplement for **Strict State Hierarchies for Graph
Homomorphism Polynomials**, by Yue Wang.

Version: **v1.0.0**. The repository contains verification code and data.

## Reproduce

Python 3.9 or later; standard library only. After downloading and extracting
the release archive, open a terminal in the extracted directory and run:

```console
python verify_manifest.py
python verify.py
python verify_manifest.py
```

Run Python normally, without `-O` or `-OO`: the checks use assertions.
No network access or external packages are needed for these commands.
Successful execution exits with status zero. `verify.py` prints a JSON report
and deterministically rewrites `results.json`; the final hash check confirms
that the regenerated result matches the archived result.

## What is checked

- The explicit graphs and bipartitions in `graphs.json` (13 graphs).
- Partition-function values from generic factor elimination on adjacency
  lists, compared with conditional gadget evaluation.
- The exact successive-difference formula for a three-member family,
  including the full integer values in `exact_values.json`.
- Selected augmented-graph Potts evaluations and explicit specializations.

Modular computations use the prime 1000003. A nonzero modular difference
certifies unequal integer evaluations. Agreement at selected parameters is
a consistency check, not a proof of a polynomial identity. The general
theorems and Tutte identities are proved in the manuscript.

The two evaluation procedures use different algorithms but share this
implementation. Their agreement is an algorithmic cross-check, not an
independent third-party mathematical review. The explicit edge lists and
target parameters allow an independent implementation to use the same inputs.

## Files

| File | Contents |
| --- | --- |
| `verify.py` | Graph construction and both evaluation procedures |
| `graphs.json` | Zero-based vertex labels, edge lists and bipartitions |
| `exact_values.json` | Target parameters and full integer values |
| `results.json` | Deterministic verification report |
| `README.txt` | Notation and parameter conventions |
| `SHA256SUMS.txt` | SHA-256 manifest for the computational supplement |
| `verify_manifest.py` | Manifest checker |

The SHA-256 manifest covers the original computational supplement; repository
metadata and the CI workflow are outside that manifest.

## Versioned download

[Release v1.0.0](https://github.com/wangyuer1-lang/markstrom-state-hierarchy-certificates/releases/tag/v1.0.0)
contains `state_hierarchy_certificates.zip` and its SHA-256 checksum.
The archive contains the seven supplement files listed above at its root.
The manuscript cites this release rather than the changing default branch.

GitHub Actions runs the checks on Python 3.9, 3.12 and 3.13.
