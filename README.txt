Strict State Hierarchies for Graph Homomorphism Polynomials
Yue Wang -- finite computational supplement

Requirements: Python 3.9 or later; standard library only. No network access.
From the extracted directory run:
    python verify.py
    python verify_manifest.py

verify.py generates the graphs, evaluates them by generic factor elimination,
compares conditional-core evaluations, checks the exact successive-gap formula,
and validates the archived edge lists and full integer values. It writes
results.json deterministically. The manifest check verifies the shipped files.

The general theorems have mathematical proofs in the manuscript. An unequal
modular evaluation certifies unequal integer values. Equal modular evaluations
alone do not certify a polynomial identity. In particular, the Tutte identities
are established in the manuscript, not inferred from samples.

FILES
verify.py: graph generation, both evaluation algorithms and all finite checks.
results.json: deterministic output, with residues modulo 1000003.
graphs.json: explicit zero-based edge lists and bipartitions for 13 graphs.
exact_values.json: full integers and input parameters for the three-member family.
SHA256SUMS.txt and verify_manifest.py: integrity checking.

NOTATION
F_rR_tT_dD has R terminals, T+D even cores and T-D odd cores.
The first R vertices are the terminals. In the paper, X_r=F_r,1,1 and
Y_r=F_r,1,0. The prefix B1_ denotes one added vertex on each bipartition side.
Source graphs are simple, unweighted, uncolored graphs; the bipartition metadata
is included for interpreting the construction and is not part of graph equality.
The target matrix is symmetric; its diagonal weights are allowed.

PARAMETERS
Non-Potts table inputs use unit activities and, with zero-based state indices,
A[i][j]=(min(i,j)+2)*(max(i,j)+3)+3*(i==j)+1.
Potts table inputs use unit activities, diagonal entries 5 and off-diagonal 1.
The integer values in exact_values.json use these same non-Potts parameters.

PROVENANCE
The two evaluation algorithms are implemented in one self-contained script.
They provide algorithmic cross-checks, not an independent third-party review.
This supplement was prepared for the present manuscript from its graph definitions.
