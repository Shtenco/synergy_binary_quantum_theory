# Findings 007 — exponent matrix dimensions, duplicate rows, and missing LLL

Source: matrix_basis/matrix_analysis.py, blob SHA e58dd577e804b07d07e6c0df1a01041639ad9100.

## 1. The current matrix is 70×6, not 55×6

Direct source parsing finds 70 row vectors in M and 70 names in NAMES. The article reports a 55×6 matrix, rank 6, nullity 49. For the current source file, numerical elimination gives rank 6 and nullity 70−6=64. This is a code/article version mismatch. The advertised 49-dimensional identity space does not describe the current checked-in matrix.

## 2. Duplicate exponent rows

The current matrix has 68 unique row vectors and two duplicate rows:

- row 57 and row 62: [-2, -3, 1/2, 0, 6, -1] (Ksi_0 and Ksi_minus);
- row 56 and row 63: [-2, 2, -2, 0, 6, -1] (Sigma_plus and Sigma_minus).

Thus each pair receives exactly the same formula from the six-dimensional exponent basis. Yet the target dictionary assigns different measured masses to each member of each pair. A single shared formula cannot distinguish the members; any claimed separate prediction requires an additional quantum number/parameter or a different row, neither of which appears in this matrix.

## 3. The repository contains several incompatible matrix versions

The article reports a 55×6 matrix (nullity 49). The current public repository contains at least these distinct versions:

- matrix_analysis.py: 70×6, 70 names, rank 6, nullity 64; this file computes a nullspace and sorts vectors but does not call LLL.
- draft/stress_kernel.py: 30×6, 30 names; it explicitly calls SymPy's .lll() on a 24×30 integer kernel basis, so this is a real LLL implementation for a different, 30-observable matrix.
- physics_constants/constant_matrix_analysis.py: 69 matrix rows but 70 names in the current file, and it does call .lll(); this is a shape/label consistency failure requiring investigation before its outputs can be trusted.
- physics_constants/matrix_analysis.md: reports a 30×6 matrix, nullity 24, and 14 short LLL invariants.

Thus LLL does exist in some repository scripts, but the cited/current matrix_analysis.py does not implement it, and the 30-row LLL result is not the article's claimed 55-row/49-nullity result. The exact 55×6 matrix and exact 29-vector LLL output used for the article have not yet been located in the current source tree.

Sorting a chosen nullspace basis by norm is not LLL reduction. A reproducible LLL claim requires a pinned matrix, exact integer lattice basis, algorithm/implementation version, and output. The repository currently mixes multiple matrices and reports different dimensions/counts.

## 4. Rank 6 does not establish a minimal physical basis

Rank 6 establishes only that the six columns of this chosen exponent matrix are linearly independent over the field used for the rank calculation. It does not establish that the six basis elements are physically fundamental, that the exponent rows were derived rather than selected, or that the nullspace relations are independent experimental predictions.

## 5. Current verdict

- Current checked-in M dimensions: 70×6, not 55×6.
- Numeric rank: 6; corresponding nullity: 64, not 49.
- Duplicate exponent rows: 2 pairs.
- LLL reduction: implemented in a 30×6 draft, absent from the 70×6 matrix_analysis.py; the article's 55×6 / 29-vector result is not yet reproduced.
- Physical interpretation of the six basis columns: NOT ESTABLISHED by rank alone.
- Article/repository version consistency: FAIL.