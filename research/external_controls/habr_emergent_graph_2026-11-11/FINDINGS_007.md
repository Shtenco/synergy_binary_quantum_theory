# Findings 007 — exponent matrix dimensions, duplicate rows, and missing LLL

Source: matrix_basis/matrix_analysis.py, blob SHA e58dd577e804b07d07e6c0df1a01041639ad9100.

## 1. The current matrix is 70×6, not 55×6

Direct source parsing finds 70 row vectors in M and 70 names in NAMES. The article reports a 55×6 matrix, rank 6, nullity 49. For the current source file, numerical elimination gives rank 6 and nullity 70−6=64. This is a code/article version mismatch. The advertised 49-dimensional identity space does not describe the current checked-in matrix.

## 2. Duplicate exponent rows

The current matrix has 68 unique row vectors and two duplicate rows:

- row 57 and row 62: [-2, -3, 1/2, 0, 6, -1] (Ksi_0 and Ksi_minus);
- row 56 and row 63: [-2, 2, -2, 0, 6, -1] (Sigma_plus and Sigma_minus).

Thus each pair receives exactly the same formula from the six-dimensional exponent basis. Yet the target dictionary assigns different measured masses to each member of each pair. A single shared formula cannot distinguish the members; any claimed separate prediction requires an additional quantum number/parameter or a different row, neither of which appears in this matrix.

## 3. The claimed LLL reduction is absent from the checked-in implementation

The source imports SymPy and NumPy. Its get_exact_invariants() calls M_sp.T.nullspace(), clears denominators, divides by the gcd, computes the residual exponent vector, and sorts by Euclidean norm. There is no call to an LLL algorithm or an LLL library in this file. Therefore the article's statement that LLL reduction selected 29 short invariants is not reproduced by this implementation as checked in.

Sorting a chosen nullspace basis by norm is not LLL reduction. Moreover, a nullspace basis is not unique; the short vectors found depend on the basis and the integer-lattice reduction method. A reproducible LLL claim requires an explicit lattice basis, algorithm, implementation/version, and output.

## 4. Rank 6 does not establish a minimal physical basis

Rank 6 establishes only that the six columns of this chosen exponent matrix are linearly independent over the field used for the rank calculation. It does not establish that the six basis elements are physically fundamental, that the exponent rows were derived rather than selected, or that the nullspace relations are independent experimental predictions.

## 5. Current verdict

- Current checked-in M dimensions: 70×6, not 55×6.
- Numeric rank: 6; corresponding nullity: 64, not 49.
- Duplicate exponent rows: 2 pairs.
- LLL reduction in source: NOT IMPLEMENTED.
- Physical interpretation of the six basis columns: NOT ESTABLISHED by rank alone.
- Article/repository version consistency: FAIL.