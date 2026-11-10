# Findings 013 — two independent 30×6 formula matrices, both LLL inputs are non-saturated

## 1. Matrix-version comparison

The repository has multiple matrices that share the same six-column idea but are not the same data:

- draft/stress_kernel.py: 30×6, basis columns use sqrt(2), sqrt(3), pi, lnK, lnN, and N with explicit N exponent −1/3.
- draft/stress_kernel_testing.py: another 30×6 copy of the same older exponent table.
- physics_constants/constant_matrix_analysis.py: first 30 rows form a newer 30×6 table in which the last primitive is N^(1/3), so its corresponding exponents are −1 rather than −1/3.
- physics_constants/constant_matrix_analysis.py then extends to 69 active rows and 69 active names; K_star_meson is commented out in both matrix and names.
- matrix_basis/matrix_analysis.py has 70 rows and 70 names; compared with the 69-row version, it includes the K_star_meson row and shifts the final neutron-lifetime / connection-energy rows.

After converting the last-column convention, the first 30 rows of the old draft matrix and the newer constant_matrix matrix differ in 24 of 30 rows. This is not just a harmless basis normalization: many coefficients of sqrt(2) and sqrt(3) are missing in the older draft. For example:

- hbar requires exponents (−2,−2,0,0,3,−1) because 1/K=1/6=(sqrt2)^−2(sqrt3)^−2; the old draft has (0,0,0,0,3,−1).
- c requires (−4,−4,1,−1,4,0) because K^−2 contributes both prime factors; the old draft has (0,0,1,−1,4,0).
- mP requires (−2,2,−1,0,−3,0) because K/4=3/2; the old draft has (0,2,−1,0,−3,0).

Therefore the saturated-kernel result in FINDINGS_012 is a correct mathematical audit of the old draft matrix, but that old matrix is not a faithful encoding of many displayed formulae. It must not be treated as the canonical physical exponent matrix.

## 2. Saturation test on the newer first-30-row matrix

The first 30 rows of physics_constants/constant_matrix_analysis.py were separately transcribed and audited. This matrix has rank 6 and rational nullity 24.

Applying the same exact method to A=2Mᵀ gives Smith diagonal (1,2,2,2,2,2). The source's row-wise denominator-clearing procedure again generates a proper sublattice, now of index:

    18014398509481984 = 2^54.

LLL with δ=3/4 on the saturated integer kernel gives 24 exact basis vectors, all with maximum absolute coefficient ≤2. The source-style row-wise basis is incomplete by the index above. Thus the saturation problem is not confined to the older draft; the newer get_kernel_vectors() implementation repeats it.

Reproducible artifact: audit_integer_kernel_corrected_matrix.py.

## 3. Interpretation

The 24 short vectors are exact exponent relations among the 30 rows as encoded in this matrix. They are not 24 independent experimental confirmations. Their physical meaning still depends on validating every exponent row against the underlying formula, units, and reference data.

## 4. Current verdict

- Old 30×6 matrix: 24/30 rows differ from the newer first-30 matrix after matching the N primitive convention.
- Old matrix integer-kernel basis: non-saturated, index ≈6.356×10^26.
- Newer first-30 matrix integer-kernel basis: non-saturated, index 2^54.
- Saturated LLL for newer first-30 matrix: 24 exact vectors, all max coefficient ≤2.
- Repository matrix versions: 30×6, 69×6, and 70×6; the 69/70 distinction is explained by the commented/active K_star_meson row, not by a names-length bug.
- The article's 55×6 / 29-vector result: still not reproduced.