# Findings 019 — full 69×6 and 70×6 integer kernels independently saturated and LLL-checked

Frozen source revisions:
- `physics_constants/constant_matrix_analysis.py`, SHA `79b5e570ab990c6143ace748881ba20203de049c`.
- `matrix_basis/matrix_analysis.py`, SHA `e58dd577e804b07d07e6c0df1a01041639ad9100`.

Independent exact computation was run on both complete matrices, not only their first 30 rows.

## 1. Matrix dimensions and exact Smith data

| Matrix | Rows × columns | Rank | Integer-kernel rank | SNF diagonal of (2M^T) |
|---|---:|---:|---:|---|
| physics_constants | 69×6 | 6 | 63 | (1,2,2,2,2,2) |
| matrix_basis | 70×6 | 6 | 64 | (1,2,2,2,2,2) |

All rational entries have denominator at most 2. Setting (A=2M^T\in\mathbb Z^{6\times n}), exact Smith normal decomposition gives (D=UAV), with (U,V) unimodular. The final (n-6) columns of (V), transposed, therefore form a saturated basis of the integer kernel of (A). In both cases (|\det V|=1), and direct exact multiplication gives (A B^T=0).

The 70-row matrix is the 69-row matrix plus the K-star row
[
(2,0,-5/2,0,6,-1),
]
inserted before the neutron-lifetime and connection-energy rows. The active row/name counts match: 69/69 and 70/70. The apparent mismatch is caused by commented-out entries and a single-quoted name in the source, not an active-array length bug.

## 2. Saturated LLL results

LLL was applied to the complete saturated integer kernels with δ=3/4. Every output row was multiplied by (A) using exact integer arithmetic and checked to be zero.

| Matrix | Saturated basis vectors | Maximum absolute coefficient | Vectors with all coefficients ≤10 | Exact direct-product checks |
|---|---:|---:|---:|---:|
| 69×6 | 63 | 2 | 63/63 | 63/63 PASS |
| 70×6 | 64 | 2 | 64/64 | 64/64 PASS |

This establishes the integer kernel for each frozen matrix and verifies every returned LLL relation. It does not establish that either matrix correctly represents the physical formulas.

## 3. The source's row-wise primitive construction is not saturated

The source's pattern—compute a rational nullspace, clear denominators separately in each vector, then divide each vector by its gcd—produces valid kernel vectors but not the complete integer lattice.

The exact index of that generated sublattice in the saturated integer kernel is:

- 69×6: (2^{179}=766247770432944429179173513575154591809369561091801088).
- 70×6: (2^{180}=1532495540865888858358347027150309183618739122183602176).

LLL on the source-style primitive basis leaves maximum coefficient 16 (69×6: 35/63 vectors have max coefficient ≤10; 70×6: 36/64). LLL on the saturated lattice reduces the maximum to 2. This is not a cosmetic improvement: the source-style basis omits an enormous index-(2^{179}) or index-(2^{180}) part of the integer relation lattice.

## 4. Formula-row provenance discrepancies found during the full-matrix inspection

The exponent matrices are not yet a clean, frozen transcription of the article's formulas:

- In `constant_matrix_analysis.py`, the active `m_qu_s` vector is ((2,2,-1/2,0,5,-1)), which encodes (K(\ln N)^5/(\sqrt\pi N^{1/3})). Its nearby source comment instead says ((\ln N)^4\pi^{7/2}/N^{1/3}), a different exponent tuple. The active vector agrees with the article formula register; the comment is stale/wrong.
- The article's Ω-c formula in the formula register is (K(\ln N)^6/(\pi^{5/2}N^{1/3})), exponent tuple ((2,2,-5/2,0,6,-1)). Both active matrices instead assign `Omega0_c` the vector ((5,1,2,0,5,-1)), corresponding to (4\pi^2\sqrt K(\ln N)^5/N^{1/3}). This is a different formula, apparently a fitted replacement, while the old formula remains in comments. It must be treated as model-version drift, not a harmless formatting difference.
- The full 70-row file includes K-star; the 69-row file comments it out. The two matrices must remain distinct versions.

## 5. Dimensional audit: the six matrix generators carry no SI dimensions

The basis is ((\sqrt2,\sqrt3,\pi,\ln K,\ln N,N^{1/3})). Every generator is dimensionless. Therefore every monomial encoded by either matrix is dimensionless, irrespective of its exponent tuple.

The matrix cannot by itself encode kg, m, s, J, C, or any other SI dimension. Comparing its monomial values directly with dimensional CODATA quantities requires an explicit dimensionful conversion scale for each independent dimension (or a complete natural-unit convention plus a fixed SI restoration map). No such scale map is part of the exponent matrix itself.

## 6. Verdict

- Exact rank and saturated integer kernels for frozen 69×6 and 70×6 matrices: **COMPUTED**.
- LLL vectors and every direct multiplication: **PASS**, 127 total exact vector checks.
- Source-style row-wise primitive kernels: **NOT SATURATED**, indices (2^{179}) and (2^{180}).
- Formula provenance: **NOT CLEAN**; at least the m_s comment and Ω-c active formula require explicit reconciliation.
- SI dimensional derivation from these six dimensionless generators: **ABSENT**.
- Physical validity of the formulas: **NOT ESTABLISHED** by the integer-kernel result.

Reproducer: `audit_integer_kernel_full_matrices.py`. It also writes JSON files containing all saturated and source-style LLL vectors and their verification metadata.
