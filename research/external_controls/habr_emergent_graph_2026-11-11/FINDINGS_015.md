# Findings 015 — the 30×6 stress-kernel research script mislabels minor gcds and basis geometry

Source: draft/stress_kernel_research.py, blob SHA 05f599a8cb76afc32e83210c4e848a7fa13e357e.

## 1. The so-called “selective SNF” is not a Smith normal form computation

The script forms M_int=6M, enumerates 6×6 row minors, and stops after the first 20 nonzero determinants. It then calls their gcd `gcd_val` and uses it to discuss the torsion order.

For this exact old 30×6 matrix:

- gcd of the first 20 nonzero maximal minors: 93,312;
- exact Smith diagonal of 6Mᵀ: (1,6,6,6,6,6);
- exact gcd of all maximal minors (product of nonzero Smith invariants): 7,776.

The sampled-minor gcd is 12 times the exact gcd. It is only an upper multiple of the gcd over all maximal minors, not the exact value. A sample whose gcd happens to be 1 would certify that the full gcd is 1, but a sample gcd greater than 1 cannot establish the exact invariant. Here it overstates the value.

More fundamentally, the script's interpretation of this determinant gcd as torsion in the group of relations is incorrect. The integer left kernel is a subgroup of Z^30 and is necessarily free abelian. The Smith invariants of the map 6Mᵀ describe the cokernel of that chosen integer matrix, not torsion in the left-kernel relation group.

## 2. The Gram-spectrum and angle analysis is basis-dependent

The script constructs one row-wise denominator-cleared rational-nullspace basis, then studies its Gram matrix, eigenvalues, condition number, and pairwise angles. These are properties of the chosen basis, not invariants of the kernel subspace or integer lattice.

An exact unimodular basis change preserves the same integer lattice but changes the Gram matrix by congruence. Recomputing the same old matrix using the saturated basis and then its LLL-reduced basis demonstrates the instability:

| Basis of the same saturated integer kernel | Gram condition number | Nearly orthogonal pairs (|cos θ|<0.01) |
|---|---:|---:|
| Smith-derived saturated basis | 1.073×10^11 | 0 |
| LLL-reduced saturated basis | 206.4 | 85 |

The lattice and its rational span are identical in these two rows. Therefore the script's reported orthogonal/parallel pair counts, Gram eigenvalue distribution, and condition number cannot be interpreted as intrinsic physical clustering or independent directions.

## 3. The script uses the same non-saturated nullspace construction

get_kernel_vectors() clears denominators separately and divides each row by its gcd. As shown in FINDINGS_012, for this old 30×6 matrix the resulting 24 vectors span a sublattice of index ≈6.356×10^26, not the full integer kernel.

Thus the Gram analysis is performed on an incomplete sublattice, and then its basis-dependent geometry is interpreted physically. Both steps need replacement.

## 4. The input matrix itself is an obsolete/inconsistent exponent version

FINDINGS_013 shows that 24 of its 30 exponent rows differ from the newer first-30-row matrix after matching the N versus N^(1/3) column convention. For example, the old hbar and c rows omit the prime-factor contributions from K=6. Consequently the script's rank/kernel/Gram results do not describe the newer formula table.

## 5. Correct mathematical workflow

1. Freeze a single formula version and verify every exponent row from the displayed formula independently.
2. Clear rational denominators globally and compute an exact Smith normal decomposition.
3. Obtain the saturated integer kernel from the unimodular transformation, not by independently primitive-scaling rational nullspace vectors.
4. Apply LLL to that saturated lattice and verify every vector exactly.
5. If studying geometry of the lattice, use basis-invariant quantities (e.g. determinant/covolume, successive minima bounds, theta series where feasible), and state the metric being used. Do not interpret Gram eigenvalues of an arbitrary basis as physical clustering.
6. Keep cokernel torsion and the free integer relation group conceptually separate.

## 6. Verdict

- First-20-minor gcd as exact SNF invariant: FAIL; 93,312 vs exact 7,776.
- Torsion interpretation as a property of the integer relation kernel: MATHEMATICALLY INCORRECT.
- Gram eigenvalue/angle clustering as intrinsic physical structure: NOT INVARIANT UNDER BASIS CHANGE.
- Integer-kernel completeness: FAIL; index ≈6.356×10^26.
- Input exponent matrix: obsolete relative to newer first-30-row matrix.
- Physical interpretation of this script's Gram/clustering outputs: NOT ESTABLISHED.