# Findings 012 — the 30×6 LLL calculation reduces a non-saturated sublattice

Source: draft/stress_kernel.py, blob SHA 0fa851abb181a121da11ae309fdc6f0de2f88f21.

## 1. Exact matrix and rational nullspace

The source matrix is exactly 30×6 with rational entries, rank 6, so the rational left-nullspace has dimension 24. This part reproduces exactly.

The source then obtains a rational nullspace basis, clears denominators separately in each basis vector, divides each vector by its coefficient gcd, and applies LLL to the resulting 24×30 integer matrix.

Every resulting row is indeed an integer kernel vector. The problem is completeness: these 24 individually primitive vectors do not span the full integer kernel lattice.

## 2. Smith-normal-form test of lattice saturation

Let A=6Mᵀ, an integer 6×30 matrix with the same integer kernel as Mᵀ. An exact Smith normal decomposition U A V = D, with U and V unimodular, gives the saturated integer kernel directly: the last 24 columns of V form a Z-basis of ker_Z(A).

Independent calculation gives:

- SNF nonzero diagonal: (1, 6, 6, 6, 6, 6);
- saturated integer kernel rank: 24;
- source's row-wise-cleared basis: 24×30 and every row is an exact kernel vector;
- index of source's generated lattice inside the full integer kernel: 635603431212363775732957584 ≈ 6.356×10²⁶.

The index is the absolute determinant of the integer change-of-basis matrix between the source's 24 rows and the saturated kernel basis. It is also the product of the nonzero Smith invariants of the source's 24×30 basis matrix.

Therefore the source procedure does not produce the full integer kernel. It produces a proper sublattice of enormous index. Validity of each printed identity is not in question; completeness and the search for short identities are.

## 3. LLL on the source sublattice versus LLL on the saturated kernel

Both calculations use SymPy LLL with δ=3/4:

| Measure | Source row-wise-cleared lattice | Saturated integer kernel |
|---|---:|---:|
| Basis vectors | 24 | 24 |
| Maximum absolute coefficient after LLL | 66 | 2 |
| Basis rows with max absolute coefficient ≤10 | 2 | 24 |

The saturated basis has 24 exact, linearly independent kernel vectors, all with maximum coefficient at most 2. The source's sublattice LLL finds only 2 such vectors under the same threshold.

This is a decisive algorithmic issue: LLL only reduces the lattice it is given. It cannot discover short integer relations excluded by a non-saturated starting lattice.

## 4. What this does and does not invalidate

- The source's printed LLL vectors remain valid exact exponent identities, provided each is checked against Mᵀv=0.
- The source's claim to have searched the full integer kernel is not justified by its construction.
- Counts of short invariants obtained from that sublattice are incomplete and cannot be treated as counts of short relations in the full integer kernel.
- A saturated kernel basis is not unique, and LLL output depends on basis conventions and δ. The exact lattice is canonical; a particular reduced basis is not.
- Even a complete integer-kernel basis establishes algebraic dependencies among the chosen formula rows, not independent empirical confirmations of the physical theory.

## 5. Reproducibility artifact

Run audit_integer_kernel.py with SymPy. It reconstructs the exact 30×6 matrix, verifies the source procedure, obtains a saturated kernel from Smith normal decomposition, computes the lattice index, performs both LLL reductions, and checks all returned basis vectors exactly.

## 6. Verdict

- Rational rank 6 and nullity 24: PASS.
- Source's 24 row-wise-cleared vectors are valid kernel vectors: PASS.
- Those vectors form the full integer kernel: FAIL; index ≈6.356×10²⁶.
- LLL search over the full integer kernel in draft/stress_kernel.py: FAIL AS IMPLEMENTED.
- Saturated integer kernel and its LLL reduction: independently reproduced; 24 basis vectors with max coefficient ≤2.
- Article's separate 55×6 / 29-vector claim: still not reproduced; this finding is specifically about the checked-in 30×6 draft matrix.