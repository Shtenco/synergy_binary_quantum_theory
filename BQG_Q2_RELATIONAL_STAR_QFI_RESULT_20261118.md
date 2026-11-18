# Q2 four-valent star QFI and gauge-null result (2026-11-18)

For one normalized four-valent spin-1/2 SU(2) intertwiner, consider the
twelve infinitesimal twist parameters (four incident half-links times three
SU(2) axes). Since the state is a singlet, the pure-state QFI is its generator
covariance. SU(2) invariance gives the exact factorization

```text
F_(i,a),(j,b) = C_ij delta_ab
C_ij = (4/3) <iota | J_i . J_j | iota>
```

The diagonal entries are one because J_i^2=3/4. The Gauss constraint
sum_i J_i |iota>=0 implies C times the all-ones vector is zero. Therefore
the full 12-parameter QFI has three exact null directions: a common SU(2)
rotation of all four incident half-links. It is not a positive-definite
metric on four independent link twists.

After fixing the fourth link's twist to zero, the induced 3-by-3 relative
metric is the leading principal block C[0:3,0:3]. It is positive definite
for generic sampled intertwiners. This gauge-fixed metric retains
state-dependent off-diagonal couplings between incident links. It must not
be replaced by four independent scalar edge conductances.

The gate checks symmetry, unit diagonal, positive semidefiniteness, the
Gauss-law null vectors, and positivity after this reference-link gauge fix.
This is a local invariant-node result; it does not yet include the full
multi-node mixed-state reduction, higher-spin sums, BQG history amplitudes,
or a continuum limit.
