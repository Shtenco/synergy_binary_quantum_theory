# Q2 three-node relational network QFI result (2026-11-18)

## Construction

A path of three four-valent SU(2)-invariant spin-1/2 intertwiners A-B-C is
formed. The A3-B0 and B3-C0 half-links are contracted with normalized
epsilon singlets. Eight uncontracted half-links remain as the boundary.
Two independent relative twists are inserted on the two shared half-links
of the middle node B, giving six local parameters (three SU(2) axes per link).

For the normalized pure boundary state, the script computes the complete
6-by-6 multiparameter QFI matrix, including cross-link terms. It also checks
the Gauss-law pullthrough identity separately for both shared half-links.

## Result and interpretation

In the tested invariant path family, both 3-by-3 diagonal link blocks are
the identity within numerical tolerance. The cross block is isotropic, c I_3. More strongly, its scalar is fixed
exactly by the intermediate node's pair Casimir:

```text
c = (4/3) <iota_B | J_B0 . J_B3 | iota_B>
```

In the two-dimensional four-valent invariant space, the projected pair
operator J_B0 . J_B3 has eigenvalues -3/4 (pair singlet) and +1/4
(pair triplet). Thus every normalized middle intertwiner obeys
-1 <= c <= 1/3. The cross term is generally nonzero and depends on the
middle-node intertwiner. The full matrix has the block form

```text
F = [ I_3   c I_3 ]
    [ c I_3 I_3   ]
```

so its eigenvalues are 1+c and 1-c, each with multiplicity three. The
Casimir interval implies nonnegative eigenvalues, with a possible zero mode
at the pure pair-singlet endpoint c=-1. This is a direct counterexample to replacing the full
network multiparameter QFI by a sum of independent two-node edge QFIs:
the diagonal local controls do not determine the full metric.

The two-link singlet-gluing probability is 1/16 for normalized invariant
four-valent node states. This follows from the maximally mixed one-leg
marginals (I/2) at each end of each singlet contraction.

## Scope boundary

This is a finite three-node spin-1/2 path calculation with pure invariant
intertwiners and two singlet contractions. It is not yet a full BQG
network measure. Boundary gauge reduction, tracing to mixed states,
higher-spin sectors, simultaneous multi-link correlations on general
graphs, the BQG history-amplitude measure, and the causal-PL pushforward
remain open. No coarse conductance is inferred by multiplying a local
QFI by a Galerkin weight.

## Exact finite numerical cross-check

For 32 seeded random normalized intertwiners (seed 20261118), the independent
local reproduction gives a cross scalar range of approximately
[-0.9536082073, 0.2698601106], maximum deviation from the Casimir formula
below 8e-16, maximum cross-block isotropy defect below 3e-16, and minimum
full-QFI eigenvalue approximately 0.0463917927. The two-link gluing
probability is 1/16. The Gauss-law pullthrough residual is below 3e-16
after using the signed identity K_shared = -sum(K_other).

These are independent local numerical checks of the stated finite algebra,
not evidence that GitHub Actions or the repository test suite has passed.
