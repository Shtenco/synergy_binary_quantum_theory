# Two-block relational edge QFI gate — result (2026-11-18)

## Construction

Take two four-valent spin-network nodes A and B. Each node carries four
spin-1/2 half-links and an arbitrary normalized vector in the two-dimensional
SU(2)-invariant intertwiner space:

`|iota_X> in Inv_SU2[(C^2) tensor4]`, for X=A,B.

The shared face is half-link A3 paired with half-link B0. Contract this pair
with the normalized invariant tensor

`|epsilon> = (|01> - |10>) / sqrt(2)`.

The six uncontracted half-links are retained. Define the unnormalized
six-qubit boundary vector and its normalized version by

```
|Psi_tilde_AB> = (<epsilon|_A3,B0 tensor I) |iota_A> |iota_B>
|Psi_AB>       = |Psi_tilde_AB> / ||Psi_tilde_AB||
rho_e          = |Psi_AB><Psi_AB|
```

A relative holonomy is inserted on B0 before the singlet contraction:
`U_B(theta) = exp(-i theta_a J^a_B0)`, with `J^a = sigma_a/2`.

This is not the total B-node Gauss generator. The Gauss constraint is

```
(J^a_B0 + J^a_B1 + J^a_B2 + J^a_B3) |iota_B> = 0.
```

Therefore the half-link generator pulls through the gluing contraction to
the Hermitian generator acting on B's three retained boundary half-links:

```
K_e^a = -(J^a_B1 + J^a_B2 + J^a_B3)
rho_e = |Psi_AB><Psi_AB|
```

Thus the requested relational state and generator are explicit on the
six-qubit boundary Hilbert space. The sign follows from the node Gauss law;
the overall sign does not affect the QFI matrix.

## Exact result

Because each four-valent intertwiner is an SU(2) singlet, the reduced state
of any one half-link is I_2/2. Its Schmidt coefficients across the split
"shared half-link | other three half-links" are both 1/sqrt(2).

Consequently:

1. The probability of the shared-face singlet contraction is exactly
   `p_e = 1/4`, independent of both intertwiner states.
2. The normalized six-leg boundary state is maximally entangled between the
   two two-dimensional spin-1/2 multiplicity sectors.
3. For the relative generators, the covariance matrix is
   `Cov(K_e^a,K_e^b) = delta_ab/4`.
4. The pure-state QFI and normalized QFI metric are

```
F_e[a,b] = 4 Cov(K_e^a,K_e^b) = delta_ab
g_e[a,b] = F_e[a,b] / 4 = delta_ab / 4
```

This is independent of the two intertwiner states in the fixed spin-1/2
channel; it is not a fitted parameter.

## Numerical gate

Reproduction:

```bash
python scripts/q2_relational_edge_qfi_gate.py \
  --seed 20261118 --trials 16 \
  --output verification_results/Q2_RELATIONAL_EDGE_QFI.json
```

The gate uses the canonical two-dimensional four-qubit intertwiner basis
from `spatial_qubit_geometry_gate.py`. For random complex superpositions in
both intertwiner spaces it checks:

- singlet-gluing probability p_e = 1/4;
- both node Gauss constraints;
- the pull-through identity from the shared-half-link twist to K_e;
- the full 3-by-3 QFI matrix F_e = I_3;
- purity and normalization of rho_e.

An independent local numerical cross-check over 16 random complex
intertwiner pairs gave p_e = 0.25, F_e = I_3, and a half-link-to-boundary
generator residual below 3e-16. This is not a report of a completed GitHub
Actions run.

## Interpretation and strict limits

This is a derived nonzero relational QFI for a specified two-node
spin-network gluing channel. It is distinct from:

- the unprojected product-state probe value F_Q/4 = 1/3;
- the vanishing QFI of a global SU(2) gauge rotation after Gauss projection;
- the combinatorial Galerkin conductances obtained by counting PL-dual
  adjacencies.

The normalized microscopic link metric I_3/4 is a candidate local
spin-1/2 edge datum. It is **not yet** the BQG coarse graph conductance.
Open tasks remain: independent boundary-gauge reduction and observable
selection; tracing or coarse-graining external legs; network correlations
when several links are glued simultaneously; superposition over spins and
intertwiners; derivation of the holonomy distribution from the BQG history
amplitude; and pushing the link QFI through the non-equitable causal-to-PL
fibre map. No spectral dimension or continuum dynamics follows from this
local result alone.
