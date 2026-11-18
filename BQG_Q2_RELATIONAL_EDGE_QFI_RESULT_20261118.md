# Two-block relational edge QFI gate — result (2026-11-18)

## Construction

Take two four-valent spin-network nodes (A,B), each with four spin-(1/2)
half-links. At each node choose an arbitrary normalized state in the
two-dimensional invariant intertwiner space
[
|iota_Xangleinmathrm{Inv}_{SU(2)}
[(\mathbb C^2)^{\otimes4}],qquad X=A,B.
]
The shared face is the fourth half-link of (A) and the first half-link of
(B). Glue it with the normalized invariant tensor
[
|\epsilon\rangle=(|01\rangle-|10\rangle)/\sqrt2.
]
The six remaining half-links are retained as the boundary system:
[
|\widetilde\Psi_{AB}\rangle
=(\langle\epsilon|_{a b}\otimes I)
|\iota_A\rangle|\iota_B\rangle,
\qquad
|\Psi_{AB}\rangle=
|\widetilde\Psi_{AB}\rangle/
\|\widetilde\Psi_{AB}\|.
]

A relative holonomy is inserted on the (B)-side shared half-link before
contraction:
[
U_B(\boldsymbol\theta)
=\exp[-i\theta_a J^a_{B,b}],
\qquad J^a=\sigma_a/2.
]
This is a half-link generator, not the total Gauss generator
J^a_{B,tot} = sum_{f=1}^4 J^a_{B,f}, which annihilates every chosen
intertwiner. The node Gauss law gives
J^a_{B,b}|iota_B> = -sum_{r=1}^3 J^a_{B,r}|iota_B>.
Consequently the inserted half-link twist pulls through the gluing map to
the Hermitian relational generator on the retained six-qubit boundary:
K_e^a = -sum_{r=1}^3 J^a_{B,r}.
The corresponding relational density operator is
rho_e = |Psi_AB><Psi_AB|, where |Psi_AB> is the normalized contracted
six-qubit boundary state defined above.

## Exact result

Because each four-valent invariant intertwiner is an SU(2) singlet, its
one-half-link reduced density matrix is (I_2/2). Equivalently, across the
split “shared half-link | other three half-links” it has two equal Schmidt
coefficients (1/\sqrt2). Therefore:

1. the singlet-gluing probability is exactly (p_e=1/4), independent of
   either node's intertwiner;
2. after normalization, the six-leg boundary state is maximally entangled
   between the two two-dimensional spin-(1/2) multiplicity sectors;
3. the relative half-link generators have covariance
   [
   \frac12\langle\{J^a,J^b\}\rangle
   -\langle J^a\rangle\langle J^b\rangle
   =\frac14\delta_{ab};
   ]
4. the pure-state QFI is
   [
   \boxed{F^{(e)}_{ab}=4\,\mathrm{Cov}(J^a,J^b)
   =\delta_{ab}},
   \qquad
   \boxed{g^{(e)}_{ab}=F^{(e)}_{ab}/4
   =\frac14\delta_{ab}}.
   ]

The value is independent of the two intertwiner states, not a fitted
parameter. The script evaluates the full six-qubit contracted boundary
state and its three tangent vectors for random complex superpositions in
both intertwiner spaces. It checks p_e=1/4, the identity QFI matrix, both
node Gauss constraints, and equality between the inserted half-link
derivative and the boundary action of K_e.

## Numerical gate

Reproduction:
```bash
python scripts/q2_relational_edge_qfi_gate.py \
  --seed 20261118 --trials 16 \
  --output verification_results/Q2_RELATIONAL_EDGE_QFI.json
```

The gate uses the repository's canonical two-dimensional four-qubit
intertwiner basis from `spatial_qubit_geometry_gate.py`. A local independent
numerical cross-check over 16 random complex intertwiner pairs gave
(p_e=0.25) and (F_e=I_3) to floating-point precision. This is a local
calculation check, not a report of a completed GitHub Actions run.

## Interpretation and strict limits

This is the first derived nonzero relational QFI in the current q=2
spin-network carrier: it comes from an explicit two-node invariant gluing
channel and a relative half-link twist. It is distinct from:

- the unprojected product-state probe value (F_Q/4=1/3);
- the vanishing QFI of a global SU(2) gauge rotation after Gauss projection;
- the combinatorial Galerkin conductances obtained by counting PL-dual
  adjacencies.

The normalized microscopic link metric (I_3/4) is a candidate local
spin-(1/2) edge datum. It is **not yet** the BQG coarse graph conductance.
The following remain open: independent boundary-gauge reduction and
observable selection; tracing or coarse-graining external legs; superposition
over spins and intertwiners; derivation of the holonomy distribution from the
BQG history amplitude; and pushing the link QFI through the non-equitable
causal-to-PL fibre map. No spectral dimension or continuum dynamics follows
from this local result alone.
