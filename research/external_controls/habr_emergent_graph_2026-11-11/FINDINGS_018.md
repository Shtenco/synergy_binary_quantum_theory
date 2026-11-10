# Findings 018 — continuum-limit and Lorentz-signature claims need additional structure

Scope: deepen the independent audit of the article's TSCO/current/metric claims in FINDINGS_005 and ledger entries J-02, J-06, J-07, T-08, S-01, S-03, S-04. This is a conditional mathematical analysis of the displayed construction, not a claim that no extension could work.

## 1. A positive Gram construction cannot produce Lorentzian signature

The proposed metric is
[
g_{\mu\nu}(x)\sim\left\langle \partial_\mu J^\alpha(x)\,
\partial_\nu J_\alpha(x)\right\rangle.
]
Suppose the contraction over α uses a positive-definite internal inner product and the expectation is a positive state/ensemble average. For any real covector (v^μ),
[
v^\mu g_{\mu\nu}v^\nu
\sim\left\langle\left\|v^\mu\partial_\mu J\right\|^2\right\rangle\ge0.
]
Thus (g) is positive semidefinite. If the gradients span all four directions, it is positive definite; if not, it is degenerate. It cannot have signature ((-+++)).

Therefore a Lorentzian metric requires at least one additional ingredient not supplied by the displayed Gram formula:
- an indefinite internal bilinear form;
- a causal/time orientation structure and a theorem connecting it to the quadratic form;
- a non-positive state pairing with carefully controlled stability;
- or an explicit analytic continuation.

If the index in (J_\alpha) is lowered with a Minkowski metric, then Lorentzian signature has already been inserted before the metric is “derived”. If it is lowered with a Euclidean metric, the no-go above applies. This is a sharp fork that the derivation must resolve.

## 2. SL(2,C) link variables do not supply a positive invariant spinor norm

For the fundamental two-component representation, a generic (g\in SL(2,\mathbb C)) is not unitary. Under Ψ→gΨ, the Hermitian form Ψ†Ψ changes to Ψ†g†gΨ. It is invariant for all Ψ only if (g^†g=I), i.e. for the SU(2) subgroup, not generic SL(2,C).

The current
[
j_{uv}=\operatorname{Im}(\Psi_u^\dagger S_{uv}\Psi_v)
]
therefore is not automatically invariant under local SL(2,C) frame changes. One must state:
1. the transformation law for Ψ at each vertex;
2. the link transformation law;
3. the invariant pairing used to contract spinor indices;
4. the reverse-edge relation and orientation convention;
5. the resulting reality and conservation properties.

The invariant antisymmetric spinor tensor ε may be relevant for SL(2,C) spinors, but it is not interchangeable with the positive Hermitian inner product Ψ†Ψ. The current formula must be re-derived in one consistent representation.

Likewise, (S_{uv}\in SL(2,\mathbb C)) alone does not imply that the graph Hamiltonian is Hermitian. A separate Hilbert-space inner product and edge-adjointness condition are needed for (U=\exp(-iH\Delta t/\hbar)) to be unitary.

## 3. The inverse-square shortcut ensemble is not automatically a local continuum discretization

For the specified three-dimensional cubic background, the number of lattice sites in a shell of radius (r) is asymptotically proportional to (r^2). With
[
q(r)=\frac{1}{6Lr^2},
]
the expected shortcut count per vertex per unit-width shell is therefore of order (1/L). Summing over (O(L)) shells gives an (O(1)) expected shortcut degree, consistent with the observed mean degree near 7.2 after adding six local edges.

But this same calculation shows why the continuum limit is not automatic: a non-vanishing fraction of shortcut opportunities lie at distances proportional to (L). After rescaling the torus to fixed unit diameter, these edges connect points separated by order-one coordinate distances. They do not become local nearest-neighbour edges as (L\to\infty).

Consequently, a normalized graph operator may converge to a nonlocal integral operator, or require a scaling that suppresses the shortcuts, rather than to a local Laplace-Beltrami operator. To claim a local continuum geometry, the theory must establish at least:
- the precise rescaling of graph Laplacian, time, fields, and edge weights;
- convergence of the discrete quadratic forms (e.g. a stated Mosco/form convergence);
- tightness/compactness of interpolated fields;
- control of the long-edge contribution, showing it vanishes or deriving the intended nonlocal continuum operator;
- convergence of distances and volume measures, not only one eigenvalue.

The shortcut kernel is defined using a Euclidean distance (r(u,v)) on the cubic lattice. Thus the current ensemble uses an embedded three-dimensional distance to decide which edges exist. That is a model input, not an emergent metric result. A genuinely relational construction would need to define the edge law without first using the continuum geometry it aims to recover.

## 4. Spectral gap alone does not establish dimension or metric convergence

The scalar λ₂ controls a global relaxation scale for the chosen graph Laplacian. It does not determine:
- local dimension;
- convergence of heat kernels at all scales;
- a four-dimensional spacetime structure;
- causal cones;
- metric signature;
- curvature or Einstein dynamics.

At minimum, a continuum claim needs joint control of heat-kernel asymptotics, volume growth, distance convergence, and a well-defined embedding/interpolation. A gap plateau in the finite random graph ensemble is not a Lorentzian metric certificate.

## 5. Required TSCO-to-continuum proof chain

A valid derivation should expose each implication, without skipping levels:
1. Define the finite-(N) state space, inner product, graph orientation, link variables, and TSCO operator domain.
2. Prove TSCO is well-defined and maps the stated state space into itself.
3. Establish existence/stability of solutions; for a memory term, either enlarge the state to include history or construct a rigorously equivalent Markovian auxiliary system.
4. Derive the discrete continuity equation from the update and its adjointness assumptions; do not replace it by pointwise zero divergence without a stationarity condition.
5. Specify a scaling family and prove a continuum convergence theorem for the dynamics and observables.
6. Derive a nondegenerate metric from the limiting observables, with a signature proof independent of any metric used in the contraction.
7. Derive causal cones and show compatibility with the limiting propagation domain.
8. Only then compare the limiting equations with relativistic field equations.

## 6. Verdict

- Positive-definite current-gradient Gram tensor has Lorentzian signature: **IMPOSSIBLE** (it is positive semidefinite).
- Lorentzian signature from the displayed metric formula: **NOT DERIVED** unless an indefinite structure is supplied; supplying it may make the signature an input.
- Generic SL(2,C) links preserve the displayed Hermitian current pairing: **FALSE** without extra structure.
- SL(2,C) links imply a Hermitian Hamiltonian/unitary TSCO: **NOT IMPLIED**.
- Inverse-square shortcut ensemble automatically converges to a local 3D/4D manifold: **NOT ESTABLISHED**; macroscopic nonlocal edges persist under the natural fixed-diameter rescaling.
- Spectral-gap behavior alone establishes Lorentzian geometry: **FALSE**.
- Full TSCO-to-Lorentzian-continuum derivation: **OPEN / UNPROVED**.
