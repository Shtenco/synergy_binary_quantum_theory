# Joint vertex-Gauss and singlet-edge null quotient (2026-11-18)

## Two exact null families

The full parameter space has \(12V\) half-link/axis coordinates. Two
different null families are present in the explicit finite contraction.

1. **Vertex Gauss generators.** At node \(v\), the same SU(2) axis rotation
   acts on all four legs. These \(3V\) generators annihilate the full QFI.
2. **Singlet-edge common rotations.** For every internal edge \(e=(u,\ell_u;v,\ell_v)\),
   a generator rotates the two contracted half-links on both endpoints by
   the same axis. Invariance of the normalized two-leg tensor \(C_{1/2}\)
   gives \(D_{u,\ell_u}+D_{v,\ell_v}=0\). These contribute \(3E\) further
   candidate null generators.

The two families are not identified with each other. Their combined
generator matrix is \(H=[G_{\rm vertex},C_{\rm edge}]\), and the joint
quotient basis is the orthogonal complement of \(\operatorname{im}H\).

## Finite graph results

For the path fixture \(V=3,E=2\):
\[
\operatorname{rank}H=3V+3E=15,\qquad
\dim\mathcal H_{\rm joint}=36-15=21.
\]
For the cycle fixture \(V=4,E=4\):
\[
\operatorname{rank}H=3V+3E=24,\qquad
\dim\mathcal H_{\rm joint}=48-24=24.
\]
The gate computes the full QFI and projective tangent, verifies both null
families, checks the rank of their joint span, and requires the joint
quotient QFI to be positive definite with matching tangent/QFI rank. Thus,
for these two finite fixtures, the earlier quotient nullities (6 and 12)
are exactly accounted for by the \(3E\) singlet-edge common rotations:
there are no additional projective stabilizers after this joint quotient.

## Interpretation and limits

This is a structural classification of null directions in the current
normalized spin-one-half singlet contraction. It does not establish that
the edge-common family is a new fundamental gauge symmetry. It may be an
internal tensor-frame redundancy of the chosen gluing representation.
A general graph with different boundary conditions may have linear
dependencies between vertex and edge generators; the implementation
computes the combined rank rather than assuming \(3V+3E\).

Files:
- \`scripts/q2_relational_joint_gauge_edge_quotient_gate.py\`
- \`tests/test_q2_relational_joint_gauge_edge_quotient_gate.py\`
- \`BQG_Q2_RELATIONAL_JOINT_GAUGE_EDGE_QUOTIENT_RESULT_20261118.md\`

CI exports \`verification_results/Q2_RELATIONAL_JOINT_GAUGE_EDGE_QUOTIENT.json\`.


## Topological rank theorem (extended audit)

The combined generator rank can be predicted exactly from graph topology.
For one Lie-algebra axis, let a_v be the coefficient of the vertex column
and b_e the coefficient of the edge column in a linear dependence.
A boundary half-link forces a_v=0. An internal edge e=(u,v) gives
a_u+b_e=0 and a_v+b_e=0, hence a_u=a_v. Therefore a connected component
with a boundary half-link has no nonzero dependence; a component with no
boundary half-links has exactly one dependence (all a_v equal and b_e=-a_v).

If c_closed is the number of connected components with no external
half-links, the exact rank is

    rank(H) = 3 * (V + E - c_closed)

and the quotient dimension is

    12V - 3 * (V + E - c_closed).

This is a theorem about the generator map, not a claim that the QFI quotient
is positive definite. A fully closed graph has no open boundary legs in the
present contraction setup, so this linear-algebra result alone does not
provide a nontrivial boundary pure-state QFI.

Tests compare the SVD rank with the formula on the open path and cycle,
a closed complete five-node graph (10 internal edges), and two disconnected
closed components plus an isolated open node.
