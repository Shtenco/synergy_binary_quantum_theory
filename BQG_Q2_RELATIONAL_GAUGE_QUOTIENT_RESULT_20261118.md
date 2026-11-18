# Global SU(2)^V gauge quotient and two-sided endpoint QFI (2026-11-18)

## Scope

This gate starts from the full half-link QFI on the existing three-node path
and four-cycle spin-1/2 tensor networks. The coordinate space has 12V real
infinitesimal rotation parameters: three axes on each of four half-links at
every vertex.

At each invariant vertex v, the Gauss generator rotates all four half-links
together. Its three parameter directions form a 12V-by-3V matrix G. The
quotient is represented by an orthonormal basis Q of the complement of this
gauge subspace, and the quotient metric is

```text
F_quot = Q^T F_half Q.
```

The gate checks F_half G = 0 numerically, symmetry and positive
semidefiniteness of the quotient, and the dimension 9V expected after
removing the 3V vertex gauge directions.

## Gauge-section independence

Two coordinate sections are compared: fixing half-link 0 and fixing
half-link 3 at every vertex. Raw eigenvalues of the two principal
submatrices are not invariant under a coordinate change and are not
compared. Instead, each section is mapped into the orthogonal quotient and
its generalized spectrum is computed relative to the induced quotient
Gram matrix. These generalized spectra must agree with the orthogonal
quotient spectrum.

This is the appropriate finite-dimensional check that the quotient result
does not depend on which half-link is selected as a coordinate reference.

## Two-sided internal-edge parameters

For every internal edge (u, leg_u, v, leg_v), the parameter map retains
three independent infinitesimal rotations at (u, leg_u) and three at
(v, leg_v). Thus there are 6E endpoint parameters. The resulting edge
endpoint QFI is the pullback of the full half-link QFI, not a replacement
for the global gauge quotient. Boundary half-links remain in the full
coordinate space because the vertex Gauss transformation acts on all four
legs, including boundary legs.

The two endpoint frames are independent coordinates here. A physical
relative-holonomy convention or edge-frame redundancy would require an
additional derived transformation law; it is not assumed by this gate.

## Interpretation boundary

The quotient can retain additional null directions after the 3V Gauss
directions are removed. Those are reported as quotient nullity, not hidden
by epsilon regularization. They may reflect stabilizers, graph topology,
boundary-state degeneracy, or further redundancies; this finite control does
not decide which interpretation is physical.

This gate does not establish a nondegenerate metric, arbitrary-graph
validity, a dynamical holonomy model, a BQG history-amplitude measure, a
causal-PL pushforward, or a continuum spectral dimension.

## Reproducibility

The executable exports one JSON record for the path and cycle fixtures.
The workflow runs the new tests and exports
`verification_results/Q2_RELATIONAL_GAUGE_QUOTIENT.json`.


## Classification of the residual null modes

A further diagnostic constructs the real-parameter tangent map of the
normalized boundary-state ray. For each of the 12V half-link/axis
parameters, the derivative of the contracted boundary state is projected
orthogonally to the complete complex state vector. This removes both the
normalization direction and an unobservable global phase.

Writing the resulting complex tangent matrix as D, its realification is

```text
R(D) = [ Re(D) ]
       [ Im(D) ] .
```

The pure-state multiparameter QFI must satisfy the exact identity

```text
F_half = 4 R(D)^T R(D).
```

A dedicated regression test checks this matrix identity, not merely equality
of numerical ranks. The tangent map is then restricted to the orthogonal
gauge quotient. Its SVD rank and nullity are compared with the quotient-QFI
rank and nullity.

The additional null vectors are therefore identified, within this finite
model, as infinitesimal stabilizers of the normalized boundary-state ray:
they generate no first-order change of the projective boundary state.
Their half-link participation by vertex and leg is recorded to expose where
the stabilizer directions live. This is a precise local classification, not
a claim that every such direction is a fundamental gauge symmetry of BQG.
The interpretation may depend on the intertwiner state, graph contraction,
and retained boundary.

The new regression checks the tangent/QFI identity and equality of
quotient tangent and QFI nullities. No numerical values are asserted here
until the workflow executes this new gate.
