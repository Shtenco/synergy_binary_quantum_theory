# Holonomy-to-boundary projective tangent pullback (2026-11-18)

## Question

The preceding finite holonomy gate computed the vertex-gauge quotient in
link-holonomy space. The graph-QFI gate computed the differential of a
normalized boundary-state ray under internal edge twists. These are different
spaces and their kernels must not be equated without an explicit map.

## Explicit probe map

For a graph with \(V\) four-valent nodes and \(E\) internal edges, order the
full half-link infinitesimal coordinates as \(12V\) real parameters. For
each oriented edge \(e=(u,l_u,v,l_v)\), the existing graph-QFI convention
inserts the three generators on its first endpoint \((u,l_u)\). Define
\[
P:\mathbb R^{3E}\to\mathbb R^{12V},\qquad
P_{12u+3l_u+a,\,3e+a}=1
\]
and all other entries zero. The actual contraction's projective differential
\(D_\partial\) then induces
\[
D_e=D_\partial P,\qquad
F_e=4\,\mathcal R(D_e)^T\mathcal R(D_e)
=P^T F_{\rm half}P,
\]
where \(\mathcal R(D)=(\Re D,\Im D)^T\).

This is not a new physical postulate: it is the pullback of the already
implemented half-link tangent under the exact edge-probe injection used by
the existing graph-QFI calculation.

## Numerical gate

The new script compares, for randomized normalized four-valent spin-1/2
intertwiners on the path and cycle fixtures:

- the edge tangent obtained by pulling back the full projective half-link
  tangent;
- the direct edge QFI from the existing graph contraction;
- the matrix identity \(F_e=P^T F_{\rm half}P\);
- tangent rank and QFI rank, and PSD up to roundoff.

It also reports the kinematic holonomy gauge-orbit rank from the previous
gate, but explicitly does not subtract that orbit from the boundary tangent.
The latter would require the simultaneous action on internal holonomies and
boundary frames, which has not yet been implemented.

## Interpretation

This closes the differential/QFI pullback for the chosen first-endpoint
probe section on the two finite graph fixtures. It does not establish that
the probe section is independent of endpoint choice, that the full vertex
gauge action is represented on the boundary Hilbert space, or that the
holonomy quotient and projective stabilizer coincide. Those are separate
remaining gates. No BQG history measure, physical coarse conductance, or
continuum limit follows.
