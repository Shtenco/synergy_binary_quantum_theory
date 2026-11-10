# BQG multi-scale channel-composition gate

## Why the five-step ancestry path is not yet a global refinement fidelity

The five one-step channel-flow matrices are finite, internally checked
representation-RG data. Their row-stochasticity establishes completeness of
the target channel resolution at each individual step. The dynamic-programming
solver then selects the path maximizing the product of the selected one-step
overlaps.

That objective is well-defined, but a product of pairwise overlaps is not
automatically equal to the overlap of the **composed microscopic refinement
map**. The one-step matrices retain only scalar traces of channel-to-channel
projector overlaps. They discard operator-valued information, in particular
off-diagonal blocks between distinct target master channels. Those blocks can
carry coherence into the next refinement step.

Consequently, the path
\[
a_0^{(3/2)}\to a_1^{(2)}\to a_1^{(5/2)}
\to a_1^{(3)}\to a_2^{(7/2)}\to a_2^{(4)}
\]
is currently a **maximum-product local ancestry candidate**, not yet a
demonstrated global multi-scale projector lineage.

## Exact finite objects

Let \(\mathcal H_j\) be the finite carrier at scale \(j\), let
\(I_j:\mathcal H_j\to\mathcal H_{j+1}\) be the exact binary refinement
isometry, and let \(P_{j,a}\) be the isolated source master-channel
projectors. Define the image projector block
\[
Q_{j,a}=I_jP_{j,a}I_j^\dagger.
\]

For target master projectors \(P_{j+1,b}\), the current one-step matrix is
\[
F_j^{ab}=\frac{1}{d_{j,a}}
\operatorname{Tr}(Q_{j,a}P_{j+1,b}),
\qquad d_{j,a}=\operatorname{Tr}P_{j,a}.
\]

If the target projectors resolve the complete image sector and \(I_j\) is
isometric, then
\[
\sum_bF_j^{ab}=1.
\]
This is a finite completeness identity, not a physical Markov law.

For two steps, compose the microscopic isometries before projecting:
\[
I_{j:j+2}=I_{j+1}I_j,
\qquad
Q_{j,a}^{(2)}
=I_{j+1}I_jP_{j,a}I_j^\dagger I_{j+1}^\dagger.
\]
The direct two-step endpoint overlap is
\[
G_{j:j+2}^{ac}
=\frac{1}{d_{j,a}}
\operatorname{Tr}\left(Q_{j,a}^{(2)}P_{j+2,c}\right).
\]

The pairwise scalar matrices predict instead
\[
(F_jF_{j+1})^{ac}
=\sum_bF_j^{ab}F_{j+1}^{bc}.
\]

This product has a precise interpretation. Define the channel replacement
map
\[
\mathcal R_{j+1}(X)
=\sum_b\frac{\operatorname{Tr}(P_{j+1,b}X)}{d_{j+1,b}}
P_{j+1,b}.
\]
Then
\[
(F_jF_{j+1})^{ac}
=\frac{1}{d_{j,a}}\operatorname{Tr}\left[
I_{j+1}\mathcal R_{j+1}(Q_{j,a})I_{j+1}^\dagger
P_{j+2,c}\right].
\]
Thus multiplying the scalar matrices is equivalent to applying an
intermediate measure-and-prepare channel: it measures the intermediate
master-channel label and replaces the state within that channel by the
normalized identity on that sector. It discards both off-diagonal channel
coherences and traceless structure internal to each channel block.

The direct composed overlap \(G_{j:j+2}\), by contrast, propagates
\(Q_{j,a}\) through the exact microscopic isometry without inserting
\(\mathcal R_{j+1}\). Equality is therefore a substantive dynamical/
representation-theoretic condition, not a consequence of row stochasticity.

The first composition defect to measure is therefore
\[
\boxed{
\Delta_{j}^{(2)}
=G_{j:j+2}-F_jF_{j+1}.
}
\]

A zero defect would support a two-step scalar transfer description for the
tested sectors. A nonzero defect is not a failure of either one-step matrix:
it shows that the scalar channel-flow matrices are not sufficient statistics
for composed refinement.

## What information can produce a defect?

For each source channel, the one-step image can be decomposed using the
target master resolution:
\[
Q_{j,a}
=\sum_{b,b'}P_{j+1,b}Q_{j,a}P_{j+1,b'}.
\]
The current matrix \(F_j^{ab}\) records only
\[
\operatorname{Tr}(P_{j+1,b}Q_{j,a}),
\]
which is a diagonal-block trace. It does not retain the off-diagonal operators
\[
P_{j+1,b}Q_{j,a}P_{j+1,b'},\qquad b\ne b',
\]
nor the full matrix structure within a degenerate channel. These terms can
contribute after the next embedding and target projection.

A useful diagnostic therefore includes both:
- diagonal one-step weights \(F_j^{ab}\);
- off-diagonal block norms and, where feasible, the actual block operators.

If symmetry forces those off-diagonal blocks to vanish, or if the next-step
pullback is scalar on every intermediate channel block, a Markov-like
factorization may follow. Neither condition should be assumed from
row-stochasticity alone.

## Required numerical gate

The next calculation should use the **same frozen binary refinement
isometries and master operators** that generated the five one-step matrices.
No post-hoc scalar fitting, channel reordering, or change of operator
normalization is allowed.

For at least the first two consecutive refinements, report:

1. The direct composed isometry \(I_{j+1}I_j\) and its isometry defect.
2. The direct endpoint matrix \(G_{j:j+2}\).
3. The matrix product \(F_jF_{j+1}\).
4. The full defect \(\Delta_j^{(2)}\), in operator and Frobenius norms.
5. The intermediate off-diagonal block norms
   \(\|P_bQ_{j,a}P_{b'}\|\), for \(b\ne b'\).
6. The endpoint projector-resolution and Gauss-projection defects.
7. The result under an independently repeated cutoff or a second adjacent
   pair of spin scales.

Only after this test should the maximum-product path be interpreted as a
coherent multi-scale ancestry rather than a convenient local path score.

## Claim boundary

The five one-step matrices and their maximum-product path remain valid finite
representation-RG results. This gate asks a stronger question about their
composition. It is not a claim that the current path is wrong, nor does a
successful two-step factorization prove continuum convergence, physical
Hamiltonian evolution, or Lorentzian spacetime emergence. The full
graph-changing habitat residuals \(E_a,F_a\), off-carrier leakage, and
independent physical master gaps remain separate gates.
