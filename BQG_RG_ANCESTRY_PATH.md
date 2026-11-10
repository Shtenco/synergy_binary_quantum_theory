# BQG ancestry-consistent RG path

The local master spectrum alone does not identify the refinement lineage. The first multiplicity branch already shows that microscopic binary blocking can preferentially enter a non-lowest master eigenchannel.

Given full channel-flow matrices

\[
\mathcal F_j^{ab}
=\frac12\operatorname{Tr}
(P^{\rm block}_{j,a}P^{\rm master}_{j+1/2,b}),
\]

the finite ancestry path is selected without fitting by maximizing

\[
\prod_j\mathcal F_j^{a_j a_{j+1}},
\]

or equivalently the sum of logarithms. The dynamic-programming implementation validates each matrix's dimensions against the source and target multiplicities, checks finiteness, row stochasticity, pass flags, and consecutive spin steps.

## Five-step result (2026-11-10)

The five consecutive matrices s2=3,4,5,6,7 pass their finite row-completeness, paired-spectrum, generated-pair and Gauss-projection checks. The maximum-product ancestry path is

\[
\boxed{
a_0^{(j=3/2)}
\to a_1^{(j=2)}
\to a_1^{(j=5/2)}
\to a_1^{(j=3)}
\to a_2^{(j=7/2)}
\to a_2^{(j=4)}.
}
\]

The step overlaps are

\[
0.9459578674471627,\quad
0.8688492971205302,\quad
0.8826884371011211,\quad
0.7836251466890334,\quad
0.866888447374884,
\]

and the cumulative product is

\[
\boxed{0.49282787562801145}.
\]

The all-low-channel path over the same five steps has product 4.694570240506031e-05, so the maximum-overlap path product is about 1.0498e4 times larger.

The lineage changes channel index at the first split and again at j=3 -> 7/2, where the target multiplicity increases. It then continues on channel 2 for the final step. Thus the correct finite object is the full matrix ancestry, not a rule that follows the lowest eigenvalue or keeps a fixed channel index.

## Next gate: direct multi-scale composition

The maximum-product path is a local one-step ancestry objective. It is not yet proven equal to the overlap of the composed microscopic refinement isometry. The next required comparison is the direct two-step endpoint matrix against the product of the two one-step matrices, including intermediate off-diagonal channel blocks. The exact definitions and acceptance outputs are specified in [BQG_MULTISCALE_CHANNEL_COMPOSITION_GATE.md](BQG_MULTISCALE_CHANNEL_COMPOSITION_GATE.md).

## What this does and does not establish

This is a finite representation-RG lineage derived from the fixed binary refinement tensor. It is not a physical-time trajectory, a probability, or a continuum theorem. The next required tests are the full graph-changing habitat residuals E_a,F_a, the off-carrier leakage, and independently defined master gaps. Only after those are available can the scale-free physical-projector refinement criterion be assessed.
