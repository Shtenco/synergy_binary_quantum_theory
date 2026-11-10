# BQG ancestry-consistent RG path

The local master spectrum alone does not identify the physical refinement
lineage.  The first multiplicity branch already shows that microscopic binary
blocking can preferentially enter a non-lowest master eigenchannel.

Given full channel-flow matrices

\[
\mathcal F_j^{ab}
=
\frac12\operatorname{Tr}
(P^{\rm block}_{j,a}P^{\rm master}_{j+1/2,b}),
\]

the finite ancestry path is selected *without fitting* by maximizing

\[
\boxed{
\prod_j\mathcal F_j^{a_j a_{j+1}}
}
\]

from the unique pre-branch source channel.

Equivalently maximize

\[
\boxed{
\sum_j\log\mathcal F_j^{a_j a_{j+1}}.
}
\]

A dynamic-programming gate computes this path exactly once all finite flow
matrices are available.

This criterion is not a physical-time probability principle.  It is a
representation-RG ancestry criterion: among discrete master eigenchannels, it
finds the lineage most compatible with the already-fixed microscopic binary
refinement tensor.

The resulting path can then be used as the candidate finite lineage on which
the graph-changing habitat residuals \(E_a,F_a\) and ultimately
\(\epsilon_j^{\rm sf}\) are evaluated.


## First four-shard ancestry prefix (2026-11-10)

The first four consecutive full channel-flow matrices, (s_2=3,4,5,6), now pass their row-completeness, paired-spectrum and Gauss-projection checks. Applying the maximum-product ancestry rule gives the finite prefix

```math
a_0^{(j=3/2)}
\to a_1^{(j=2)}
\to a_1^{(j=5/2)}
\to a_1^{(j=3)}
\to a_2^{(j=7/2)}.
```

The step overlaps are

```math
0.9459578674471627,\quad
0.8688492971205302,\quad
0.8826884371011211,\quad
0.7836251466890334,
```

with cumulative product

```math
0.5685020686576171.
```

For comparison, the all-low-channel path has cumulative product (4.732255132750176\times10^{-5}) over the same four steps. The maximum-overlap ancestry prefix is therefore about (1.20\times10^4) times larger.

The path is a finite representation-RG lineage, not a physical-time trajectory or probability. The fifth shard (s_2=7) is still required for the complete five-step path; after that, the lineage must still pass the graph-changing (E_a,F_a) habitat-refinement and scale-free master-residual tests.
