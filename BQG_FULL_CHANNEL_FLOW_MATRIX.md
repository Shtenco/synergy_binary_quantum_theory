# BQG full multiplicity-channel flow matrix

The previous low-to-low overlap test temporarily assumed that the physically
relevant local master channel is always the lowest eigenchannel.  The first
multiplicity branching \(j=3/2\to2\) falsifies that shortcut: binary
microscopic refinement overlaps only about \(5.4\%\) with the lowest \(j=2\)
master channel and therefore about \(94.6\%\) with the other channel.

The correct finite representation-RG object is therefore the **full channel
flow matrix**

\[
\boxed{
\mathcal F_j^{ab}
=
\frac12
\operatorname{Tr}
\left(
P_{j\to j+1/2,a}^{\rm block}
P_{j+1/2,b}^{\rm master}
\right).
}
\]

Here:

- \(a\) labels every isolated source master multiplicity channel;
- \(b\) labels every target master multiplicity channel;
- \(P^{\rm block}_{a}\) is the target \([2,2]\) copy generated from source
  channel \(a\) by the exact binary-ancilla symmetric blocking tensor.

For each fixed source channel \(a\),

\[
\boxed{
\sum_b \mathcal F_j^{ab}=1
}
\]

within numerical precision, because the target master channels form an
orthogonal resolution of the full target \([2,2]\)-isotypic sector.

This matrix is **not** a physical Markov transition probability.  It is an
RG-compatibility/fidelity matrix.

Its purpose is to identify ancestry-consistent channel trajectories without
the unsupported rule "always choose the lowest local master eigenvalue".

The next canonical finite question is whether the sequence of flow matrices
contains a high-overlap lineage across increasing \(j\), and whether that
lineage stabilizes.


## First completed full-flow shard (2026-11-10)

The (s_2=3), (j=3/2\to2) matrix has now been computed by the full all-channel gate and passed its internal checks:

```math
\mathcal F_{3/2\to2}
=
\begin{pmatrix}
0.05404213255283678 & 0.9459578674471627
\end{pmatrix}.
```

The row sum is (0.9999999999999994), the row-stochasticity defect is (5.55\times10^{-16}), target paired-eigenvalue spreads are below (5.4\times10^{-15}), and the Gauss-projection leakage is (1.67\times10^{-31}).

The numerical result is persisted at `docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_3.json`. It confirms that the first source channel flows predominantly to target channel 1 (the higher local master eigenvalue), not target channel 0. The remaining four consecutive scale shards are required before the global ancestry path can be selected.


The (s_2=4), (j=2\to5/2) shard has also passed. Its complete matrix is

```math
\mathcal F_{2\to5/2}
=
\begin{pmatrix}
0.965612879727417 & 0.03438712027258314\\
0.13115070287946956 & 0.8688492971205302
\end{pmatrix}.
```

Both row sums differ from one by at most (2.22\times10^{-16}); target paired-eigenvalue spreads are below (4.5\times10^{-15}); Gauss-projection leakage is (2.30\times10^{-31}).

Combining the first two matrices, the maximum-overlap ancestry path is so far

```math
j=3/2:\ a_0
\longrightarrow
j=2:\ a_1
\longrightarrow
j=5/2:\ a_1,
```

with step overlaps (0.9459578674) and (0.8688492971). Thus the dynamically inherited branch switches to the upper channel at the first split and remains on the upper channel at the next scale. The low-to-low shortcut would instead multiply (0.0540421\times0.9656129\), which is not the microscopic ancestry path.

The (s_2=4) result is persisted in `docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_4.json`. The (s_2=5,6,7) shards remain in progress.
