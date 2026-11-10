# BQG full multiplicity-channel flow matrix

The finite representation-RG object is the full channel-flow matrix

\[
\mathcal F_j^{ab}
=\frac12\operatorname{Tr}
\left(P_{j\to j+1/2,a}^{\rm block}
P_{j+1/2,b}^{\rm master}\right),
\]

where a runs over every isolated source master multiplicity channel and b over every target master channel. The target projectors resolve the full target [2,2]-isotypic sector, so each source row should sum to one. This is an RG-compatibility/projector-overlap matrix, **not** a physical Markov transition probability.

## Five consecutive finite flow matrices

All five shards for s2=3,4,5,6,7, i.e. j=3/2 -> 2 -> 5/2 -> 3 -> 7/2 -> 4, completed with passed=true in workflow run [42797583866](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42797583866). The five JSON inputs are persisted under docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_{3,4,5,6,7}.json.

\[
\mathcal F_{3/2\to2}
=\begin{pmatrix}0.0540421326&0.9459578674\end{pmatrix}.
\]

\[
\mathcal F_{2\to5/2}
=\begin{pmatrix}
0.9656128797&0.0343871203\\
0.1311507029&0.8688492971
\end{pmatrix}.
\]

\[
\mathcal F_{5/2\to3}
=\begin{pmatrix}
0.8397696212&0.1602303788\\
0.1173115629&0.8826884371
\end{pmatrix}.
\]

\[
\mathcal F_{3\to7/2}
=\begin{pmatrix}
0.0010798725&0.9126526548&0.0862674727\\
0.0118507795&0.2045240738&0.7836251467
\end{pmatrix}.
\]

\[
\mathcal F_{7/2\to4}
=\begin{pmatrix}
0.9920365891&0.0045659110&0.0033975000\\
0.0019135244&0.8233707897&0.1747156859\\
0.0060241961&0.1270873565&0.8668884474
\end{pmatrix}.
\]

Across the five matrices:
- maximum row-stochasticity defect: 1.11e-15;
- maximum target paired-eigenvalue spread: 2.13e-14;
- maximum generated-pair spread: 2.22e-16;
- Gauss-projection leakage: between 1.67e-31 and 2.89e-30.

These checks support the internal finite-matrix construction and completeness of the channel resolution at these five steps. They do not establish continuum convergence or a physical transition law.

## Maximum-overlap ancestry path

Dynamic programming over the full matrices, with no lowest-eigenvalue assumption, gives

\[
\boxed{
a_0^{(3/2)}
\to a_1^{(2)}
\to a_1^{(5/2)}
\to a_1^{(3)}
\to a_2^{(7/2)}
\to a_2^{(4)}.
}
\]

The five step overlaps are

\[
0.9459578674471627,\quad
0.8688492971205302,\quad
0.8826884371011211,\quad
0.7836251466890334,\quad
0.866888447374884.
\]

Their cumulative product is

\[
\boxed{0.49282787562801145}.
\]

The rejected shortcut that always selects channel index zero has product

\[
4.694570240506031\times10^{-5},
\]

so the maximum-overlap lineage product is about 1.0498e4 times larger on this five-step interval.

The ancestry switches channel index at the first branching 0 -> 1, and again when the target multiplicity grows at j=3 -> 7/2, where 1 -> 2. It then remains on channel 2 for the last step. This is a finite, data-derived lineage, not a claim that channel indices themselves are invariant labels across changing multiplicities.

The next decisive test is not another overlap product: it is to evaluate the full graph-changing forward/return habitat residuals E_a,F_a, the off-carrier leakage, and independently defined master gaps along this candidate lineage. Only then can the scale-free physical-projector refinement criterion be evaluated.

## Claim boundary

This result is a five-step finite representation-RG overlap calculation. It is not a physical-time trajectory, a probability, a continuum fixed point, a proof of a physical projector, or evidence for an emergent Lorentzian spacetime. The finite matrix checks and the still-open physical/continuum gates must remain separate.
