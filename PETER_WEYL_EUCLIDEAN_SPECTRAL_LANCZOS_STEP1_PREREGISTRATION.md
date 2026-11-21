# Preregistration: direct Euclidean master-Lanczos continuation from Q1

Status: **frozen before the actual Q1 rank, A1 spectrum or B2 residual is observed.**

## Target

The projected-heat `mu2` workflow constructs the first actual Euclidean master-Lanczos recurrence

\[
M_E Q_0=Q_0A_0+Q_1B_1,
\]

with

\[
Q_0=V_0,
\qquad A_0=V_0^\dagger M_EV_0,
\qquad B_1^\dagger B_1=R_1.
\]

No rank of `Q1` is preregistered.

For every emitted orthonormal `Q1` column, compute the same fixed-cutoff master action

\[
Z_1=M_EQ_1=\sum_{v=0}^4 H_v^E H_v^E Q_1
\]

with the same explicitly Hermitian Euclidean node operator and no additional post-action tolerance prune.

Then

\[
\boxed{A_1=Q_1^\dagger M_EQ_1=Q_1^\dagger Z_1.}
\]

The next residual is

\[
\boxed{
R_2^{state}
=Z_1-Q_0B_1^\dagger-Q_1A_1.
}
\]

Full reorthogonalization against all previous retained Lanczos blocks is mandatory before the rank decision.  Its Gram matrix is

\[
\boxed{G_{R_2}=R_2^\dagger R_2=B_2^\dagger B_2\succeq0.}
\]

No expected eigenvalue, rank or termination outcome is frozen.

## Finite spectral-history decision

If a propagated numerical error certificate proves

\[
\|R_2^{state}\|\le\epsilon_{term},
\]

then the fixed-regulator cyclic space

\[
\operatorname{span}\{Q_0,Q_1\}
\]

is invariant under `M_E`, and the two-block Jacobi matrix

\[
J_1=
\begin{pmatrix}
A_0&B_1^\dagger\\
B_1&A_1
\end{pmatrix}
\]

closes the Euclidean boundary history exactly at that regulator:

\[
V_0^\dagger f(M_E)V_0
=E_0^\dagger f(J_1)E_0
\]

for every finite-spectrum Borel function `f`, including the master heat kernel and zero spectral projector.

If `B2` is nonzero, factor the rank-revealed residual into `Q2 B2` and repeat the same recurrence.  Krylov depth is increased because the operator generates a new certified direction, not because a desired zero has not yet appeared.

## Regulator boundary

A finite closure at `Jmax=5/2` is a theorem about the declared truncated Euclidean reference operator only.  It is **not** by itself continuum/refinement closure.  The same spectral quantities must subsequently be tested across a preregistered regulator sequence.  Positive master eigenvalues remain constraint-violation spectrum, not particle masses or physical frequencies.

## Full BQG boundary

Even an exact Euclidean finite-regulator spectral closure does not emit the physical BQG projector.  The production Lorentzian master and quantum `HH <-> D_target`/explicit `D_target` gate remain mandatory.  The same Lanczos architecture is then reused on the complete positive `M_BQG` without changing the downstream relational-source machinery.
