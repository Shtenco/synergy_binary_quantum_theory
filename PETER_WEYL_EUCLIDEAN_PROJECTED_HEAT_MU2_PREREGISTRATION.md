# Preregistration: actual Euclidean projected-heat second moment

Status: **frozen before any `mu2` value is observed.**

## Physical target

For the orthonormal q=2 boundary injection `V0` and the frozen finite Euclidean
normal master

\[
M_E=\sum_{v=0}^4 H_v^{E\dagger}H_v^E,
\]

the next actual projected-heat moment is

\[
\boxed{\mu_2=V_0^\dagger M_E^2V_0.}
\]

The existing implementation `k5_peter_weyl_safe_hda_column.py` defines every
Euclidean node Hamiltonian by explicit Hermitian symmetrization of the oriented
primitive and its adjoint. Therefore, in this frozen finite reference
convention,

\[
H_v^{E\dagger}=H_v^E
\]

up to the declared numerical implementation tolerances, and hence

\[
\boxed{M_EV_0=\sum_v H_v^E H_v^E V_0.}
\]

Define the 32 master-image columns

\[
Y_i=\sum_{v=0}^4 H_v^E H_v^E|b_i\rangle.
\]

Then

\[
\boxed{\mu_2=Y^\dagger Y.}
\]

No value, rank, eigenvalue, gap or physical-zero conclusion is preregistered.

## Numerical convention

The history-targeted producer does **not** use `compose_on_sparse`, because that
helper applies an additional final `1e-8` pruning after the second hit. Instead,
for every input `b_i` and node `v`, it evaluates

```text
first  = PW.apply_H_cached_state({b_i:1}, v, Jmax2=5)
second = PW.apply_H_cached_state(first,     v, Jmax2=5)
```

and sums the five `second` states without any new tolerance pruning beyond the
thresholds already internal to the frozen reference implementation.

This makes `Y` a direct numerical realization of the frozen Hermitian reference
operator, not a composition of two separately retained `1e-8` column maps.

## Mandatory self-adjointness / provenance check

The same `Y` columns must obey

\[
V_0^\dagger Y=V_0^\dagger M_EV_0=M_{EE}^{unpruned}.
\]

The reusable Euclidean packet stores a retained `M_EE` together with a rigorous
operator-norm perturbation upper bound

\[
\epsilon_{prune}\ge
\|M_{EE}^{unpruned}-M_{EE}^{retained}\|_2.
\]

Therefore production acceptance requires

\[
\boxed{
\|V_0^\dagger Y-M_{EE}^{retained}\|_2
\le \epsilon_{prune}+\epsilon_{num}.
}
\]

The `mu2` gate MUST fail closed if the Euclidean source packet does not contain
that pruning certificate. An older packet without the bound is insufficient.

## Projected-heat consequence

When accepted, the actual Euclidean short-heat expansion advances from first to
second order:

\[
\boxed{
V_0^\dagger e^{-\tau M_E}V_0
=I-\tau\mu_1+\frac{\tau^2}{2}\mu_2+O(\tau^3),
}
\]

where `mu1=M_EE` is already measured independently from the 160 one-hit columns.

This still does **not** justify

\[
e^{-\tau\mu_1}
\]

as the projected heat kernel. Equality would require invariance of the boundary
subspace under the full master, which is not assumed.

## Claim boundary

This is an Euclidean-sector finite-regulator history moment. It is not the full
BQG physical projector, because the Lorentzian master blocks and actual quantum
`HH <-> Dtarget` habitat closure remain open. It is also not physical time,
`omega`, dark matter, dark energy, or a source-dressed `W_BQG`.
