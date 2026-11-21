# Preregistration: reusable Euclidean K2-even constraint-generated carrier

Status: **frozen before the 800 two-hit columns are evaluated.**

## Target

Let the existing certified q=2 boundary block be

\[
V_0=(|b_0\rangle,\ldots,|b_{31}\rangle)
\]

and let the already computed one-hit Euclidean columns be

\[
|g_{vi}\rangle=H^E_v|b_i\rangle,
\qquad v=0,\ldots,4.
\]

The next even constraint-generated family is

\[
\boxed{
|q_{wvi}\rangle=H^E_w H^E_v|b_i\rangle,
\qquad w,v=0,\ldots,4,\ i=0,\ldots,31.
}
\]

There are 800 labelled columns before rank reduction.

No expected rank, gap, zero count or projector overlap is frozen.

## Why this is the first enlarged even carrier relevant to the boundary heat problem

The q=2 all-j=1/2 boundary is even under doubled-spin parity. One Euclidean Hamiltonian action is odd. The positive Euclidean normal master

\[
M_E=\sum_w H_w^{E\dagger}H_w^E
\]

is parity even, so the master-Krylov sequence generated from `V0` remains in the even sector. Consequently the one-hit odd K1 carrier is useful operator data, but it does not itself enlarge the even boundary carrier on which

\[
V_0^\dagger e^{-\tau M_E}V_0
\]

lives.

The two-hit family `q_(wvi)` is the first actual representation-dressed even carrier that can overlap the original boundary and enlarge its finite trial habitat.

This statement does **not** identify the two-hit family with `M_E V0`, because the production master uses `H^dagger H`; no self-adjointness shortcut is assumed.

## Reuse rule

The 160 one-hit states MUST be loaded from the certified reusable Euclidean packet produced by workflow run `33970844680`. They are not recomputed inside every K2 shard.

Each independent shard is labelled by `(target_node=w, source_node=v)` and evaluates only the 32 second actions

\[
H_w^E |g_{vi}\rangle.
\]

The operator, `Jmax=5/2`, ordering and pruning convention are identical to the frozen Euclidean HH/K1-Ritz calculation.

## Rank diagnostics

Let `Q` denote the 800-column map made from the two-hit states. Define

\[
G_2=Q^\dagger Q.
\]

Let

\[
X=V_0^\dagger Q
\]

be the exact boundary-return block. Since `V0` is orthonormal, the Gram of the component orthogonal to the boundary is

\[
\boxed{
G_{2,\perp}=G_2-X^\dagger X.
}
\]

The measured number of genuinely new even directions beyond the boundary is

\[
r_{2,new}=\operatorname{rank}G_{2,\perp},
\]

and the total measured even carrier dimension through two constraint hits is

\[
\boxed{
\dim K_{even,\le2}=32+r_{2,new}.
}
\]

Rank is a scientific output, not a PASS target.

## Pass/fail

Implementation PASS requires only:

- exact 25 x 32 labelled coverage;
- all source one-hit packets pass their existing gates;
- all second-action amplitudes are finite;
- `Jmax=5/2` is respected;
- `G2` and `G2_perp` are Hermitian positive semidefinite within frozen numerical tolerance;
- boundary-return amplitudes reconstructed from sparse states agree with the combined Gram block;
- no missing/duplicate `(w,v,i)` labels.

## What this does not compute

This carrier does not by itself give

\[
\mu_2=V_0^\dagger M_E^2V_0,
\]

because `M_E` contains adjoints and the two-hit forward family is not silently identified with `M_E V0`.

It also does not emit `P_phys`, a Lorentzian master, an HDA certificate, a physical-time propagator, scalar dark matter, or FLRW dark energy.

Its purpose is narrower and physical: create the first reusable **even enlarged habitat data** needed by the already-frozen block-Krylov/projected-heat producer.
