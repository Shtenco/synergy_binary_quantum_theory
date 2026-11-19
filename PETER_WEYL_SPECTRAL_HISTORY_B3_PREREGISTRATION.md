# Peter–Weyl spectral-history B3 preregistration

Status: **next exact finite Euclidean shell preregistered; no closure outcome assumed.**

This calculation continues the already certified parity block-Lanczos chain for

\[
H=H^{sine}_{E,0}+H^{sine}_{E,1}
\]

on the complete 32-dimensional all-\(j=1/2\) logical K5 seed sector.

It is distinct from the positive-master route based on
\(\mathbb M=\sum_v H_v^\dagger H_v\).  The present calculation asks a narrower and cheaper question: **does the exact finite two-vertex Euclidean constraint Krylov chain close after the two hopping blocks already measured in the repository?**

## 1. Frozen previous evidence

The certified higher-shell artifact from workflow run `31852849936` gives

\[
P H P=0,
\]

\[
B_1^\dagger B_1=K=P H^2P,
\]

and

\[
B_2^\dagger B_2=\Lambda
=K^{-1/2}(P H^4P-K^2)K^{-1/2}.
\]

The recorded reconstruction errors are approximately

```text
||B1^dag B1-K||      = 1.616e-13
||B2^dag B2-Lambda|| = 1.645e-13
```

with both `K` and `Lambda` full rank on the 32-dimensional seed block.

The source artifacts contain every exact sparse column

\[
a_i=H|i\rangle,
\qquad
b_i=H^2|i\rangle.
\]

Thirty-one columns come from source run `31844567559`; column 28 was independently repaired in run `31852849936`.  The B3 extension records a SHA256 for each inherited source column before acting on it.

## 2. One additional operator application

For every seed basis vector compute only

\[
\boxed{c_i=H b_i=H^3|i\rangle.}
\]

Then the new moment is

\[
\boxed{H_6=P H^6P=C^\dagger C}
\]

where `C` is the column map formed from the `c_i`.

No large Ritz basis is required.

## 3. Exact B3 identity

Use the parity block recurrence

\[
H Q_0=Q_1B_1,
\]

\[
H Q_1=Q_0B_1^\dagger+Q_2B_2,
\]

\[
H Q_2=Q_1B_2^\dagger+Q_3B_3.
\]

Because the parity grading kills diagonal Lanczos blocks in this chain,

\[
H^2Q_0=Q_0K+Q_2B_2B_1.
\]

The component of `H^3 Q0` along `Q1` is

\[
C_1=Q_1^\dagger H^3Q_0
=B_1^{-1}H_4,
\qquad H_4=P H^4P.
\]

A mandatory internal identity is

\[
\boxed{B_1^{-1}H_4=B_1K+\Lambda B_1.}
\]

Subtract this already-explained component:

\[
R_3=H^3Q_0-Q_1C_1.
\]

Then

\[
R_3^\dagger R_3
=H_6-C_1^\dagger C_1.
\]

Since

\[
R_3=Q_3B_3B_2B_1,
\]

and `B1`, `B2` are full rank, the next hopping Gram is exactly

\[
\boxed{
B_3^\dagger B_3
=(B_2B_1)^{-\dagger}
\left(H_6-C_1^\dagger C_1\right)
(B_2B_1)^{-1}.
}
\]

This formula is evaluated by `scripts/peter_weyl_spectral_history_b3_gate.py`.

## 4. Preregistered stopping rule

The outcome is not chosen in advance.

### Closure

If the resolved rank is zero,

\[
\boxed{\operatorname{rank}(B_3^\dagger B_3)=0,}
\]

then this finite Euclidean two-vertex seed Krylov history closes at depth 2.

### Non-closure

If

\[
\operatorname{rank}(B_3^\dagger B_3)>0,
\]

then the finite history is **not** closed at depth 2.  This is not a failed calculation: it is evidence that the operator generates a genuinely new shell and the next exact hopping must be computed.

The code therefore separates `passed` (algebra/provenance/regulator checks) from `science_status` (closed versus open after B3).

## 5. Frozen numerical rank rule

Let

\[
S=\max(\|\operatorname{spec}(B_3^\dagger B_3)\|_\infty,1).
\]

The resolved numerical rank uses

\[
\lambda>10^{-10}S.
\]

Negative eigenvalues beyond the separately frozen floating-point PSD tolerance fail the gate rather than being clipped into a desired closure result.

## 6. Regulator guard

The inherited calculation used doubled-spin cutoff `Jmax2=5`, i.e.

\[
J_{max}=5/2.
\]

The third-hit output is accepted as cutoff-safe only if its observed maximum spin is at most `2.0`.  Thus the third action remains one half-unit below the truncation wall.  If the output reaches the wall, the result is not promoted as a regulator-safe B3 certificate.

## 7. Claim boundary

Even a zero `B3` would establish only exact closure of this **finite Euclidean two-vertex constraint Krylov history**.  It would not establish:

- the full positive-master physical zero-sector history;
- the complete Lorentzian/global constraint family;
- regulator removal or the continuum rigging map;
- physical time/frequency;
- connected `W[J]`;
- gauge-reduced `Gamma[g]`;
- a physical graviton propagator or cosmological prediction.

Those remain separate fail-closed physicalization gates.
