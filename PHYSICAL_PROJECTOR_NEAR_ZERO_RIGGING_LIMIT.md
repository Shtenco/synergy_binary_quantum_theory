# Physical projector from an asymptotic near-zero master sector

Status: **operator-limit theorem/protocol for the physicalization frontier; no BQG near-zero sector is asserted by this document.**

## 1. Why exact zero at every finite regulator is too strong

For an exactly anomaly-free finite constraint representation, the master constraint may possess an exact zero eigenspace and

\[
P_0=\lim_{\tau\to\infty}e^{-\tau\mathbb M}
\]

is immediate.

In a refinement programme, however, the regulated quantum constraint algebra can close only asymptotically. A physical continuum zero sector can then appear as eigenvalues that approach zero with the regulator even if no exact finite-regulator zero exists.

Therefore the implication

```text
finite M_epsilon has lambda_min > 0
-> physical sector is empty
```

is not valid without a regulator-limit analysis.

At the same time, choosing an arbitrary low-eigenvalue window after seeing the spectrum is not acceptable.

## 2. Spectral-separation criterion

Let the positive regulated master have ordered eigenvalues

\[
0\le\lambda_1(\epsilon)\le\cdots\le\lambda_N(\epsilon).
\]

A rank-`r` asymptotic zero sector is a candidate only if a preregistered refinement sequence satisfies at least

\[
\boxed{
\lambda_r(\epsilon)\to0
}
\]

and

\[
\boxed{
\frac{\lambda_r(\epsilon)}{\lambda_{r+1}(\epsilon)}\to0.
}
\]

The second condition is the important one: it provides a scale separation between a collapsing low cluster and the rest of the spectrum.

A merely small eigenvalue with no increasing separation is not a physical-projector certificate.

## 3. A non-arbitrary heat-kernel schedule

Assume for the moment

\[
0<\lambda_r<\lambda_{r+1}.
\]

Choose

\[
\boxed{
\tau_\epsilon
=\frac1{\sqrt{\lambda_r(\epsilon)\lambda_{r+1}(\epsilon)}}.
}
\]

Then

\[
\tau_\epsilon\lambda_r
=\sqrt{\frac{\lambda_r}{\lambda_{r+1}}}
\to0,
\]

while

\[
\tau_\epsilon\lambda_{r+1}
=\sqrt{\frac{\lambda_{r+1}}{\lambda_r}}
\to\infty.
\]

For the low cluster,

\[
\max_{i\le r}|e^{-\tau\lambda_i}-1|
\le 1-e^{-\tau\lambda_r},
\]

and for the high sector

\[
\max_{i>r}e^{-\tau\lambda_i}
\le e^{-\tau\lambda_{r+1}}.
\]

Hence, if `P_r(epsilon)` denotes the spectral projector onto the first `r` modes,

\[
\boxed{
\left\|e^{-\tau_\epsilon\mathbb M_\epsilon}-P_r(\epsilon)\right\|
\le
\max\left(
1-e^{-\sqrt{\lambda_r/\lambda_{r+1}}},
\,e^{-\sqrt{\lambda_{r+1}/\lambda_r}}
\right)
\to0.
}
\]

The heat-kernel separation therefore follows from a measured eigenvalue ratio rather than an arbitrary spectral threshold.

If the low sector is already exactly zero at finite regulator, revert to the exact-gap theorem with `tau Delta -> infinity`.

## 4. Eigenvalue collapse is not enough

Even with the spectral ratio above, the physical subspace itself must stabilize under the declared refinement embeddings.

Let

\[
I_{\epsilon\to\epsilon'}
\]

be the frozen comparison/embedding map between regulator Hilbert spaces. Then a physical-projector claim requires a Cauchy-type condition such as

\[
\boxed{
\left\|
P_r(\epsilon')
-I P_r(\epsilon)I^\dagger
\right\|_{relevant}
\to0.
}
\]

The precise norm/domain must be frozen for the production refinement family.

A low eigenspace that rotates without convergence is not a continuum physical Hilbert space merely because its eigenvalues approach zero.

## 5. Boundary overlap is the relevant observable

For a frozen boundary block `B_epsilon`, define

\[
G_r(\epsilon)
=B_\epsilon^\dagger P_r(\epsilon)B_\epsilon.
\]

A useful physical boundary sector requires

\[
\boxed{
G_r(\epsilon)\to G_{phys}
}
\]

with stable nonzero support.

The corresponding heat approximation is

\[
G_{\tau}(\epsilon)
=B_\epsilon^\dagger e^{-\tau_\epsilon\mathbb M_\epsilon}B_\epsilon.
\]

The same schedule must stabilize source-dressed matrices

\[
A_{a,\tau}
=B^\dagger e^{-\tau\mathbb M/2}O_a e^{-\tau\mathbb M/2}B.
\]

It is not sufficient for `G_tau` alone to look stable while source matrices drift.

## 6. Rank selection must not be post-hoc

The integer `r` may not be chosen simply because one observed spectrum has an attractive gap.

A production claim must instead use one of the following fail-closed procedures:

1. `r` is fixed by an independently derived symmetry/constraint-sector prediction before the refinement scan; or
2. the entire ordered set of adjacent ratios `lambda_i/lambda_{i+1}` is reported on every refinement, and a persistent isolated cluster is treated as a discovered candidate that must then survive new held-out refinement levels without changing the selection rule.

The first refinement levels can discover a candidate rank; later levels must act as held-out tests.

## 7. Relation to the existing joint cutoff theorem

For the fixed-input Euclidean HH family, the repository already proves that the Peter-Weyl cutoff error is exactly zero once

\[
J_{max}\ge5/2.
\]

Thus a near-zero trend in that declared family cannot be blamed on further `Jmax` truncation above the safe wall.

For the full Lorentzian HH support the conservative safe wall is larger and must be respected independently.

For genuinely growing collective/refinement input spins, a separate schedule

\[
J_{max}(b)\ge j_{in}(b)+r_{hits}/2
\]

must accompany the refinement limit.

## 8. Production diagnostics

For every regulator/refinement level report, without hiding failed levels:

- all low master eigenvalues relevant to the candidate cluster;
- `lambda_r/lambda_{r+1}`;
- chosen `tau_epsilon` from the frozen formula;
- heat-to-spectral-projector error;
- boundary Gram `B^dag P_r B` and `B^dag exp(-tau M) B`;
- principal angles / projector distance between embedded low subspaces on successive refinements;
- source-dressed operator convergence;
- sensitivity to the positive master metric `G^{AB}`.

## 9. Falsifiers

The asymptotic physical-projector interpretation fails if any of the following persists:

1. `lambda_r` does not approach zero;
2. `lambda_r/lambda_{r+1}` does not decrease toward zero;
3. the discovered rank changes on held-out refinement without a derived reason;
4. the low eigenspace does not converge under the frozen embeddings;
5. the boundary overlap tends to zero;
6. source matrices do not converge under the same sequence;
7. the limiting low sector changes when only the positive master metric `G^{AB}` is changed;
8. the effect disappears once the Peter-Weyl cutoff is raised above its proven safe wall.

## 10. Scientific boundary

This protocol allows a continuum physical sector to emerge from asymptotically vanishing regulator defects without demanding an exact zero at every finite step. It does not license fitting a spectral window to cosmological data and does not by itself establish any BQG physical state, graviton, dark matter or dark energy.
