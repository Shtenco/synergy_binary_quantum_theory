# Background absolute normalization and the vacuum-energy observable

Status: **physicalization normalization firewall. No BQG cosmological constant is derived here.**

## 1. The central distinction

Two different uses of a generating functional require different normalization discipline.

### Connected fluctuations on one fixed background

For normalized expectation values and connected correlators it is natural to use a normalized physical density matrix on a fixed background,

\[
\operatorname{Tr}\rho_{phys}=1.
\]

Then

\[
Z_{corr}[J]
=\operatorname{Tr}\rho_{phys}e^{J\cdot O}
\]

satisfies `Z_corr[0]=1` by construction.

This is appropriate for susceptibilities such as

\[
\frac{\delta^2W}{\delta J_A\delta J_B}.
\]

### Homogeneous gravitational background / vacuum energy

The zero-source gravitational effective action depends on the **relative absolute physical amplitude of different background geometries**.

If one separately imposes

\[
Z[g;0]=1
\]

for every geometry `g`, all geometry-dependent zero-source free-energy information is erased by definition.

That procedure cannot be used to determine a cosmological term.

## 2. Correct background object

Let

\[
\mathcal Z_{phys}[g]
\]

be the theory-specific physical boundary/history amplitude with one regulator and measure normalization convention held fixed across the declared background family.

Only a geometry-independent overall multiplicative constant is physically irrelevant at this stage.

A safe observable is an amplitude ratio to a frozen reference geometry `g_ref`:

\[
\boxed{
\Delta W[g;g_{ref}]
=-\log\frac{\mathcal Z_{phys}[g]}{\mathcal Z_{phys}[g_{ref}]}
}
\]

for a Euclideanized/positive amplitude convention, or the corresponding phase/effective-action relation in the Lorentzian construction.

A geometry-independent normalization cancels in the ratio, while a geometry-dependent extensive contribution survives.

## 3. Why this matters for a cosmological term

Under a common length rescaling

\[
\ell\to\lambda\ell,
\]

the four-dimensional Einstein-Hilbert curvature contribution scales schematically as

\[
S_{EH}\sim\lambda^2,
\]

while a cosmological four-volume contribution scales as

\[
S_{\Lambda}\sim\lambda^4.
\]

Therefore a background amplitude of the form

\[
\mathcal Z(\lambda)
\propto
\exp[-c_2\lambda^2-c_4\lambda^4+\cdots]
\]

contains distinguishable curvature-like and vacuum-volume-like components:

\[
\boxed{
-\log\mathcal Z
=c_2\lambda^2+c_4\lambda^4+\cdots.
}
\]

Normalizing `Z(lambda)` independently to unity at every `lambda` destroys both terms and therefore cannot constitute evidence that `c4=0`.

## 4. Extensive history / refinement version

For a homogeneous family with `N_4` equivalent history cells, a vacuum-like contribution should be studied through an extensive limit such as

\[
\boxed{
f_{hist}
=-\lim_{N_4\to\infty}
\frac1{N_4}
\log\mathcal Z_{phys}(N_4)
}
\]

or its refinement-corrected analogue.

If the physical four-volume per cell is derived,

\[
V_4=N_4v_4,
\]

then an extensive constant contribution supplies a candidate vacuum energy density only after the common physical scale and measure normalization are fixed.

The raw number `f_hist` is not yet an observed cosmological constant.

## 5. Projector normalization caveat

Rigging maps and heat-kernel projectors can carry regulator-dependent normalization factors.

A legal renormalization may remove a factor that is independent of the physical background geometry and source configuration.

It may **not** remove a factor whose logarithm scales with physical four-volume or another background invariant, because that is precisely the information the gravitational effective action is meant to determine.

Therefore every future projector normalization must be classified as

```text
regulator-only / geometry-independent -> removable convention
geometry-dependent extensive factor   -> physical effective-action candidate
source-dependent factor                -> physical response candidate
```

before normalization.

## 6. Relation to the local j=1 volume source

The exact local collective-volume positive control has

\[
\Gamma_V(p_0)=0
\]

at its normalized zero-source state.

This is a Legendre-transform normalization statement for a finite normalized trace. It is **not** evidence for zero vacuum energy and must never be compared to the cosmological constant.

The local volume susceptibility and the global zero-source background free-energy density are different observables.

## 7. FLRW production target

The homogeneous physicalization branch should preserve an unnormalized/ratio-valued object

\[
\mathcal Z_{phys}[a,N]
\]

with one frozen measure convention and form

\[
\Gamma_{FLRW}[a,N]
\]

up to a geometry-independent additive constant only.

Then lapse and scale-factor variations determine the effective background density and pressure after the gauge/source ordering rules are respected.

Only after

\[
\rho_{hist}(a)
\]

is derived may one infer

\[
w_{hist}(a)
=-1-\frac13\frac{d\ln\rho_{hist}}{d\ln a}.
\]

## 8. Strong falsifiers

A claimed BQG dark-energy term is rejected if

1. it disappears or changes arbitrarily under one fixed geometry-independent normalization convention;
2. it exists only because each background was assigned a different measure normalization;
3. the inferred vacuum contribution is actually a finite normalized-trace constant such as `Gamma_V(p0)=0`;
4. a supposed `lambda^4` term is not stable under refinement/extensive scaling;
5. the same physical history measure used for scalar/TT observables cannot reproduce the background amplitude;
6. an evolving `w(a)` is fitted before `rho_hist(a)` is derived.

## 9. Claim boundary

This rule preserves the observable that could contain a cosmological-volume term. It does not predict its sign, magnitude or even its existence. The current canonical statement remains: the physical BQG cosmological constant/dark-energy sector is open.
