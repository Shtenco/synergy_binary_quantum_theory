# Physical scalar IR discriminators

Status: **exact interpretation guard for future derived kernels; not a dark-sector model**.

This document freezes the distinctions that must be applied after a theory-specific physical scalar Hessian is finally obtained.

The target object is the **gauge/constraint-reduced physical-history kernel**

\[
K_{\rm red}(\omega,\mathbf k),
\]

not an unreduced lapse/shift matrix, a constraint resolvent in `z`, a smoothing exponent, a static local source Hessian, or a fitted cosmological response function.

---

## 1. Three different infrared structures

### A. Analytic noncritical response

If a scalar eigenchannel has

\[
K(0,k)=K_0+K_2k^2+O(k^4),
\qquad K_0\neq0,
\]

then

\[
K^{-1}(0,k)
=\frac1{K_0}+O(k^2).
\]

It is finite at `k=0` and does not by itself generate a long-range `1/k^2` response.

### B. Constraint/Poisson long-range response

If after legitimate constraint reduction a source-response channel has

\[
K(0,k)=Z_k k^2+O(k^4),
\qquad Z_k\neq0,
\]

then

\[
\boxed{K^{-1}(0,k)\sim\frac1{Z_k k^2}}.
\]

This is the structural behavior required for a long-range Poisson-like response. It does **not** by itself imply a new propagating particle: in GR the scalar Newtonian potentials are constrained response variables.

### C. Additional propagating scalar branch

A genuine extra scalar degree of freedom requires a physical zero of the frequency-dependent reduced kernel, for example

\[
K(\omega,k)
=Z_\omega\omega^2-Z_k k^2-m^2+\cdots
\]

with

\[
\boxed{K(\omega_*(k),k)=0}
\]

on an allowed physical sheet after gauge reduction. A dark-matter interpretation additionally requires the correct residue, stability, sound speed, abundance and universal gravitational/lensing response.

These three cases must never be conflated.

---

## 2. Why observer smoothing is not a dark-matter propagator

`OBSERVER_SCALE_SMOOTHING.md` obtains approximately

\[
\delta g_{RMS}\sim b^{-2}
\]

from

\[
N(b)\sim b^4,
\qquad
\delta g_{RMS}\sim N^{-1/2},
\]

under the explicit assumption of zero-mean, sufficiently weakly correlated microscopic contributions.

That result is a coarse-reconstruction/self-averaging law. It is not a physical connected scalar two-point function. The same document explicitly notes that long-range correlations can change the exponent and separately rejects the shortcut from smoothing to a quantum vacuum power spectrum.

Therefore no inference

```text
b^-2 smoothing -> 1/k^2 gravity -> dark matter
```

is allowed.

A long-range scalar response must instead arise from the actual physical projector/history kernel, for example through a constraint zero mode, a critical collective channel, or an additional physical pole.

---

## 3. Pure-GR degree-count guard

The repository's exact canonical discriminator starts from the connection/triad phase space with 18 phase dimensions per spatial point. For ordinary GR,

\[
18-2(3_G+3_D+1_H)=4
\]

physical phase dimensions, hence

\[
\boxed{N_{config}^{GR}=2}.
\]

These are the two local gravitational tensor/helicity configuration modes.

Consequently, if the future BQG physical reduction has the same independent first-class rank/reducibility structure as pure GR, a new propagating local scalar gravitational degree of freedom cannot simultaneously be claimed without explaining where the canonical count changed.

If an extra scalar pole survives physical reduction, the project must explicitly identify one of the following:

1. a changed independent constraint rank/reducibility structure;
2. an added collective physical field/carrier not contained in pure metric GR;
3. an effectively nonlocal history degree of freedom whose canonical counting differs from the local ADM count.

A pole extracted before this accounting is not accepted as emergent dark matter.

---

## 4. Modified gravity does not require a new scalar particle

It is possible for the physical source response to differ from GR while retaining only the GR propagating tensor degree count.

If after reduction

\[
\mu_{\rm BQG}(a,k)\neq1
\quad\text{or}\quad
\Sigma_{\rm BQG}(a,k)\neq1
\]

but there is no additional stable physical pole, the correct interpretation is a modified gravitational constraint/response law.

This possibility is especially important for dark-matter phenomenology: an enhanced nonrelativistic potential without a new pole is not particle dark matter. It must still pass the same lensing, growth, cluster and early-universe tests with one universal metric/source coupling.

---

## 5. Dark-energy guard

HDA closure does not determine the cosmological constant. In the repository's ADM family the cosmological term cancels from the hypersurface-deformation bracket, while `c_DW=1/2` and the relative kinetic/curvature normalization are selected.

The existing Regge/EH bridge tests the curvature action class. Under uniform four-dimensional length scaling,

\[
\int\sqrt g R\to\lambda^2\int\sqrt g R,
\]

whereas a cosmological-volume term scales as

\[
\int\sqrt g\to\lambda^4\int\sqrt g.
\]

Therefore neither successful HDA closure nor the existing Regge curvature reconstruction fixes a vacuum-energy coefficient.

A BQG dark-energy claim requires the **same physical history effective action** to generate a homogeneous volume-like contribution in `Gamma_FLRW[a,N]`. Only then may lapse and scale-factor variations be used to derive `rho_hist(a)` and `p_hist(a)`.

---

## 6. Frozen decision tree for a future scalar result

Given a derived physical `K_red(omega,k)`:

```text
Is the calculation after physical projection and gauge/constraint reduction?
  no  -> no physical scalar interpretation
  yes -> continue

Does static source response contain a k^-2 inverse-kernel law?
  no  -> no long-range Poisson-like scalar response in that channel
  yes -> long-range gravitational response exists; continue

Is there an additional frequency-dependent physical pole?
  no  -> modified/GR-like constraint response, not new dark particle
  yes -> audit canonical degree count and residue/stability

Does the extra pole survive regulator/refinement and have positive physical residue?
  no  -> reject dark-mode interpretation
  yes -> test sound speed, clustering, abundance and lensing

Does homogeneous Gamma_FLRW contain an independently generated volume-like term?
  no  -> no derived dark-energy density from this calculation
  yes -> derive rho_hist(a), p_hist(a), w_hist(a) before opening cosmology data
```

---

## 7. Falsifiable infrared outputs

Once the production history calculation exists, the minimal zero-fit scalar output should include

\[
\boxed{
\begin{aligned}
&K_{red}(\omega,k),\\
&\text{all physical scalar pole locations and residues},\\
&\lim_{k\to0} k^2 G_{\Psi T}(0,k),\\
&\lim_{k\to0} k^2 G_{(\Phi+\Psi)T}(0,k),\\
&\mu_{\rm BQG}(a,k),\\
&\Sigma_{\rm BQG}(a,k),\\
&\rho_{hist}(a),\ p_{hist}(a),\ w_{hist}(a),\\
&\text{constraint rank/reducibility by refinement scale}.
\end{aligned}}
\]

All of them must be generated from one frozen physical projector/history and one metric/source normalization. None may be selected from observational data after the fact.
