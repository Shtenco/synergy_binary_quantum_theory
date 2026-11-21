# Spectral small-world import audit — 2026-09-06

Status: **methodological import audit; no external numerical coincidence is promoted to a BQG prediction.**

External article reviewed:

- https://habr.com/ru/articles/1046730/

Mathematical convergence reference used for the lattice-sum check:

- https://mathworld.wolfram.com/EpsteinZetaFunction.html

## 1. Executive verdict

The external small-world theory contains a genuinely useful mathematical toolbox, but its strongest numerical claims mix that toolbox with substantial model freedom and at least two derivational gaps. BQG should import the **spectral/operator machinery**, not the fitted or assumed numerical layer.

Useful imports:

```text
operator spectrum
heat kernel / heat trace
spectral measure
spectral zeta diagnostics
gap scaling
spectral dimension diagnostics
refinement / RG comparison
held-out numerical checks
```

Rejected imports:

```text
assume d=3 and then count the result as a derivation of d=3
fix K=6 from the same dimensional premise and reuse K in later “independent” checks
choose a Watts–Strogatz/small-world graph before deriving BQG connectivity
identify a master-constraint eigenvalue with a particle mass or physical frequency
fit graph-channel weights to known interaction scales
use target-scale agreement as a substitute for a zero-fit derivation
silently regularize a divergent lattice sum
replace zeta(4) by zeta(2) without a derived asymptotic step
```

The BQG replacement is stronger: derive the spectral graph from the actual positive BQG master operator through its cyclic block-Krylov quotient.

---

## 2. The main degrees-of-freedom problem in the article

The article explicitly starts with effective spatial dimension three, states an optimal local degree range, writes

\[
pK=N^{-1/3},
\]

and then selects

\[
K=6.
\]

It subsequently uses the same `d=3`, `K=6`, small-world ansatz and related spectral assumptions to obtain a value of `ln N` near 280 and to describe the numerical agreement of two formulae as an independent consistency check.

That is not a clean independence structure. Shared upstream assumptions induce correlated downstream agreement.

For BQG the firewall is:

\[
\boxed{\text{connectivity, spectral dimension and coordination must be outputs of the frozen operator/refinement construction.}}
\]

They may not be chosen to make a later physical number work.

---

## 3. Concrete lattice-sum problem

The article uses

\[
\sum_{\mathbf n\in\mathbb Z^3\setminus\{0\}}\frac{1}{|\mathbf n|^2}
=6\,\zeta(2).
\]

As an ordinary three-dimensional lattice sum, the left-hand side is divergent.

For a positive-definite quadratic form in `n` dimensions, the defining Epstein-zeta series converges in the half-plane

\[
\operatorname{Re}s>\frac n2.
\]

Here `n=3` and `s=1`, so

\[
1<\frac32.
\]

Therefore the displayed equality is not an equality of ordinary convergent sums.

An analytic continuation / zeta regularization can be introduced, but then it must be declared explicitly, with the continuation prescription and normalization derived and held fixed. It cannot be silently substituted for the divergent lattice sum.

BQG avoids this problem at finite regulator. For the derived finite quotient spectrum

\[
J_V u_a=\lambda_a u_a,
\]

define only the positive finite sum

\[
\boxed{
\zeta_V^{(+)}(s)
=\sum_{\lambda_a>0}\operatorname{tr}(W_a)\lambda_a^{-s}.
}
\]

No analytic continuation is needed until a genuine regulator-removal limit requires one.

---

## 4. Concrete `zeta(4) -> zeta(2)` gap

The article writes a nonlocal contribution of the form

\[
\sum_{n=1}^{\infty}\frac{p}{n^2}\,\ln K\,\frac1{n^2}
= p\ln K\,\zeta(4),
\]

and immediately afterward states that in the leading order in `p`

\[
\zeta_{\rm nonlocal}(1)\approx p\ln K\,\zeta(2).
\]

No transformation shown between those two equations converts `zeta(4)` into `zeta(2)`.

A valid replacement would require an additional derived measure, degeneracy, Jacobian, renormalization rule, or asymptotic resummation. Without that missing derivation the substitution is not accepted by BQG.

---

## 5. What the spectral apparatus becomes inside BQG

Let the actual regulated constraints be `C_A` and let

\[
\mathbb M_G=C_A^\dagger G^{AB}C_B\ge0.
\]

For a frozen seed block `V`, generate

\[
\mathcal K_r(\mathbb M_G,V)
=\operatorname{span}\{V,\mathbb M_GV,\ldots,\mathbb M_G^rV\}.
\]

The derived block graph is

\[
J_r=Q^\dagger\mathbb M_GQ.
\]

Its topology is not selected phenomenologically. Every new block appears only because the actual operator generates a new linearly independent direction.

The finite exact closure certificate is

\[
\boxed{B_{r+1}=0}
\]

or equivalently a certified zero next residual after full reorthogonalization.

Then

\[
\boxed{
V^\dagger f(\mathbb M_G)V
=E_0^\dagger f(J_r)E_0
}
\]

for every function on the finite spectrum, including

\[
e^{-\sigma\mathbb M_G}
\]

and the isolated zero spectral projector.

This is an operator theorem, not numerical pattern matching.

---

## 6. Constraint spectral dimension, not spacetime dimension

From the seed-weighted spectral measure one may define

\[
P_V(\sigma)
=\frac1{d_V}\sum_a e^{-\sigma\lambda_a}\operatorname{tr}W_a
\]

and

\[
\boxed{
d_{s,\mathbb M}^{(V)}(\sigma)
=-2\frac{d\ln P_V}{d\ln\sigma}.
}
\]

This is useful for comparing finite regulators and for detecting crossovers in the constraint spectrum.

The canonical name remains **constraint spectral dimension**.

Forbidden shortcut:

```text
constraint spectral dimension ~= 3
therefore physical spacetime dimension is derived
```

A spacetime interpretation requires a separate theorem tying the master spectral measure to the physical history effective kernel / derived geometry under refinement.

---

## 7. Current actual Euclidean BQG status after repository recovery

The recovered research line is a direct descendant of the current canonical `main` and contains the complete spectral-history production architecture.

Important already-computed finite facts include the 32-dimensional q=2 Euclidean boundary master at `Jmax=5/2`:

```text
rank      = 32
nullity   = 0
lambda_min = 9.651811183254074
lambda_max = 14.48385071910081
```

Therefore the bare 32D boundary carrier contains no common Euclidean master zero.

The next actual history calculation is not a guessed spectral model. It is

\[
Y_i=\mathbb M_E b_i,
\qquad i=0,\ldots,31,
\]

followed by

\[
\mu_1=V_0^\dagger Y,
\qquad
\mu_2=Y^\dagger Y,
\]

\[
\boxed{
R_1=\mu_2-\mu_1^\dagger\mu_1=B_1^\dagger B_1.
}
\]

The implementation to factor `R1 -> Q1,B1` exists and is preregistered.

However, as of this audit the historical GitHub Actions run did **not** complete all 32 raw master-image columns. Only one old artifact (`input_index=13`) was recoverable from that run. Hence the actual numerical rank of `Q1`, the actual `B1`, `A1`, and `B2` are not yet claimed.

This corrects any earlier wording that suggested that the production `Q1` numerical result had already been emitted.

---

## 8. Recovered real master-image column #13

The preserved production artifact has

```text
schema                         = BQG_EUCLIDEAN_MASTER_IMAGE_COLUMN_V1
input_index                    = 13
Jmax                           = 2.5
sparse output support          = 2273
output norm                    = 14.336696317044131
max spin reached               = 1.5
all sparse amplitudes finite   = true
spin cutoff respected          = true
state SHA256                   = ffbea26837aa99dac9951b0db7b143c9a4ac90480c9665db5ccfb5fa7740d918
```

This is a genuine microscopic datum generated by

\[
\mathbb M_E b_{13}=\sum_{v=0}^{4}(H_v^E)^\dagger H_v^E b_{13}
=\sum_{v=0}^{4}(H_v^E)^2b_{13}
\]

in the frozen explicitly Hermitian Euclidean convention.

It is not enough by itself to determine `mu2`; all 32 columns are required.

---

## 9. Exact next decision tree

The production spectral-history calculation is now frozen as follows.

### Step A — finish all 32 Euclidean master images

Compute

\[
Y=(Y_0,\ldots,Y_{31}).
\]

No post-hoc pruning threshold may be selected from the observed `R1` spectrum.

### Step B — calculate the first true spectral edge

\[
\mu_2=Y^\dagger Y,
\qquad
R_1=\mu_2-\mu_1^\dagger\mu_1.
\]

If `R1` is certified zero, the seed closes at depth zero.

If not, factor its positive support and emit actual `Q1,B1`.

### Step C — apply the same master to Q1

\[
Z_1=\mathbb M_EQ_1,
\qquad
A_1=Q_1^\dagger Z_1,
\]

\[
R_2=Z_1-Q_0B_1^\dagger-Q_1A_1.
\]

After full reorthogonalization:

- if `R2=0` within propagated numerical error, the two-block graph closes the finite Euclidean seed history exactly;
- otherwise factor `R2=Q2B2` and continue.

### Step D — regulator sequence

A finite `Jmax=5/2` closure is not continuum closure. Repeat the preregistered observables over increasing cutoffs/refinements and test normalized spectral/projector stability.

### Step E — full physical BQG master

The Euclidean reference master alone is insufficient for physical history. The complete positive master must incorporate the production constraint packet, including the Lorentzian/HDA requirements or an explicitly certified independent `D_target` construction on the same habitat.

### Step F — source-complete history

Only then use the source-complete seed and the existing legal order

\[
P_{\rm phys}\to Z[J]\to W[J]\to\Gamma[g]\to\Gamma^{(2)}_{\rm metric}.
\]

No master spectral parameter is relabeled as physical frequency.

---

## 10. Final scientific judgment

The external article is most valuable to BQG as a reminder that a theory can and should be interrogated through its **entire operator spectrum**, heat history, gap structure and refinement flow rather than a few hand-picked scalar coincidences.

BQG can use that apparatus in a substantially cleaner way because its spectral graph need not be postulated:

\[
\boxed{
\text{actual constraints}
\to\mathbb M_{BQG}
\to\text{derived Krylov graph}
\to d\Sigma_V
\to\text{finite exact history/projector}
}
\]

when the residual actually closes.

The strongest defensible goal is therefore not to reproduce the external article's `N`, masses or constants. It is to make the BQG finite and refinement histories **spectrally complete, zero-fit, source-complete, and falsifiable**.
