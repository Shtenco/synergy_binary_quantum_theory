# BCQG collective HDA intertwining theorem

**Status:** exact finite-dimensional compression identity + conditional collective-HDA inheritance theorem.  This is C3 infrastructure, not a claim that the production collective residual has already been measured.

## 1. Exact algebra identity under an isometric block map

Let

\[
W:\mathcal H_{eff}\to\mathcal H,
\qquad W^\dagger W=I,
\qquad P=WW^\dagger,
\qquad Q=I-P.
\]

For any two microscopic operators `A,B` define their direct compressed operators

\[
A_{eff}=W^\dagger A W,
\qquad
B_{eff}=W^\dagger B W.
\]

Insert `I=P+Q` in the microscopic commutator.  Exactly,

\[
W^\dagger[A,B]W
=[A_{eff},B_{eff}]
+W^\dagger A Q B W
-W^\dagger B Q A W.
\]

Therefore

\[
\boxed{
[A_{eff},B_{eff}]-W^\dagger[A,B]W
=-W^\dagger A Q B W+W^\dagger B Q A W.
}
\]

For Hermitian `A,B`, define the operator leakage norms

\[
\eta_A=\|QAW\|,
\qquad
\eta_B=\|QBW\|.
\]

Since

\[
W^\dagger A Q=(QAW)^\dagger,
\]

the exact identity implies the deterministic bound

\[
\boxed{
\|[A_{eff},B_{eff}]-W^\dagger[A,B]W\|
\le2\eta_A\eta_B.
}
\]

Thus an exactly invariant carrier inherits the algebra exactly, while a nearly invariant carrier has a quantitatively controlled algebra error.  No GR target enters this result.

## 2. Application to the HDA

Let the microscopic/regulator operators obey

\[
[H[N],H[M]]
=i\hbar D[\beta_{micro}(N,M)]+E_{micro}[N,M],
\]

where

\[
\beta_{micro}^a
=Q_{micro}^{ab}(N\partial_bM-M\partial_bN).
\]

Define direct compressed generators

\[
H_{eff}[N]=W^\dagger H[N]W,
\qquad
D_{proj}[\beta]=W^\dagger D[\beta]W.
\]

The collective metric is independently measured from the BCQG flux response.  Let

\[
\beta_{eff}^a
=Q_{eff}^{ab}(N\partial_bM-M\partial_bN)
\]

and let `D_eff[beta_eff]` be the held-out collective diffeomorphism target built from that coarse metric.

Define three independent defects:

### microscopic projected residual

\[
\delta_{micro}
=\|W^\dagger E_{micro}W\|;
\]

### compression/intertwining leakage

\[
\delta_{comp}
=\|[H_{eff}[N],H_{eff}[M]]-W^\dagger[H[N],H[M]]W\|;
\]

### structure-function blocking defect

\[
\delta_{str}
=\|D_{proj}[\beta_{micro}]-D_{eff}[\beta_{eff}]\|.
\]

Then, by the triangle inequality and the exact compression identity,

\[
\boxed{
\|[H_{eff}[N],H_{eff}[M]]-i\hbar D_{eff}[\beta_{eff}]\|
\le
\delta_{micro}+2\eta_N\eta_M+\hbar\delta_{str}.
}
\]

This is the required microscopic-to-collective HDA bridge.

A collective HDA fixed point therefore follows if, in the frozen refinement family,

\[
\frac{\delta_{micro}}{\|D_{eff}\|}\to0,
\qquad
\frac{\eta_N\eta_M}{\|D_{eff}\|}\to0,
\qquad
\frac{\hbar\delta_{str}}{\|D_{eff}\|}\to0.
\]

The preregistered direct collective bracket remains the final held-out check; this theorem explains exactly which mechanisms can make it fail.

## 3. Why the operator-first route survives a matrix-valued coarse metric

The production route implementation is not restricted to a scalar expectation-value metric.  It accepts an arbitrary finite matrix-valued positive geometry block `Q^{ab}` and defines

\[
\Omega(k)=\sqrt{Q^{ab}k_ak_b}/\epsilon,
\qquad
R[N]=\frac12\{N,\Omega\}.
\]

No pairwise commutativity of the matrices `Q^{ab}` is assumed.  For each momentum the Hermitian matrix

\[
A(k)=Q^{ab}k_ak_b
\]

is diagonalized spectrally and its positive square root is used before any geometry expectation value is taken.

Therefore C2 may provide a noncommuting finite coarse geometry block without changing the definition of the route-normal operator.  What C3 must test is the coarse structure-function/intertwining defect, not a fictitious commutativity assumption.

## 4. Hermitian v1.3 geometry and the old HDA scaling theorem

The physical v1.3 local geometry is

\[
\boxed{G_v=-\frac23E_v-\frac{32}{9}S_v},
\qquad
S_v=-\frac i2(L_{raw,v}-L_{raw,v}^\dagger).
\]

The fixed-cutoff HDA composition estimates depend on locality/boundedness of the finite geometry block and on the operator-first route scaling, not on the historical shorthand `+(32i/9)L_raw` being valid globally.  Hermitian completion is a bounded linear antisymmetrization at fixed safe cutoff, so the previously established hierarchy remains the relevant conditional scaling architecture:

\[
C_{cross}/D=O(\epsilon),
\qquad
C_{GG}/D=O(\epsilon^2),
\]

with the declared simultaneous-cutoff path supplying the earlier positive exponents when its support assumptions hold.

This statement does **not** replace the corrected Lorentzian V2 amplitude/covariance rerun (C1), because the finite constants and actual blocked Krylov space depend on the corrected operator.

## 5. C3 reduced certificate

After C1 and a C2 Schur/gap fixed point, C3 no longer requires treating every old IR target as an unrelated numerical miracle.  It requires four measured quantities on held-out lapses/perturbations:

1. projected microscopic HDA residual `delta_micro`;
2. Hamiltonian carrier leakages `eta_N,eta_M`;
3. structure-function blocking defect `delta_str` using the independently measured coarse metric;
4. direct collective bracket residual as an overdetermined check of the bound.

Additionally verify:

- the three shifts generated by held-out lapse pairs act independently on the coarse metric;
- the scalar Hamiltonian action/gradient is nonzero on the regular carrier;
- C2 leaves no unclassified gapless sector;
- the projected constraint bracket matrix has no second-class remainder.

If these quantities converge with refinement, the restricted IR universality theorem forces the leading ADM/Einstein Hamiltonian up to `G,Lambda`.

## 6. What this does and does not close

This theorem removes a conceptual gap: **HDA inheritance is controlled by leakage and structure-function matching rather than assumed.**

It does not turn a queued production run into data.  Internal gravitational closure still requires green C1, production C2, and the measured C3 residuals above.
