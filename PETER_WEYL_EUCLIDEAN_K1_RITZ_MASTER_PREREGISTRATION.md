# Preregistration: Euclidean master Ritz spectrum on the first constraint-generated Peter-Weyl habitat

Status: **frozen before the K1 Gram rank and before any K1 master Ritz eigenvalue is observed.**

## 1. Scientific question

The bare q=2 boundary carrier is not invariant under the Euclidean node constraints. The next finite question is therefore not whether a zero vector exists in the original 32 states, but whether allowing the first genuine representation-dressed layer lowers the positive normal-constraint master toward a zero sector.

Define

\[
\mathcal K_1
=\operatorname{span}\{B,H_v^EB\}_{v=0}^4.
\]

The present calculation computes the Euclidean master Ritz spectrum on this measured subspace without constructing a dense ambient Peter-Weyl Hamiltonian.

No expected minimum eigenvalue, rank, nullity or gap ratio is frozen.

## 2. First generated layer

Use the 160 labelled outgoing columns

\[
|g_\alpha\rangle
=H_v^E|b_i\rangle,
\qquad
\alpha=(v,i),
\]

and their Gram

\[
\boxed{
G_{\alpha\beta}
=\langle g_\alpha|g_\beta\rangle.
}
\]

Let

\[
G=U\Lambda U^\dagger
\]

and retain only the numerical support selected by the same frozen relative Gram tolerance used in the K1 span gate.

An orthonormal basis for the generated odd layer is represented algebraically by

\[
Q_1=A U_r\Lambda_r^{-1/2},
\]

where `A` is the column map whose columns are the `g_alpha`.

## 3. Euclidean master on the odd K1 layer

The identity-metric Euclidean normal master is

\[
\mathbb M_E
=\sum_{w=0}^4H_w^{E\dagger}H_w^E.
\]

Define

\[
\boxed{
D_{\alpha\beta}
=\langle g_\alpha|\mathbb M_E|g_\beta\rangle
=\sum_{w=0}^{4}
\langle H_w^Eg_\alpha|H_w^Eg_\beta\rangle.
}
\]

The exact Ritz matrix on the orthonormal generated layer is

\[
\boxed{
M_{odd}^{K1}
=\Lambda_r^{-1/2}U_r^\dagger D U_r\Lambda_r^{-1/2}.
}
\]

Equivalently its eigenvalues solve the generalized problem

\[
D c=\lambda G c
\]

on `supp(G)`.

## 4. Parity block decomposition

The original boundary `B` is even under the doubled-spin grading and one Euclidean action is odd. Since

\[
\mathbb M_E=\sum_wH_w^\dagger H_w
\]

is even,

\[
\boxed{
\langle B|\mathbb M_E|Q_1\rangle=0.
}
\]

Therefore the K1 Ritz master is block diagonal:

\[
\boxed{
M_E|_{\mathcal K_1}
=
M_B\oplus M_{odd}^{K1}.
}
\]

The full K1 Ritz spectrum is simply the union of the bare-boundary five-node master spectrum and the generated-odd spectrum.

## 5. Regulator wall

Every matrix element of `D` contains two Euclidean Hamiltonian actions from an all-`j=1/2` input.

The exact K5 HH reachability theorem freezes

\[
\boxed{J_{max}=5/2}
\]

as sufficient and saturated for this Euclidean two-H family. No larger cutoff is chosen from the forthcoming Ritz spectrum.

The current reference Peter-Weyl implementation uses amplitude pruning at `1e-8`. This is retained for the first pass.

**Critical rule:** any apparent Ritz eigenvalue near the numerical/pruning floor must be rerun with tighter amplitude thresholds before being interpreted as a physical near-zero mode. The first pass may discover a candidate; it may not certify continuum zero physics at its own numerical floor.

## 6. Sharded calculation

For each target node `w`, independently compute

\[
D^{(w)}_{\alpha\beta}
=\langle H_wg_\alpha|H_wg_\beta\rangle.
\]

Then

\[
\boxed{D=\sum_{w=0}^{4}D^{(w)}.}
\]

Each shard also recomputes `G` from the same 160 first-layer columns. The aggregator requires the five independently produced `G` matrices to agree within the frozen numerical tolerance before summing the second-layer master blocks.

This duplication is computationally expensive but gives a useful deterministic cross-run regression for the first-layer data.

## 7. Outputs frozen before result

Report:

- `rank(G)` and generated odd dimension;
- cross-shard `G` consistency error;
- PSD/Hermiticity errors of every `D^(w)` and summed `D`;
- eigenvalues of `M_odd^K1`;
- bare-boundary master eigenvalues reconstructed from diagonal blocks of `G`;
- union K1 spectrum;
- minimum boundary eigenvalue `lambda_B,min`;
- minimum generated-layer Ritz eigenvalue `lambda_K1,odd,min`;
- full K1 minimum;
- ratio
  \[
  r_{dress}=\lambda_{K1,min}/\lambda_{B,min};
  \]
- per-target-node second-action support statistics;
- maximum spin reached.

## 8. Pass/fail versus scientific outcome

Implementation PASS requires:

- all five target-node shards are present;
- first-layer Grams agree;
- `G`, every `D^(w)` and `D` are Hermitian PSD within tolerance;
- generalized/whitened Ritz matrix is Hermitian PSD;
- all required columns were evaluated.

The following are **not** pass criteria:

```text
Ritz minimum decreases
Ritz minimum is zero
Ritz minimum is nonzero
particular generated-layer rank
particular dressing ratio
```

All are scientific outputs.

## 9. Interpretation tree

### If an exact robust zero appears

Then the first Euclidean representation-dressed habitat contains a common normal-constraint zero direction. It must next be checked against

- tighter amplitude thresholds;
- tangential/diffeomorphism constraints;
- Lorentzian constraint dressing;
- projector/refinement stability.

It is not yet a physical scalar or dark mode.

### If the minimum decreases but remains positive

This is evidence that representation dressing moves the boundary carrier toward the Euclidean constraint surface. It motivates K2/block-Lanczos and the near-zero refinement analysis, but does not establish a zero sector.

### If the minimum does not decrease

The first dressing layer is insufficient. The result is retained and the next layer is justified only if independent physical/constraint structure requires it; one may not add depth merely until a desired zero appears.

## 10. Full-physics boundary

This is a **normal Euclidean constraint** Ritz calculation. `FULL_CONSTRAINT_FAMILY_PHYSICAL_MASTER.md` still requires the regulator-safe tangential/diffeomorphism sector or a theorem that it vanishes on the same asymptotic low sector. Lorentzian `H_L`, physical history, source dressing, scalar cosmology and dark-sector interpretation remain separate gates.
