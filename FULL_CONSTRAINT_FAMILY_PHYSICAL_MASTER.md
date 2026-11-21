# Full constraint family required by the BQG physical master

Status: **physicalization scope correction, refined by the exact HDA-kernel theorem.** The current five-node Euclidean master remains a normal-constraint diagnostic, not by itself the full refinement-level physical projector.

## 1. Constraint hierarchy

The Peter-Weyl spin-network basis used by the current K5 calculations is already Gauss reduced at every four-valent node. Thus the local SU(2) Gauss constraint is implemented kinematically in the declared basis.

The remaining gravitational structure contains

- node normal/Hamiltonian constraints `H_v`;
- the commutator-derived HDA action
  \[
  D^{comm}_{vw}:=-i[H_v,H_w];
  \]
- and, whenever the finite regulator defines one independently, a target tangential/diffeomorphism operator `D_I^target` against which the HDA commutator is compared.

These last two objects must not be conflated.

## 2. Exact theorem: a commutator block is redundant on an exact common Hamiltonian kernel

If

\[
H_v|\psi\rangle=0\qquad\forall v,
\]

then identically

\[
[H_v,H_w]|\psi\rangle=0\qquad\forall v,w.
\]

Therefore

\[
\boxed{
\ker\sum_v H_v^\dagger H_v
=
\ker\left[
\sum_v H_v^\dagger H_v
+
\sum_{v<w}(D^{comm}_{vw})^\dagger D^{comm}_{vw}
\right].
}
\]

So adding the commutators of the **same** `H_v` as a second master block cannot further select an exact common zero state. `scripts/hda_kernel_redundancy_gate.py` verifies this algebraically/numerically on a noncommuting finite control.

This is an exact operator statement and does not require the continuum HDA to be anomaly-free.

## 3. Why finite HDA defects still matter

The BQG HDA test is stronger than the identity above because it compares the quantum commutator with a separately defined target tangential action:

\[
D^{comm}_{vw}
\stackrel{IR}{\longrightarrow}
D^{target}[\beta_{vw}(q,N,M)].
\]

At finite regulator the measured defect

\[
\Delta^{HDA}_{vw}
=D^{comm}_{vw}-D^{target}_{vw}
\]

need not vanish exactly. Hence

\[
D^{comm}_{vw}|\psi\rangle=0
\]

does **not** imply

\[
D^{target}_{vw}|\psi\rangle=0
\]

unless the defect itself vanishes on the candidate physical sector.

Thus there are two legal routes:

### Route A — explicit target tangential master

If regulator-safe independent `D_I^target` matrices are available on the same habitat, use

\[
\boxed{
\mathbb M_{full}
=
\sum_{v,w}H_v^\dagger G_H^{vw}H_w
+
\sum_{I,J}(D_I^{target})^\dagger G_D^{IJ}D_J^{target}.
}
\]

### Route B — Hamiltonian exact/low projector plus HDA-target certificate

Use the Hamiltonian master

\[
\mathbb M_H=\sum_v H_v^\dagger H_v
\]

and separately require, for every target tangential generator,

\[
\boxed{
\|D_I^{target}P_{phys}^H\|=0
}
\]

for an exact finite zero sector, or

\[
\boxed{
\|D_I^{target}P_{low}^H\|\to0
}
\]

under refinement for a near-zero/rigging sector.

Only after that certificate may `P_H` be promoted to the full physical projector.

## 4. Near-zero warning is genuinely stronger

For approximate states it is not enough that

\[
\langle\psi|\mathbb M_H|\psi\rangle\to0.
\]

`scripts/hda_kernel_redundancy_gate.py` contains a Hermitian counterexample with

\[
\langle M_H\rangle=2\epsilon^2\to0,
\]

while

\[
\|[H_1,H_2]\psi\|=1
\]

because a constraint operator norm grows as `1/epsilon`.

Therefore even the commutator residual itself must be controlled in a near-zero/refinement construction unless a uniform operator bound proves the implication. The independent `D^target` residual remains mandatory as well whenever the regulator defines that target separately.

## 5. Lorentzian normal master

On one common parity-complete habitat define

\[
\boxed{
H_v(\lambda_L)=H_v^E+\lambda_L H_v^L.
}
\]

The normal master is the matrix pencil

\[
\boxed{
M_H(\lambda_L)
=M_{EE}+\lambda_L M_{EL}+\lambda_L^2M_{LL},
}
\]

with

\[
M_{EE,ij}=\sum_v\langle E_i^{(v)}|E_j^{(v)}\rangle,
\]

\[
M_{LL,ij}=\sum_v\langle L_i^{(v)}|L_j^{(v)}\rangle,
\]

\[
M_{EL,ij}=\sum_v\left(
\langle E_i^{(v)}|L_j^{(v)}\rangle+
\langle L_i^{(v)}|E_j^{(v)}\rangle
\right).
\]

On the pure even q=2 boundary the mixed block is forbidden by doubled-spin parity, but on an enlarged even+odd candidate habitat it can be nonzero and must be retained.

The relative Lorentzian coefficient/order convention is theory input and may not be fitted from DM, DE, lensing or gravitational-wave data.

## 6. Production data rule: preserve outgoing columns

The physical master depends on complete outgoing images, not on direct logical returns. Every expensive calculation must therefore preserve

\[
|E_i^{(v)}\rangle=H_v^E|b_i\rangle,
\qquad
|L_i^{(v)}\rangle=H_v^L|b_i\rangle
\]

in a common sparse Gauss/Peter-Weyl basis before any projection.

The Lorentzian aggregate now serializes its complete Gauss outgoing column so it can be reused directly in Gram/master assembly. Future Euclidean/Lorentzian column jobs must follow the same rule. Recomputing a column merely because a downstream Gram is needed is prohibited workflow design.

## 7. Suggested block-Krylov organization

Starting from the q=2 boundary block `B`, retain separately labelled generated vectors:

\[
\mathcal K_1=\operatorname{span}\{B,H_vB\},
\]

\[
\mathcal K_2\supset\operatorname{span}\{B,H_vB,H_wH_vB\}.
\]

The same two-H data used for HDA closure can measure both target/commutator residuals and Ritz master spectra. Constraint labels must remain explicit: summing node Hamiltonians before squaring can create accidental cancellations that are absent from

\[
\sum_v\|H_v\psi\|^2.
\]

## 8. Physical-projector decision tree

At each regulator/refinement level:

1. construct complete separately labelled `H_v^E` and `H_v^L` images on one habitat;
2. freeze the Lorentzian coefficient/order convention;
3. assemble `M_H(lambda_L)` from complete outgoing columns;
4. inspect exact/near-zero spectral separation using `PHYSICAL_PROJECTOR_NEAR_ZERO_RIGGING_LIMIT.md`;
5. verify commutator/HDA residuals and the independent `D^target` residual when such a target is declared;
6. compare low-subspace projectors under refinement embeddings;
7. compute the q=2 boundary overlap;
8. feed the resulting `P_BQG` into the already existing relational/source stack;
9. only then evaluate scalar, TT and FLRW sectors of the same `Gamma_BQG`.

## 9. Consequence for scalar cosmology

A low mode of an unreduced normal master can still be

- a target-tangential/gauge mode because of finite HDA defect;
- a regulator anomaly;
- or a genuine physical scalar response.

Only the last survives the target-HDA/refinement checks and then the scalar source/gauge reduction. No normal-master eigenvalue by itself is dark matter or dark energy.

## 10. Claim boundary

The five-node Euclidean boundary master and the K1 outgoing-span/Ritz calculations remain valuable because they measure the true normal-constraint geometry and the minimal generated habitat. The full theory-specific projector is obtained only after the Lorentzian constraint family and the declared HDA-target/refinement conditions are incorporated. No new projector, Hilbert-space formalism, source functional or observable dictionary is required downstream; those components already exist and are reused.
