# Preregistration: regulator-safe logical master Gram construction

Status: **frozen integration protocol before the full Lorentzian logical result is known**.

The direct-return question

\[
P_{log}H_LP_{log}=0\;?
\]

is scientifically useful but is not the physical-constraint question. A constraint may map a logical state entirely out of the logical sector while remaining nonzero. The next production-oriented finite object is therefore the **Gram/master matrix of complete outgoing states**.

---

## 1. Frozen logical input domain

Use the complete canonical all-`j=1/2` K5 Gauss basis

\[
\mathcal H_{log}=\operatorname{span}\{\psi_i\}_{i=1}^{32}.
\]

No pair partial trace is taken before the master matrix is assembled.

The first calculation uses one local normal-constraint node `v=0`; the all-node extension is a separate mandatory step.

---

## 2. Constraint images, not direct returns

For every input basis state compute the complete sparse Peter--Weyl images

\[
|E_i^{(v)}\rangle=H_{E,v}^{sine}|\psi_i\rangle,
\]

\[
|L_i^{(v)}\rangle=H_{L,v}^{raw}|\psi_i\rangle,
\]

using the same exact ordering, zero-aware volume convention and regulator walls frozen in the existing Euclidean and Lorentzian gates.

The Lorentzian image is the complete 24-term epsilon-oriented sum from `PETER_WEYL_LORENTZIAN_LOGICAL_RETURN_PREREGISTRATION.md`. It is **not** projected back to `H_log` before inner products are formed.

---

## 3. Frozen one-node Gram blocks

Define

\[
M_E^{(v)}{}_{ij}
=\langle E_i^{(v)}|E_j^{(v)}\rangle,
\]

\[
M_L^{(v)}{}_{ij}
=\langle L_i^{(v)}|L_j^{(v)}\rangle,
\]

and the mixed block

\[
X_{EL}^{(v)}{}_{ij}
=\langle E_i^{(v)}|L_j^{(v)}\rangle.
\]

Each diagonal Gram is positive semidefinite by construction.

The doubled-spin grading is already frozen:

```text
H_E : odd
H_L : even
```

on the even all-`j=1/2` input sector. Therefore the Euclidean and Lorentzian output vectors live in orthogonal parity sectors and the preregistered exact target is

\[
\boxed{X_{EL}^{(v)}=0}
\]

up to numerical sparse-basis tolerance.

This parity target is independent of whether `P_log H_L P_log` is zero or nonzero.

---

## 4. Overall Lorentzian coefficient does not select the kernel

For a local real coefficient `lambda != 0`, write the formal normal constraint

\[
G_v=H_{E,v}^{sine}+\lambda H_{L,v}.
\]

On the declared parity-separated logical input,

\[
G_v^\dagger G_v
\quad\longrightarrow\quad
M_E^{(v)}+|\lambda|^2M_L^{(v)}
\]

because the mixed Gram vanishes.

For positive `|lambda|^2`,

\[
\boxed{
\ker\left(M_E^{(v)}+|\lambda|^2M_L^{(v)}\right)
=\ker M_E^{(v)}\cap\ker M_L^{(v)}.
}
\]

Therefore the **zero-sector dimension and vectors** do not require fitting `beta` or an overall Lorentzian normalization. Spectral gaps away from zero do depend on the relative coefficient and are not interpreted physically until the quantum normalization convention is frozen.

---

## 5. Hard numerical diagnostics

For each Gram matrix report:

- Hermiticity defect;
- minimum eigenvalue;
- rank/nullity under a declared relative tolerance;
- smallest positive eigenvalue;
- condition number on support;
- trace and Frobenius norm.

For the mixed block report

\[
\frac{\|X_{EL}\|_F}
{\sqrt{\|M_E\|_F\|M_L\|_F}}
\]

and require parity suppression below the frozen numerical tolerance appropriate to the existing sparse engine.

No eigenvalue is to be interpreted as a particle mass or physical frequency.

---

## 6. All-node master is mandatory

The actual K5 finite normal-constraint master uses the five local constraints,

\[
\boxed{
M_{normal}
=\sum_{v=0}^{4}
\left(M_E^{(v)}+M_L^{(v)}\right)
}
\]

for unit positive constraint-space metric when only the kernel is being tested.

More generally any strictly positive constraint-space metric `G^{AB}` has the same common-zero condition if the complete declared constraint family is used.

One-node nullity is therefore **not** the finite physical nullity.

A graph-relabeling/S5 reconstruction of the five node matrices is allowed only after its unitary basis map is derived and checked against at least one directly calculated nontrivial column at another node. Otherwise all five nodes must be computed directly.

---

## 7. Relation to the existing 32D Euclidean master code

`scripts/peter_weyl_master_32_gate.py` already implements the correct core technology:

```text
full 32 logical inputs
 -> complete outgoing sparse states
 -> full 32x32 Gram A^dagger A
 -> eigen/rank audit
 -> only afterwards any partial trace
```

The new production path reuses that ordering but replaces the historical Euclidean pair operator by the declared node-local Euclidean/Lorentzian constraint images.

---

## 8. Restricted-domain claim boundary

Even an exact all-five-node kernel inside `H_log` establishes only

\[
\ker M_{normal}\cap\mathcal H_{log}.
\]

It does **not** prove that the full physical Hilbert space is exhausted by `H_log`. A true finite physical projector on a graph-changing Peter--Weyl habitat must allow the full declared finite domain/codomain needed by the constraints and their zero sector.

Consequently:

- nonzero restricted nullity gives genuine finite candidate physical boundary states inside the microscopic logical carrier;
- zero restricted nullity means no exact normal-constraint state exists **inside that carrier**, not that the full theory has no physical states.

---

## 9. Source dressing after the master

If a nontrivial finite zero sector is obtained, the next calculation is not a resolvent fit. It is a source-deformed physical-boundary amplitude, schematically

\[
Z_{phys}[J]
=\langle\Psi_{out}|
P_{phys}\,\mathcal T\exp(J\cdot O_{rel})
|\Psi_{in}\rangle,
\]

or an explicitly declared finite heat-kernel/rigging approximation with the projector limit controlled.

The scalar sources must then include the derived collective volume/scale carrier and transported metric-shape directions. Connected interblock Hessians are taken from

\[
W=\log Z_{phys},
\]

not from the constraint resolvent.

Only after this step can the programme legitimately approach a physical `Gamma_scalar^(2)`.