# End-to-end control: face qubits -> Einstein

Status: **single-data-path finite Euclidean control; not microscopic geometrogenesis**.

## Purpose

Previous checks validated separate arrows.  This test forces the first reconstruction arrows to operate on one and the same dataset.  The reconstruction code receives six qubit density matrices per sample point, one on each independent oriented coordinate 2-face.  No metric, tetrad or connection is passed downstream.

For the positive control, the qubits encode the self-dual two-form field of the unit Euclidean four-sphere in stereographic coordinates.  After encoding, the original geometric field is discarded and reconstructed only through Pauli expectation values.

The chain is

\[
\boxed{
\rho_{\mu\nu}
\to B^i
\to \Delta_{simp}
\to g_U
\to A_B:\ D_AB=0
\to F(A_B)
\to (\bar F,\operatorname{Tr}F).
}
\]

The generic vacuum Einstein condition in this language is **not** that the tracefree self-dual block vanish.  Rather,

\[
\bar F^{ij}=0,
\qquad
\operatorname{Tr}F=-\Lambda,
\]

while the tracefree part of $F^{ij}$ is self-dual Weyl curvature.  The present positive control is $S^4$, whose Weyl tensor happens to vanish.

## Qubit encoding

For each independent face $(\mu\nu)$,

\[
\rho_{\mu\nu}
=\frac12\left[I+rac{B^i_{\mu\nu}}{4}\sigma_i\right].
\]

The factor 4 keeps every Bloch vector inside the unit ball for the unit-$S^4$ stereographic patch used here.  Reconstruction uses only

\[
B^i_{\mu\nu}=4\,\operatorname{Tr}(\rho_{\mu\nu}\sigma^i).
\]

This is an oracle-encoded control of composability, not a claim that a frozen rewrite rule has generated the qubits.

## Positive control: unit $S^4$

Across five independent sample points the single-path result is

| defect / observable | worst or reconstructed value |
|:--|--:|
| qubit -> $B$ decoding error | `0.0` |
| simplicity defect | `2.1274680382e-16` |
| Urbantke conformal metric error | `3.3306690739e-16` |
| $D_AB$ residual | `6.1920200345e-17` |
| anti-self-dual curvature defect | `8.6321136850e-09` |
| reconstructed $\Lambda=-\langle\operatorname{Tr}F\rangle$ | `2.999999897308107` |
| relative error vs exact unit-$S^4$ $\Lambda=3$ | `3.4230631070e-08` |
| spatial constancy defect of $\operatorname{Tr}F$ | `6.3528924421e-09` |

Thus the same qubit data reproduce

\[
\boxed{\Lambda_{rec}=2.9999998973}
\]

without passing the value of $\Lambda$ to the reconstruction chain.

## Negative control: metric-compatible but non-Einstein

The second dataset is also encoded only as face qubits.  It comes from a smooth conformally-flat, simple and nondegenerate $B$ field with

\[
\Omega(x)=\exp\left[0.15(x_0x_1+0.30x_2^2-0.20x_3)\right].
\]

It passes the early arrows:

- qubit decoding: machine precision;
- simplicity: machine precision;
- Urbantke metric: machine precision;
- compatible connection: machine precision.

But it fails exactly at the Einstein-curvature arrow:

\[
\boxed{
\Delta_{ASD}=0.73723\ldots0.74174.
}
\]

The minimum over the five points is

\[
\boxed{\Delta_{ASD}^{nonE}=0.737225989.}
\]

Therefore the end-to-end chain is not merely recognizing a smooth metric or a simple $B$ field.

## What is now mathematically connected

![Маршрут к физическому TT и IR action](docs/figures/10-continuum-tt.svg)


The finite control establishes composability of

\[
\boxed{
\text{face qubit data}
\to B
\to \text{simple metric sector}
\to \text{compatible connection}
\to \text{Einstein curvature test}.
}
\]

This is stronger than juxtaposing separate tests because all downstream quantities are reconstructed from the same qubit inputs.

## What remains open

The first physical arrow is still absent:

\[
\boxed{
\text{frozen local causal/frame rewrite dynamics}
\dashrightarrow
\{\rho_f,U_e,K_2\}\text{ in the required scaling phase}.
}
\]

The present test deliberately does not hide that gap.  It also does not establish:

- dimension emergence without the coordinate 2-face control scaffold;
- Lorentzian reality conditions;
- quantum measure/unitarity;
- two gapless physical graviton modes in the microscopic ensemble;
- matter/chirality/anomalies;
- experimental validity.

The Regge/Fierz--Pauli/EH/Ward calculations remain an **independent second route** from the reconstructed metric sector to continuum GR, rather than a logically necessary step in the Plebański route.

## Reproduction

```bash
python scripts/qubit_to_einstein_end_to_end.py \
  --B-step 2e-5 \
  --A-step 2e-4 \
  --output verification_results/qubit_to_einstein_end_to_end.json
```
