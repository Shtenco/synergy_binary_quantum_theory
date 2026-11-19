# Signed operator-first full-HDA certificate

Status: **conditional asymptotic closure theorem + finite spin-changed route evidence**.

This certificate upgrades the older fixed-cutoff composition argument to the current production definitions:

\[
H_E^{sine}=(T-T^\dagger)/(2i),
\]

\[
G_v=-\frac23E_v+\frac{32i}{9}L_{raw,v}
\qquad(\beta=\hbar=1),
\]

and

\[
R_{op}[N]=\frac12\left\{N,\sqrt{\hat Q^{ab}\hat P_a\hat P_b}\right\}.
\]

No Euclidean or Lorentzian coefficient is fitted to the HDA result.

## 1. Geometry sector is bounded at fixed safe cutoff

The exact hit-depth audit places a complete Lorentzian HH pair inside the declared wall

\[
J_{max}=13/2
\]

for the frozen all-`j=1/2` seed. At a fixed finite Peter-Weyl cutoff, the local signed geometry operators `G_v` therefore act on finite-dimensional regulator-safe sectors and have finite norm.

The exact signed channel expansion is already frozen:

\[
\boxed{
[G_0,G_1]
=\frac49 EE-\frac{64i}{27}(EL+LE)-\frac{1024}{81}LL.
}
\]

The physical-sine `EE` channel has independently completed:

```text
support = 514
raw norm = 2.879453814704955.
```

The brute-force `EL/LE/LL` workers reached their six-hour runner wall and were cancelled before writing artifacts. This is a computational timeout, not a physics failure; those channel-resolved calculations remain an independent finite falsifier.

## 2. Operator-first route is valid after genuine geometry change

The route operator is not expectation-valued. It uses the positive spectral square root of the complete fixed-spin geometry block.

A dedicated gate applied `R_op` to **five distinct genuine higher-spin sectors reached by `H_E^sine`**. It passed every positivity and HDA-scaling condition.

Initial all-`j=1/2` sector:

```text
Delta_R(1/64) = 8.205159710207802e-7
p_R           = 0.9999594708960342
Q error vs independent 4x4 construction = 0.
```

Five spin-changed sectors give

```text
endpoint defects:
1.9931725228e-6
2.4783541485e-6
2.2654862408e-6
3.6365781785e-6
9.3706514344e-7

fitted epsilon exponents:
0.9998813243
0.9999820816
0.9999698910
0.9999362611
0.9999583732.
```

The minimum checked spectral eigenvalue is `-1.07e-14`, i.e. zero within floating-point error; all positivity checks pass.

Provenance:

```text
workflow run: 31858615323
head:         905358f76699370aa13e017bd852c45b696b3e5f
artifact:     9244277324
digest:       sha256:c1af8de00183fddf328f6bdfba386e2320b842e10d3de98d90ad150b0876213c
```

Thus the production route operator used *after* a spin-changing geometry action is now independently tested, not inferred from the old expectation-first surrogate.

## 3. Exact cancellation of the dangerous mixed `1/epsilon` term

Use the frozen habitat family

\[
N=\bar N+\epsilon n,
\qquad
M=\bar M+\epsilon m,
\]

and the WKB scaling

\[
\Omega_Q=\epsilon^{-1}\widetilde\Omega_Q.
\]

For a geometry transition `Q0 -> Qg`, write

\[
\Delta R_M
=\frac{\bar M}{\epsilon}\Delta\widetilde\Omega+\Delta S_m,
\]

\[
\Delta R_N
=\frac{\bar N}{\epsilon}\Delta\widetilde\Omega+\Delta S_n.
\]

At node `v`,

\[
N_v=\bar N+\epsilon n_v,
\qquad
M_v=\bar M+\epsilon m_v.
\]

Then **before taking any matrix element**,

\[
\begin{aligned}
N_v\Delta R_M-M_v\Delta R_N
={}&\bar M n_v\Delta\widetilde\Omega
-\bar N m_v\Delta\widetilde\Omega\\
&+\bar N\Delta S_m-\bar M\Delta S_n\\
&+\epsilon(n_v\Delta S_m-m_v\Delta S_n).
\end{aligned}
\]

The apparent `1/epsilon` contribution cancels algebraically. Therefore a bounded local `G_v` gives an absolute mixed cross of order `O(1)`.

The route diffeomorphism target has norm

\[
\|D\|=O(\epsilon^{-1}),
\]

so

\[
\boxed{C_{G\times R}/D=O(\epsilon).}
\]

This argument is operator-first: `Delta Omega` is a finite matrix difference between positive spectral square-root blocks. No expectation value of `Q` is inserted before the square root.

## 4. Pure geometry channel

For the two node lapse coefficients,

\[
\begin{aligned}
N_0M_1-N_1M_0
={}&\epsilon\big[\bar N(m_1-m_0)+\bar M(n_0-n_1)\big]\\
&+\epsilon^2(n_0m_1-n_1m_0).
\end{aligned}
\]

There is no zeroth-order term. At fixed safe cutoff `[G0,G1]` is bounded, hence the absolute pure-geometry contamination is `O(epsilon)`. Relative to `D=O(epsilon^-1)`,

\[
\boxed{C_{GG}/D=O(\epsilon^2).}
\]

This conclusion does not depend on the eventual finite norm of the `EL/LE/LL` sum; those numbers determine finite calibration, not the asymptotic power.

## 5. Full fixed-cutoff HDA conclusion

Combining the tested operator-first route residual with the two exact composition identities gives

\[
\boxed{
\Delta_{full}
\le
\Delta_{R,op}
+C_1\epsilon+C_2\epsilon^2
\longrightarrow0
}
\]

at every fixed regulator-safe finite cutoff on the declared two-node WKB habitat.

This is the current **full signed operator-first HDA candidate theorem**.

It is stronger than the older architecture certificate in two respects:

1. the Euclidean ordering is the physically used `H_E^sine`, with a preregistered finite HDA PASS;
2. the route square root remains operator-first on genuine spin-changed sectors and has directly measured near-unit regulator exponents there.

## 6. Conditional simultaneous cutoff

The separately frozen norm envelope is

\[
C_{G\times R}/D
=O(\epsilon J_{max}^{13/2}),
\]

\[
C_{GG}/D
=O(\epsilon^2J_{max}^{13}).
\]

For

\[
J_{max}\sim\epsilon^{-\alpha},
\]

both vanish when

\[
0<\alpha<2/13.
\]

BCQG v1 uses the explicit interior path

\[
\boxed{\alpha=1/8},
\]

which gives

\[
\boxed{C_{G\times R}/D=O(\epsilon^{3/16})},
\]

\[
\boxed{C_{GG}/D=O(\epsilon^{3/8})}.
\]

This is a conditional diagonal continuum path, **not** a uniform arbitrary-path theorem.

## 7. What the still-running finite calculation can falsify

The exact off-shell finite calculation remains scientifically valuable. It must assemble

\[
[H_E^{sine}+(1+\beta^2)H_L+R_{op},
 H_E^{sine}+(1+\beta^2)H_L+R_{op}]
\]

with frozen signed coefficients and compare the five `epsilon` points against the exact diffeomorphism target.

A finite FAIL must be retained and can identify:

- a hidden scalar-projection/order problem in `H_L`;
- an unexpectedly large finite `EL/LE/LL` coefficient;
- an error in the full `G x R_op` implementation;
- or a habitat-specific obstruction not visible in the asymptotic norm theorem.

But an Actions timeout is not such a FAIL.

## Certificate statement

> **Conditional on the declared finite-cutoff domains and on the frozen polynomial joint-cutoff norm envelope, BCQG now has an operator-first, signed, physically sine-ordered full-HDA asymptotic closure theorem. The remaining heavy `EL/LE/LL` computation is a channel-resolved finite falsifier/calibration rather than the only mathematical bridge to closure.**

Executable symbolic gate: `scripts/full_operator_first_hda_theorem_gate.py`.

Spin-changed route evidence: `verification_results/PETER_WEYL_OPERATOR_ROUTE_SPINCHANGED_BLOCKS.json`.
