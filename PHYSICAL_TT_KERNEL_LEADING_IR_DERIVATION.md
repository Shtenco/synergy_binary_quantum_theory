# Leading TT IR form from the registered ADM/HDA sector

Status: **conditional structural derivation; not the theory-specific physical BQG 1PI kernel.**

This note records one narrow implication of the already registered local two-derivative ADM/HDA sector. It adds no clock, microscopic operator, transfer coefficient or history postulate.

The question is:

\[
\boxed{\text{what leading TT form follows if the registered BQG ADM/HDA premises hold?}}
\]

The answer is the ordinary massless two-polarization relativistic **ADM/HDA reference form**, up to one overall normalization. The separate physicalization chain is still required to derive the actual connected BQG `Gamma_phys^(2)`.

## 1. Registered ADM/HDA premise

Inside the declared local family

\[
H[N]=\int d^3x\,N\left[
A\frac{\pi_{ab}\pi^{ab}-c\pi^2}{\sqrt q}
-B\sqrt q(R-2\Lambda)
\right],
\]

the existing structural HDA result gives

\[
\boxed{c=\frac12,\qquad AB=1.}
\]

The overall normalization remains free at this level.

## 2. TT reduction

For orthonormal TT polarizations `h_A`,

\[
H_{TT}^{(2)}
=\sum_A\int d^3x\left[
A\pi_A^2+\frac B4(\partial_i h_A)(\partial_i h_A)
\right].
\]

Hamilton's equation gives

\[
\pi_A=\frac{\dot h_A}{2A},
\]

and the canonical Legendre transform gives

\[
\boxed{
S_{TT,ADM}^{(2)}
=\frac1{4A}\sum_A\int d^4x
\left[\dot h_A^2-AB(\nabla h_A)^2\right].
}
\]

With `AB=1`,

\[
\boxed{
S_{TT,ADM}^{(2)}
=\frac1{4A}\sum_A\int d^4x
\left[\dot h_A^2-(\nabla h_A)^2\right].
}
\]

Therefore, within these premises,

\[
\boxed{m_{TT}=0,\qquad N_{TT}=2,\qquad c_T=1.}
\]

## 3. Conditional kernel notation

The corresponding Euclidean ADM/HDA quadratic kernel is

\[
\boxed{
K^{ADM}_{E,TT}(\omega_E,\mathbf k)
=Z_T(\omega_E^2+\mathbf k^2)I_2+O(\partial^4),
\qquad Z_T=\frac1{2A}.
}
\]

With the conventional GR parametrization `A=16 pi G`, `Z_T=1/(32 pi G)` in this polarization convention.

The Lorentzian reference form is

\[
\boxed{
K^{ADM}_{TT}(\omega,\mathbf k)
=Z_T[-(\omega+i0)^2+\mathbf k^2]I_2+O(\partial^4).
}
\]

`K_ADM` is used instead of claiming this object is already the theory-specific physical `Gamma_BQG^(2)`.

## 4. What is not derived here

This calculation does not provide:

- a theory-specific BQG rigging map or boundary-history amplitude;
- a physical inner product;
- `Z_phys[J_g]` or `W_phys[J_g]`;
- the connected physical metric 1PI effective action;
- the microscopic value of the overall Newton/residue normalization;
- physical quartic frequency structures;
- the physical six-Wilson vector.

Hence the repository machine flag remains

```text
physical_TT_kernel_frozen = false
```

until the physical source/history chain is actually completed.

## 5. Constraint spectral parameters are not omega

The rule remains

\[
\boxed{z_{constraint}\ne\omega_{physical}}
\]

unless an independent physical-history/time construction establishes the relation.

No master gap, Lanczos shell eigenvalue, C8 character angle or reduced-control frequency is used here as physical `omega`.

## 6. Relation to the microscopic signed-G bridge

The separate five-block calculation uses

\[
G_{frozen}=-\frac23H_E^{sine}-\frac{32}{9}S,
\qquad
S=-\frac i2(L_{raw}-L_{raw}^\dagger),
\]

and can produce

\[
\mathbf c_{micro,spatial}^{BQG}
\]

from the actual signed microscopic constraint after Schur reduction, metric transport and TT projection.

That quantity is a valuable theory-specific microscopic spatial precursor. It is deliberately **not** named `c_BQG_IR` until the physical connected generating functional has been derived.

## 7. Six-dimensional quartic spatial sector

The algebraic theorem remains:

\[
\boxed{\dim \mathcal W^{S4,parity-even}_{TT,k^4}=6.}
\]

Thus a completed physical spatial quartic kernel can be written

\[
\Gamma^{(4)}_{TT,spatial}
=Z_Ta_*^2\sum_{r=1}^6c_r^{IR}W_r(\mathbf k).
\]

The basis and exact six-observable extraction matrix are already closed. The physical values `c_r^IR` are not supplied by this ADM/HDA leading-order argument.

## 8. Exact status

```text
registered ADM/HDA -> conditional O(partial^2) TT shape: DERIVED
actual signed five-block microscopic spatial precursor: UPSTREAM HEAVY OPERATOR CALC OPEN
physical projector/history: OPEN
connected physical Gamma_TT^(2): OPEN
physical c_BQG_IR: OPEN
one common physical scale: OPEN
blind external comparison: NOT YET PERFORMED
```

No new microscopic entity is required by this statement, and no experimental confirmation is claimed.