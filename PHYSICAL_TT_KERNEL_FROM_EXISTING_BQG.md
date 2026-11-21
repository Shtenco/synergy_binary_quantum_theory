# Conditional ADM/HDA leading TT IR form from existing BQG structure

Status: **the registered local two-derivative ADM/HDA premises imply the standard leading TT form; the theory-specific physical 1PI BQG kernel remains open.**

This file introduces no new microscopic degree of freedom, clock, history operator or transfer coefficient. It records what follows conditionally from the already registered local ADM/HDA sector and, just as importantly, what does **not** follow until the physical projector/history and connected generating functional are constructed.

## 1. Registered conditional input

Inside the declared local two-derivative ADM family,

\[
H[N]=\int d^3x\,N\left[
A\frac{\pi_{ab}\pi^{ab}-c\pi^2}{\sqrt q}
-B\sqrt q(R-2\Lambda)
\right]+O(\partial^4),
\]

the existing HDA result selects

\[
\boxed{c=\frac12,\qquad AB=1.}
\]

This is a structural/IR statement under the registered hypotheses. It is not by itself a derivation of the BQG physical history measure or physical inner product.

## 2. TT reduction of that ADM family

For two orthonormal TT polarizations `h_A`, `A=1,2`,

\[
H_{TT}^{(2)}
=\sum_A\int d^3x\left[
A\pi_A^2+\frac B4(\partial_i h_A)(\partial_i h_A)
\right]+O(\partial^4).
\]

Since

\[
\dot h_A=2A\pi_A,
\qquad
\pi_A=\frac{\dot h_A}{2A},
\]

the ordinary canonical Legendre transform of this selected ADM model gives

\[
\boxed{
S_{TT,ADM}^{(2)}
=\frac1{4A}\sum_A\int d^4x
\left[\dot h_A^2-AB(\nabla h_A)^2\right].
}
\]

Using `AB=1`,

\[
\boxed{
S_{TT,ADM}^{(2)}
=\frac1{4A}\sum_A\int d^4x
\left[\dot h_A^2-(\nabla h_A)^2\right]
+O(\partial^4).
}
\]

Therefore the **conditional leading ADM/HDA TT reference form** is massless, two-polarization and relativistic:

\[
\boxed{m_{TT}=0,\qquad N_{TT}=2,\qquad c_T^2=1.}
\]

No constraint spectral parameter is used as physical frequency in this derivation.

## 3. Reference Hessian notation

If this local ADM action is written in Euclidean quadratic-Hessian notation,

\[
S_{E,TT}^{(2)}
=\frac12\int\frac{d\omega_Ed^3k}{(2\pi)^4}
\,h_A(-p)\,K^{ADM}_{E,AB}(p)\,h_B(p),
\]

then

\[
\boxed{
K^{ADM}_{E,TT}(\omega_E,\mathbf k)
=Z_T(\omega_E^2+\mathbf k^2)I_2+O(\partial^4),
\qquad Z_T=\frac1{2A}.
}
\]

With the conventional GR parametrization `A=16 pi G`, this becomes `Z_T=1/(32 pi G)` in the stated polarization normalization.

Likewise the conditional Lorentzian ADM/HDA reference form is

\[
\boxed{
K^{ADM}_{TT}(\omega,\mathbf k)
=Z_T[-(\omega+i0)^2+\mathbf k^2]I_2+O(\partial^4).
}
\]

The notation `K_ADM` is intentional: this result must not be silently promoted to the theory-specific physical BQG `Gamma_phys^(2)`.

## 4. Why the physical machine gate remains open

The repository's legal physicalization order remains

\[
\boxed{
\{C_A\}
\to M_G
\to P_{phys}/\eta
\to Z_{phys}[J_g]
\to W_{phys}[J_g]
\to\Gamma_{phys}[g]
\to\Gamma^{(2)}_{metric,phys}(\omega,\mathbf k)
\to\Gamma^{(2)}_{TT,phys}.
}
\]

The candidate theory has not yet completed the theory-specific refinement/rigging/boundary-history step and therefore has not yet derived the connected physical BQG 1PI metric kernel.

Thus the machine truth must remain

```text
physical_TT_kernel_frozen = false
```

until the final source-derived object exists.

## 5. Constraint resolvents remain constraint resolvents

The prohibition remains exact:

\[
\boxed{z_{constraint}\ne\omega_{physical}}
\]

unless an independent physical-time/history construction derives an identification.

The same applies to:

- master-constraint gaps;
- Lanczos shell eigenvalues;
- C8 character angles;
- reduced-control clock labels.

They may be microscopic inputs. They are not physical frequencies by renaming.

## 6. Relation to the five-block signed-G calculation

The microscopic spatial calculation now being built uses the separately frozen signed gravitational constraint

\[
\boxed{
G_{frozen}=-\frac23H_E^{sine}-\frac{32}{9}S,
\qquad
S=-\frac i2(L_{raw}-L_{raw}^\dagger).
}
\]

Its centered five-block chain is

```text
G_frozen
 -> complete common P+Q basis
 -> zero-energy Schur precursor
 -> K_q
 -> measured M_hq
 -> K_h(k)
 -> TT
 -> c_micro_spatial_BQG
```

This chain is useful because it produces a theory-specific microscopic spatial fingerprint without adding phenomenological coefficients. But it is still a **constraint-derived microscopic precursor**, not the physical connected 1PI kernel.

The extractor therefore deliberately outputs

```text
c_micro_spatial_BQG = <six numbers, if all guards pass>
c_BQG_IR            = null
```

until physicalization is completed.

## 7. Six-dimensional quartic shell

The algebraically closed parity-even tetrahedral spatial TT sector remains six-dimensional:

\[
\Gamma^{(4)}_{TT,spatial}
\sim a_*^2\sum_{r=1}^6c_r W_r(\mathbf k).
\]

The six-basis theorem and extractor are already complete. What remains open is the **status of the coefficients**:

1. the signed five-block constraint calculation can produce a microscopic spatial six-vector;
2. the theory-specific connected physical history must determine the physical pole six-vector;
3. only the latter may be frozen as `c_BQG_IR` and carried into a blind external comparison.

Independent `omega^4` and `omega^2 k^2` quartic structures are also outputs of the eventual physical effective action and are not fixed by the spatial six-vector alone.

## 8. Absolute normalization

HDA fixes the relative leading kinetic/gradient structure in the declared ADM family, not the microscopic physical value of the overall normalization `A` or equivalently `G`.

Therefore both remain distinct:

```text
conditional leading ADM/HDA TT shape: DERIVED WITHIN REGISTERED IR PREMISES
microscopic absolute Newton/residue scale: OPEN
physical BQG connected TT 1PI kernel: OPEN
physical quartic TT six-vector: OPEN
```

## 9. Exact current statement

The strongest claim supported by this derivation is:

> Given the registered local two-derivative ADM/HDA premises, the TT reduction has the ordinary massless two-polarization relativistic leading form up to one overall normalization.

It is **not**:

> BQG has already derived its physical interacting graviton 1PI kernel from the microscopic quantum history.

That stronger statement remains a target of the physicalization bridge. No experimental confirmation is claimed.