# q=2 -> exact causal-volume dimension-three fixed point

Status: **exact asymptotic scaling theorem for the frozen q=2 route rewrite, cross-checked independently by PL topology and diffusion/scaling observables**.

This note closes the specific “dimension gradient” that had previously appeared only numerically as

```text
... -> 2.9517 -> 2.99385 -> 2.999229782 -> ...
```

The limit is not guessed from a fit. Once `q=2` has been selected independently by the local binary homogeneity equation

\[
q+2=2^q,
\]

the frozen route rewrite has an analytic growth law whose unique infrared causal-volume exponent is exactly three.

---

## 1. Independent microscopic input: q=2

The local selector is

\[
q+2=2^q.
\]

For integers `q>=1` its unique solution is

\[
\boxed{q=2}.
\]

Therefore the number of route midpoints created per active causal edge is

\[
B=2^q=4.
\]

This value is fixed before any dimension measurement below.

---

## 2. Exact recursive growth law

![Комбинаторная размерность и динамический масштаб](docs/figures/04-emergent-dimensions.svg)


In one generation every active causal edge produces `B` intermediate route vertices and two causal child edges per route. Therefore

\[
E_{g+1}^{\rm active}=2B\,E_g^{\rm active}.
\]

For `q=2`,

\[
\boxed{2B=8}.
\]

The causal graph distance between the original boundary vertices doubles per generation,

\[
L_{g+1}=2L_g.
\]

Starting from one active edge with two endpoint vertices, the exact vertex count is

\[
N_g
=2+B\sum_{r=0}^{g-1}(2B)^r
=2+B\frac{(2B)^g-1}{2B-1}.
\]

For `B=4`,

\[
\boxed{N_g=\frac{4\,8^g+10}{7}}.
\]

No continuum fit is used in this identity.

---

## 3. The finite-step dimension gradient

Define the one-step causal-volume exponent by comparing consecutive generations,

\[
d_g
=\frac{\log(N_g/N_{g-1})}{\log(L_g/L_{g-1})}
=\log_2\frac{N_g}{N_{g-1}}.
\]

Substituting the exact `q=2` counts gives

\[
\boxed{
d_g
=3+\log_2\left(
1-\frac{35}{16\,8^{g-1}+40}
\right).
}
\]

Hence for every finite `g>=2`,

\[
d_g<3,
\]

and because the positive correction in the logarithm decreases monotonically,

\[
d_{g+1}>d_g.
\]

Therefore

\[
\boxed{d_g\nearrow3}.
\]

Representative values are

| generation `g` | `d_g` |
|---:|---:|
| 2 | 2.662965012722429 |
| 3 | 2.951744831392779 |
| 4 | 2.993853015664851 |
| 5 | 2.999229782139151 |
| 6 | 2.999903693848493 |
| 7 | 2.999987961279020 |
| 8 | 2.999998495152814 |

The previously frozen held-out causal-volume value

\[
\boxed{d_H=2.999229782139151}
\]

is therefore not an isolated numerical coincidence: it is the `g=5` point of this exact monotone sequence.

---

## 4. The fixed point is exactly three

As `g -> infinity`, additive boundary terms in `N_g` become irrelevant and

\[
\frac{N_g}{N_{g-1}}\to 2B=8.
\]

Consequently

\[
\boxed{
d_*^{\rm causal-volume}
=\log_2(8)
=3.
}
\]

Equivalently, for general frozen `q`, the same route rule has

\[
d_*^{\rm causal-volume}=\log_2(2^{q+1})=q+1.
\]

Because `q=2` was selected independently, the dimension is not inserted by choosing a three-dimensional lattice.

---

## 5. Why this is stronger than one fitted exponent

The project now has three logically different reasons that point to the same spatial dimension.

### A. Exact local topology

The q=2 route shell is an octahedral

\[
S^2.
\]

A vertex of a combinatorial three-manifold has an `S^2` link.

### B. Exact global PL completion and stability

The canonical minimal+flag completion is the boundary of the 16-cell,

\[
M^3\cong S^3,
\]

with

```text
(V,E,F,T) = (8,24,32,16)
Betti     = (1,0,0,1)
```

and recursive barycentric refinement preserves all vertex/edge/face link conditions through the tested levels

```text
16 -> 384 -> 9216 tetrahedra.
```

### C. Dynamical/scaling observables

The frozen route calculation separately gives

\[
\boxed{d_H=2.999229782139151},
\qquad
\boxed{z=0.998281156\ \text{(rounded project value)}}.
\]

In the repository code, the quantity historically labelled `ds_slice_holdout` is already the ratio

\[
\boxed{
d_{\rm eff}^{\rm slice}
\equiv\frac{d_H}{z}
\simeq3.004393867.
}
\]

This notation matters: `3.004393867` has **already been divided by `z` once**. It must not be divided by `z` again.

Together with the exact causal-volume limit,

\[
\boxed{
D_{\rm topo}=3,
\qquad
d_{\rm causal-volume}\to3,
\qquad
d_{\rm eff}^{\rm slice}\simeq3.00439,
\qquad z\simeq1.
}
\]

These are not the same estimator repeated three times.

---

## 6. From three spatial dimensions to 3+1 scaling

The frozen route code defines the one-causal-time history exponent as

\[
\boxed{
d_{\rm eff}^{\rm history}
=1+\frac{d_H}{z}
=1+d_{\rm eff}^{\rm slice}.
}
\]

With the frozen held-out values,

\[
\boxed{
d_{\rm eff}^{\rm history}
\simeq4.004393867.
}
\]

This is the internally consistent project convention.

A separate continuum heat-kernel formula is sometimes written

\[
d_s^{\rm history}=1+\frac{D}{z},
\]

where `D` denotes the spatial volume/Hausdorff dimension. If one calls `D/z` itself a spatial spectral exponent, then the equivalent formula is simply

\[
d_s^{\rm history}=1+d_s^{\rm slice}.
\]

The two forms are the same statement with different notation; applying `/z` to `d_H/z` a second time would double-count the dynamical exponent.

Thus the strengthened dimension chain is

\[
\boxed{
q=2
\to S^2\ \text{local link}
\to M^3\cong S^3\ \text{PL phase}
\to d_{\rm causal-volume}=3\ \text{fixed point}
\to d_H/z\simeq3.00439
\to z\simeq1
\to 3+1\ \text{IR scaling}.
}
\]

---

## 7. Exact scope and non-claims

The theorem proved here is the asymptotic **causal-volume exponent of the frozen q=2 route rewrite** relative to its causal depth scale.

It does not by itself prove that every graph metric, every possible nonflag global gluing, or every interacting quantum ensemble has Hausdorff dimension three. That is why the repository keeps the PL-link, diffusion/scaling and Hodge/two-form gates independent.

Likewise, the existence/stability of the canonical `S^3` PL completion is not a theorem that the bare causal graph uniquely forces that gluing.

The scientific strength comes from agreement of independent observables, not from relabelling one exponent as several different results.

---

## 8. Reproduction

```bash
python scripts/q2_dimension3_fixed_point_gate.py \
  --max-generation 10 \
  --output verification_results/Q2_DIMENSION3_FIXED_POINT.json
```

The executable checks the closed-form identity, monotonic convergence from below and the exact fixed point `d*=3`.
