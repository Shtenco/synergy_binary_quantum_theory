# Preregistration: full five-node Euclidean Peter-Weyl master on the q=2 boundary carrier

Status: **frozen before the first five-node spectrum is evaluated.**

## 1. Question

The previously certified `32 x 32` logical Gram used the combined two-node image

\[
(H^E_0+H^E_1)|b_i\rangle.
\]

That calculation was already enough to show that the bare q=2 logical carrier is not closed under the Euclidean constraint. It is not, however, the canonical positive master for the complete five-node K5 constraint family.

The present gate asks the sharper finite question:

> Does the full set of five regulated Euclidean node constraints possess any common zero direction **inside the bare 32-dimensional all-`j=1/2` q=2 boundary carrier**?

No expected rank, nullity, eigenvalue or symmetry multiplicity is frozen in advance.

## 2. Frozen boundary basis

Use the existing ordered basis

\[
\mathcal B=\{|b_i\rangle\}_{i=0}^{31}
\]

returned by `scripts/k5_peter_weyl_safe_hda_column.py::basis_full_jhalf()`:

- all ten K5 links have `j=1/2`;
- each of the five four-valent nodes is in one of its two Gauss-singlet intertwiner states;
- hence `dim B = 2^5 = 32`.

No basis vector is selected or removed after seeing the result.

## 3. Frozen Euclidean operators and cutoff

For every K5 node

\[
v=0,1,2,3,4,
\]

use the existing orientation-covariant Hermitian Peter-Weyl Euclidean node operator implemented by

```text
PW.apply_H_cached_state(state, v, JMAX2)
```

in `scripts/k5_peter_weyl_safe_hda_column.py`.

Freeze

\[
\boxed{J_{max}=5/2.}
\]

This cutoff is not selected from the forthcoming spectrum. `K5_HH_REACHABLE_SPACE.md` proves that `Jmax=5/2` is safe for the present Euclidean `HH` family starting from the all-`j=1/2` K5 boundary; one `H` lies strictly inside the same wall and the two-H reachability actually saturates `j=5/2`.

The same amplitude pruning convention already used by the safe Peter-Weyl HDA code is retained. No new coefficient threshold is tuned from the master spectrum.

## 4. Master definition

For every boundary column and node compute the full outgoing sparse state

\[
|a_{v i}\rangle=H^E_v|b_i\rangle.
\]

The canonical identity-metric Euclidean master compressed to the boundary is

\[
\boxed{
M_E^{(B)}{}_{ij}
=
\sum_{v=0}^{4}
\langle a_{v i}|a_{v j}\rangle
=
B^\dagger\left(\sum_v H_v^{E\dagger}H_v^E\right)B.
}
\]

No cross-node terms

\[
H_v^\dagger H_w,\qquad v\ne w
\]

are included. This is deliberate: the master constraint with identity positive metric is a sum of positive node norms, not the Gram of the summed Hamiltonian.

For any vector `c` in the 32D carrier,

\[
c^\dagger M_E^{(B)}c
=
\sum_v\|H_v^E Bc\|^2\ge0.
\]

Therefore an exact zero eigenvector of this compressed master is equivalent, within the declared boundary carrier, to

\[
H_v^E Bc=0
\quad\forall v.
\]

## 5. Diagnostics frozen before the result

Report:

1. Hermiticity error of `M_E^(B)`;
2. minimum/maximum eigenvalues;
3. numerical rank and nullity at the declared relative tolerance;
4. smallest positive eigenvalue and condition number on support;
5. per-node Gram trace, Frobenius norm and image support statistics;
6. maximum direct return norm to the all-`j=1/2` boundary after one `H_v`;
7. pairwise node-trace spread as an S5/K5 permutation-covariance diagnostic;
8. the spectrum itself.

The direct return is expected to vanish from the already proved doubled-spin parity theorem, but it is rechecked independently here.

## 6. Pass/fail criteria

The scientific gate passes implementation consistency if:

- the master is Hermitian within numerical tolerance;
- it is positive semidefinite within numerical tolerance;
- the first-action direct return to the even all-`j=1/2` boundary is below the frozen amplitude tolerance;
- all five node families are evaluated;
- the five per-node traces agree to numerical covariance tolerance.

**Rank/nullity is explicitly not a pass criterion.**

Both outcomes are scientifically retained:

```text
nullity > 0
  -> there exists a common five-node Euclidean zero direction already inside
     the bare q=2 boundary carrier; characterize it before adding H_L.

nullity = 0
  -> the complete Euclidean constraint family has no physical zero vector
     inside the bare carrier; enlarged Peter-Weyl dressing is mandatory.
```

## 7. Claim boundary

Even if `nullity=0`, this does **not** prove that the enlarged Euclidean master has no zero sector. Conversely, if `nullity>0`, this does not close the Lorentzian physical projector.

The legal next object is the boundary overlap with the enlarged master projector,

\[
B^\dagger P_0 B,
\]

computed on a constraint-generated block-Krylov habitat and then extended from Euclidean to Euclidean+Lorentzian constraints.

No dark matter, dark energy, cosmological constant, physical frequency or propagating scalar is inferred from this finite master spectrum.
