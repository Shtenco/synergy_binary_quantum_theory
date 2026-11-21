# Repo-wide dependency graph for the first physical BQG scalar kernel

Status: **integration map; no new cosmological parameter is inserted here**.

The purpose of this document is to prevent the scalar/cosmology programme from duplicating machinery that already exists elsewhere in the repository, and to prevent constraint-space quantities from being silently renamed physical observables.

The target is

\[
\boxed{
\Gamma_{\rm scalar}^{(2)}(\omega,\mathbf k)
\longrightarrow
\rho_{\rm hist}(a),\Phi(a,k),\Psi(a,k),
\mu_{\rm BQG}(a,k),\Sigma_{\rm BQG}(a,k)
}
\]

with every function produced by one frozen BQG history/source construction rather than selected phenomenologically.

---

## 1. The existing repository already contains most of the kinematic/operator ladder

### A. Microscopic geometry and local sources

Existing exact pieces:

```text
q=2 route/Walsh tetrahedral geometry
 -> four-face Gauss singlet logical qubit
 -> logical X,Z shape observables
 -> exact shape-to-metric Jacobian
 -> relational C8/projector positive control
 -> exact local source Z(jX,jZ)
 -> exact local Gamma_shape
```

The newly merged scalar closure establishes

\[
\Gamma_{\rm shape}(s)
=s\operatorname{artanh}s+\frac12\log(1-s^2),
\]

but also the constructive no-go that the `X/Z` tangent is trace-free and lacks a linear FLRW volume mode.

### B. Collective volume/conformal carrier

Existing representation growth supplies

```text
j=1/2 shape carrier
 -> symmetric q=2 blocking
 -> j=1 four-valent intertwiner
 -> first non-scalar absolute-volume spectrum
```

with

\[
\operatorname{spec}V_{j=1}=\{0,3^{1/4},3^{1/4}\}.
\]

The same repository independently identifies common radial flux scaling as the DeWitt conformal direction. Therefore the missing scalar is not an arbitrary new field: the first structurally motivated volume/scale carrier appears at collective `j=1`.

### C. Real Lorentzian constraint amplitudes

Existing Peter--Weyl code already contains the sequence

\[
V
\to H_E^{\rm sine}
\to K=[V,H_E^{\rm sine}]
\to C_e(V),C_e(K)
\to \operatorname{Tr}_{aux}[C(K)C(K)C(V)].
\]

The ordered triple is a genuine sparse state-to-state amplitude at the preregistered `Jmax=7/2`, not support counting. Exact representation theory reduces scalar-relevant rank words to

\[
000,\;011,\;101,\;110,\;111.
\]

The branch `research/physical-scalar-kernel-dynamics` now preregisters the missing **24-term epsilon-oriented node sum** and its return to the all-`j=1/2` logical sector.

### D. Constraint algebra / lapse / diffeomorphism structure

Already available:

- exact dual-K5 lapse cochain `omega_vw=N_v M_w-N_w M_v`;
- regulator-safe finite three-node graph-changing HDA scaling;
- flux/habitat diffeomorphism target;
- inverse-volume-free densitized HDA target;
- DeWitt/HDA uniqueness and conformal-sign controls.

These are essential gauge/constraint data, but the lapse cochain is **not** itself the physical Newtonian potential `Psi`.

### E. Physical projector architecture

Already available:

\[
\mathbb M_G=\sum_{A,B}C_A^\dagger G^{AB}C_B\ge0,
\qquad
P_{phys}=\mathbf 1_{\{0\}}(\mathbb M_G),
\]

with

\[
\ker\mathbb M_G=\bigcap_A\ker C_A
\]

at finite regulator. This is the correct place to combine Euclidean, Lorentzian, Gauss/diffeomorphism-reduced normal constraints before source dressing.

### F. Source / effective-action architecture

The repository already has the exact positive-control arrow

\[
P_{rel}
\to O_{rel}
\to Z[J]
\to W[J]
\to \Gamma^{(2)},
\]

but with the deliberately nonphysical `R=J` C8 clock/system control. The production calculation must replace that control by the actual graph-changing physical projector/boundary amplitude.

### G. Interblock momentum architecture

The nearest shared-face programme has already reduced one reciprocal transfer to exactly six real `S3` amplitudes. `S4` then transports one canonical pair to the four tetrahedral neighbor directions.

The tetrahedral second moment is exactly

\[
\sum_a n_a^i n_a^j=\frac43\delta^{ij},
\]

so a derived reciprocal scalar nearest-neighbor transfer has an isotropic leading `k^2` symbol. The fourth moment retains the controlled cubic/tetrahedral memory used by the six-Wilson TT programme.

The missing object is not the momentum algebra. It is the **microscopic connected transfer amplitude with physical history weighting**.

---

## 2. The legal production chain

The shortest no-shortcut path is now

```text
1. finish full epsilon-oriented H_L amplitudes
2. build the finite declared constraint family C_A
3. form the finite master constraint M
4. isolate P_phys / controlled heat-kernel projector
5. freeze physical boundary/history normalization
6. define relational scalar source insertions
   - volume/scale
   - lapse-response probe
   - transported shape/shear
7. compute connected two-block source Hessians from W=log Z
8. transport the canonical shared-face result over S4 neighbors
9. form the low-k moment symbol
10. assemble the scalar ADM quadratic kernel
11. eliminate constraint/gauge variables by a Schur complement
12. transform the reduced metric response to Bardeen/Weyl observables
13. couple the same conserved matter source to dynamics and lensing
14. only then define mu_BQG, Sigma_BQG and search for physical scalar poles
15. evaluate Gamma on the homogeneous FLRW family to obtain rho_hist(a), p_hist(a)
```

No constraint resolvent parameter `z` enters this chain as physical frequency.

---

## 3. Minimal scalar quadratic carrier

A useful continuum notation for organizing the finite calculation is

\[
\chi=(c,s),
\]

where the constraint-like scalar variables are schematically

\[
c=(\delta N,B)
\]

and the spatial scalar metric variables are schematically

\[
s=(\zeta,E)
\]

or their discrete volume/shape counterparts.

The exact physical-history 1PI Hessian, once derived, can be block written as

\[
\Gamma^{(2)}_{\rm scalar}
=
\begin{pmatrix}
K_{cc}&K_{cs}\\
K_{sc}&K_{ss}
\end{pmatrix}.
\]

Before interpreting any pole, lapse/longitudinal-shift constraints must be solved or integrated out. When `K_cc` is invertible on the declared reduced domain, the reduced scalar kernel is the Schur complement

\[
\boxed{
K_{red}
=K_{ss}-K_{sc}K_{cc}^{-1}K_{cs}.
}
\]

If `K_cc` has gauge null directions, the correct constrained/pseudoinverse reduction must be used only after the corresponding gauge/Dirac structure is frozen.

A zero of an unreduced lapse/shift block is **not** evidence for a dark-matter particle.

---

## 4. Relation to the existing finite Dirac reduction

`K5_FINITE_DIRAC_REDUCTION.md` already demonstrates, in a regulator-unsafe `Jmax=1/2` control, that tangential-like commutators reduce the 32-dimensional fully active logical sector to a two-dimensional common scalar sector and the normal constraint then selects one branch.

This is valuable structural evidence for the ordering

```text
tangential/diffeomorphism reduction
 -> normal constraint
 -> physical scalar sector
```

but it must **not** be copied as the production scalar kernel because that finite test lies below the Peter--Weyl safe HDA wall. The production reduction must be repeated on the regulator-safe graph-changing habitat with the completed Lorentzian constraint.

---

## 5. How the physical source response is obtained

After reduction, couple one universally normalized conserved matter/probe source `T` through the same emergent metric dictionary used by gravity. Schematically,

\[
K_{red}\,s=-J_T.
\]

Then

\[
s=-K_{red}^{-1}J_T.
\]

Only after the discrete-to-continuum gauge-invariant map is frozen may one identify the resulting combinations with `Phi` and `Psi`.

The definitions to be matched are

\[
-k^2\Psi
=4\pi G a^2\mu_{\rm BQG}(a,k)\rho\Delta,
\]

\[
-k^2(\Phi+\Psi)
=8\pi G a^2\Sigma_{\rm BQG}(a,k)\rho\Delta.
\]

`mu` and `Sigma` therefore cannot be independently fitted. They are two projections of one reduced inverse kernel and one source coupling.

---

## 6. Modified gravity versus emergent dark degree of freedom

The completed calculation has two distinct possible outcomes.

### Outcome A: modified constraint response

If `K_red` has only the GR scalar constraint structure but its source response differs from GR,

\[
\mu_{\rm BQG}\ne1
\quad\text{and/or}\quad
\Sigma_{\rm BQG}\ne1,
\]

without an additional stable physical pole, then the effect is **modified gravity**, not a new dark particle.

### Outcome B: additional physical scalar pole

If after physical projection and gauge reduction

\[
\det K_{red}(\omega,k)=0
\]

contains an additional branch not gauge-equivalent to the GR constraint sector, then it becomes a candidate emergent dark degree of freedom. It must still pass:

- positive physical residue;
- no ghost;
- no gradient instability;
- sufficiently small effective sound speed for structure formation if it is to mimic cold dark matter;
- consistent coupling to the same Weyl potential used in lensing;
- refinement and regulator stability.

A nonzero finite `P H_L P` matrix element alone does not satisfy these conditions.

---

## 7. Background / dark-energy branch

The homogeneous background must be obtained independently from the same physical effective action,

\[
\Gamma_{phys}[g]\big|_{FLRW}
=\Gamma_{FLRW}[a,N].
\]

Its lapse variation supplies the effective history density and its scale-factor variation supplies pressure. Covariant conservation then gives

\[
\boxed{
w_{hist}(a)
=-1-\frac13\frac{d\ln\rho_{hist}}{d\ln a}.}
\]

A constant term is vacuum-like only if the physical projector/history calculation actually produces it. The static local `Gamma_V` and the finite `S4` curvature reconstruction constant are not cosmological constants.

---

## 8. Connection to the TT branch

The TT and scalar sectors must share the same physicalization front end:

\[
\boxed{
\{C_A\}
\to\mathbb M
\to P_{phys}
\to Z[J_g]
\to W[J_g]
\to\Gamma[g].}
\]

Only the projections differ:

```text
TT projection      -> gravitational-wave poles and six quartic Wilson coefficients
scalar reduction   -> Phi/Psi, growth, lensing and possible dark scalar poles
FLRW restriction   -> rho_hist(a), p_hist(a), H(a)
```

This is scientifically important: a future scalar result may not use one history weighting while the TT result uses another.

---

## 9. Current exact, finite and open nodes

### Exact / reusable now

- q=2 logical shape and metric Jacobian;
- local relational-source positive-control algebra;
- exact local `Gamma_shape`;
- collective `j=1` nontrivial volume carrier;
- DeWitt conformal direction/signature control;
- real sine-ordered `K=[V,H_E]` amplitude machinery;
- covariant `C(V)` and `C(K)` legs;
- one real ordered `K-K-V` triple;
- scalar rank-channel selection rules;
- noncommuting 24-term epsilon assembler identity;
- exact finite master-projector theorem;
- dual lapse cochain and HDA support geometry;
- shared-face `S3` six-amplitude reduction;
- tetrahedral `k^2/k^4` moment tensors;
- complete six-dimensional TT quartic dictionary.

### Being computed in the current branch

- full 24-term sine-Lorentzian epsilon sum;
- direct logical-return witness `P H_L P`.

### Genuinely open after that

- full logical/collective Lorentzian matrix and Hermitian constraint convention;
- regulator-safe Lorentzian HDA on the same graph-changing habitat;
- full finite physical constraint family and zero-sector projector on that habitat;
- physical boundary/history state and source insertion prescription;
- connected shared-face scalar history amplitudes;
- physical lapse/shift/volume/shape scalar Hessian;
- discrete gauge-invariant/Bardeen map;
- universal conserved matter-source coupling;
- `rho_hist(a)`, `Phi`, `Psi`, `mu_BQG`, `Sigma_BQG` from one frozen output.

---

## 10. Falsification discipline

The scalar programme fails or must change architecture if any of the following persists under the declared refinement sequence:

1. the completed Lorentzian constraint cannot be made compatible with the already-frozen HDA/DeWitt principal structure;
2. the physical projector has no controlled refinement/rigging limit;
3. connected scalar interblock source response remains zero, so no physical spatial scalar kernel forms;
4. derivative-order `k^2` anisotropy survives in the leading metric cone;
5. the reduced scalar kernel contains ghosts/gradient instabilities;
6. dynamics and lensing require different source normalizations;
7. a proposed DM-like pole cannot cluster or lens consistently;
8. a proposed DE-like background term is not produced by the same physical `Gamma`.

The goal is not to force dark matter or dark energy to appear. The goal is to ask the completed BQG physical history what scalar response it actually generates.