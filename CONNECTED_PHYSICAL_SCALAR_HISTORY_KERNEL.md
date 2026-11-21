# Connected physical scalar history kernel: separating projection depth, relational time and spatial transfer

Status: **production architecture / no-shortcut correction. No physical BQG omega or interblock coefficient is claimed here.**

## 1. Three parameters that must never be conflated

The current repository contains three mathematically useful but physically different parameters:

1. master heat-kernel depth
   \[
   \tau
   \quad\text{in}\quad
   e^{-\tau\mathbb M};
   \]
2. constraint-resolvent parameter
   \[
   z
   \quad\text{in Feshbach / Krylov resolvents};
   \]
3. physical or relational history separation
   \[
   \Delta t_{rel}.
   \]

Only the third may eventually be Fourier-conjugate to a physical frequency.

Therefore

\[
\boxed{
\tau\ne t_{phys},
\qquad
z\ne\omega,
\qquad
\omega\leftrightarrow\Delta t_{rel}\text{ only after a physical clock/history map is derived.}
}
\]

This distinction is mandatory for both the scalar and TT programmes.

## 2. Physical history state / amplitude first

Let the regulated graph-changing constraints define the physical history weighting through an exact or asymptotic master/rigging construction.

The boundary q=2 carrier `B` is mapped to a dressed physical history sector. Symbolically,

\[
B
\xrightarrow{\text{rigging / projector}}
\mathcal H_{phys}^{boundary/history}.
\]

A relational clock or boundary-history label must then be derived or frozen independently of the projector depth.

The existing C8 Page--Wootters construction is an exact architectural positive control for this ordering, but its declared `R=J` step is not the production gravitational evolution operator.

## 3. Relational scalar observables

Use gauge-invariant/relational scalar observables on coarse block `b`, for example

\[
\mathcal O_A(t,b),
\]

where the physical scalar carrier may include

- collective volume/conformal information;
- transported scalar shape/shear information;
- matter/probe density variables when a conserved source is included;
- any additional gauge-invariant scalar combination derived from the physical constraint reduction.

Lapse and longitudinal shift are not inserted here as propagating physical observables. They belong to the off-shell/background-field constraint reduction described separately in `SCALAR_SOURCE_GAUGE_ORDERING_CORRECTION.md`.

## 4. The central connected object

For a physical history state/density matrix, define the connected relational two-point kernel

\[
\boxed{
C_{AB}(\Delta t;b,c)
=
\langle
\delta\mathcal O_A(t,b)
\,\delta\mathcal O_B(t+\Delta t,c)
\rangle_{phys,c}
}
\]

with the appropriate time ordering / relational ordering frozen by the physical history construction.

For noncommuting equal-history source insertions, the zero-source Hessian of a simple exponential source functional yields the corresponding symmetrized covariance. For genuine frequency dependence, however, the ordering at separated relational history labels must be derived; a static symmetrized covariance is not enough.

## 5. Spatial nearest-neighbor reduction already supplied by the repository

For a canonical shared-face pair, the exact stabilizer is `S3`. The reciprocal parity-even six-edge transfer decomposes into two real symmetric `2 x 2` multiplicity matrices,

\[
T_{A_1}^{even}
=\begin{pmatrix}a_1&m_1\\m_1&f_1\end{pmatrix},
\]

\[
T_E^{even}
=\begin{pmatrix}a_E&m_E\\m_E&f_E\end{pmatrix}\otimes I_E.
\]

Thus one canonical nearest-neighbor physical source correlator requires only six real transfer functions for each relational-frequency/shell argument.

The `S4` tetrahedral action transports that canonical pair to the remaining three faces.

## 6. Momentum transform

For translationally/coarse-homogeneous relational data, define

\[
C_{AB}(\Delta t,\mathbf k)
=
\sum_{\delta}
C_{AB}(\Delta t,\delta)
\,e^{i\mathbf k\cdot\mathbf r_\delta}.
\]

Reciprocity gives an even cosine expansion.

Because the tetrahedral neighbor normals satisfy

\[
\sum_a n_a^in_a^j=\frac43\delta^{ij},
\]

the leading scalar nearest-neighbor term is isotropic:

\[
\boxed{
C_{AB}(\Delta t,\mathbf k)
=C^{(0)}_{AB}(\Delta t)
+C^{(2)}_{AB}(\Delta t)k^2
+O(k^4).
}
\]

The microscopic coefficient is still a physical-history output. Symmetry fixes its tensor form, not its value.

## 7. Relational-frequency transform

Only after the history label has a physical relational interpretation define

\[
\boxed{
C_{AB}(\omega,\mathbf k)
=
\sum_{\Delta t}
e^{i\omega\Delta t}
C_{AB}(\Delta t,\mathbf k)
}
\]

or the corresponding continuum transform.

A Euclidean/imaginary-time construction would require an independently justified analytic-continuation prescription. No such continuation is implied by the master heat kernel.

## 8. From connected response to the 1PI scalar kernel

On the non-gauge physical scalar support,

\[
\Gamma_{scalar}^{(2)}(\omega,k)
=C_{phys}^{-1}(\omega,k)
\]

with the appropriate matrix inverse/pseudoinverse on the frozen support.

In the off-shell ADM route, the corresponding background-field Hessian must first be reduced by the Hamiltonian/diffeomorphism constraints; the final gauge-invariant kernel should agree with the relational route on common observables.

This gives a powerful future internal consistency test:

\[
\boxed{
\Gamma_{scalar,\,relational}^{(2)}
\stackrel{?}{=}
\Gamma_{scalar,\,Dirac/Schur}^{(2)}
}
\]

after matching conventions and physical normalization.

## 9. Dark-sector discriminator

Once the physical kernel exists, separate three cases.

### Constraint-like long-range response

\[
K_{red}(0,k)\sim Z_k k^2
\quad\Rightarrow\quad
K_{red}^{-1}(0,k)\sim\frac1{Z_k k^2}.
\]

This can modify gravity without adding a propagating scalar particle.

### Genuine additional physical mode

A new branch requires

\[
\det K_{red}(\omega,k)=0
\]

after physical projection and gauge reduction, with positive residue and stable dispersion.

### Short-range analytic response

\[
K_{red}(0,k)=K_0+O(k^2),
\qquad K_0\ne0,
\]

gives no long-range Poisson enhancement.

## 10. What a cold-dark-matter interpretation would additionally require

An extra scalar branch is not automatically dark matter. It must also produce, from the same history kernel,

- positive physical norm/residue;
- no gradient instability;
- sufficiently small effective sound speed on structure-forming scales;
- correct background dilution or abundance;
- dynamics and lensing through the same metric response;
- regulator/refinement stability.

## 11. Dark-energy/background sector remains the zero-momentum restriction of the same theory

The homogeneous branch must be derived from the same physical history measure, not from a separate fitted function:

\[
\Gamma_{phys}[g]\big|_{FLRW}
=\Gamma_{FLRW}[a,N].
\]

Its variations determine `rho_hist(a)` and `p_hist(a)`. A volume-like term can be called vacuum-like only if the physical history actually generates the required homogeneous contribution.

## 12. Minimal next production data structure

The next real scalar-history computation should therefore store a tensor with explicit labels

```text
projector/refinement level
relational history separation Delta t
source carrier A,B
block separation / shared-face class
S3 multiplicity channel
```

and only later transform

```text
Delta t -> omega
block displacement -> k.
```

This prevents projection depth, constraint spectral variables and physical spacetime arguments from being mixed during data reduction.
