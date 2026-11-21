# BQG spectral-history graph: exact finite cyclic closure of the master history

Status: **finite spectral-history closure theorem + hardened fail-closed implementation. The spectral layer itself is not authorized to emit a physical BQG projector; physical promotion requires the independent evidence-linked certificate in `BQG_PHYSICAL_SPECTRAL_HISTORY_CERTIFICATE.md`.**

## 1. Purpose

The useful mathematical part of spectral-graph models is the operator/spectral machinery itself:

- positive self-adjoint graph/master operators;
- spectral measures;
- heat kernels;
- zeta functions and return probabilities;
- Krylov/Lanczos compression;
- reconstruction of matrix functions from a cyclic spectral measure.

BQG does **not** import a guessed Watts--Strogatz graph, a chosen node count `N`, a prescribed degree `K`, a nonlocal-link probability `p`, fitted interaction-channel weights, or an identification of a master-constraint eigenvalue with a particle mass.

Instead, the spectral graph is **derived from the already-defined BQG master constraint**.

The finite spectral chain is

```text
actual regulated BQG master M_BQG >= 0
        -> declared boundary/source seed V
        -> block Krylov cyclic subspace K(M_BQG,V)
        -> block-Jacobi / spectral quotient J_V
        -> matrix-valued spectral measure dSigma_V(lambda)
        -> V^dag f(M_BQG) V
        -> candidate zero spectral weight / finite projector content
```

Physical promotion is a separate chain:

```text
finite spectral result
  + complete full-master/HDA certificate on the same habitat
  + refinement/rigging-map certificate
  -> PHYSICAL_PROJECTOR_HISTORY certificate
  + source-dressed convergence certificate
  -> connected W[J] -> Gamma_BQG
```

No dark-matter, dark-energy, mass, coupling, or cosmological datum enters either construction.

---

## 2. Finite cyclic spectral-closure theorem

Let `M=M^dag >= 0` be a finite regulated master constraint and let the seed block `V` have orthonormal columns.

Define the cyclic block-Krylov space

\[
\mathcal K_r(M,V)=\operatorname{span}\{V,MV,M^2V,\ldots,M^rV\}.
\]

Apply rank-revealing block Lanczos/Gram--Schmidt, obtaining orthonormal blocks `Q_n`. In a non-deflated presentation,

\[
M Q_n = Q_{n-1}B_n^\dagger+Q_nA_n+Q_{n+1}B_{n+1},
\]

with `A_n=A_n^dag` and off-diagonal blocks `B_n` determined by residual factorization.

The finite block-Jacobi matrix

\[
J_r=Q^\dagger M Q,
\qquad Q=(Q_0,\ldots,Q_r),
\]

is the **BQG spectral-history graph** on this seed. Its block vertices are Krylov layers, diagonal weights are `A_n`, and adjacent block-edge weights are `B_n`.

If the next residual vanishes,

\[
\boxed{B_{r+1}=0,}
\]

then `K_r` is invariant under `M`. Consequently, for every bounded Borel function `f` on the finite spectrum,

\[
\boxed{
V^\dagger f(M)V
=E_0^\dagger f(J_r)E_0.
}
\]

Therefore the finite regulated seed history closes **exactly**:

\[
\boxed{
H_V(\sigma)
=V^\dagger e^{-\sigma M}V
=E_0^\dagger e^{-\sigma J_r}E_0,
}
\]

and the zero-projector overlap is

\[
\boxed{
G_{V,0}
=V^\dagger\mathbf 1_{\{0\}}(M)V
=E_0^\dagger\mathbf 1_{\{0\}}(J_r)E_0
=\lim_{\sigma\to\infty}H_V(\sigma).
}
\]

Here `sigma` is the master heat/projector-flow parameter. It is **not** relational proper time and is **not** the physical frequency `omega`.

This theorem says that the quotient contains the entire finite seed-visible functional calculus. It does not say that `M` has already been certified as the complete physical BQG master or that a continuum/refinement physical sector exists.

---

## 3. Matrix-valued spectral measure

Let

\[
J_r u_a=\lambda_a u_a.
\]

Define the positive matrix weights

\[
W_a=E_0^\dagger u_a u_a^\dagger E_0.
\]

Then the seed spectral measure is

\[
\boxed{
d\Sigma_V(\lambda)
=\sum_a W_a\,\delta(\lambda-\lambda_a)d\lambda,
}
\]

and

\[
V^\dagger f(M)V
=\int f(\lambda)d\Sigma_V(\lambda).
\]

In particular,

\[
\mu_n=V^\dagger M^nV
=\int \lambda^n d\Sigma_V(\lambda).
\]

Thus the BQG master moments are moments of the same derived spectral graph that supplies the finite heat history and candidate zero weight.

---

## 4. Moment/Hankel construction

A useful independent construction of the same finite cyclic quotient starts from certified moments.

For depth `r`, form

\[
\mathcal H_r=[\mu_{i+j}]_{i,j=0}^{r},
\qquad
\mathcal H_r^{(1)}=[\mu_{i+j+1}]_{i,j=0}^{r}.
\]

`H_r` is the Gram matrix of the raw block Krylov family

\[
(V,MV,\ldots,M^rV).
\]

Whitening the numerical support of `H_r` and projecting `H_r^(1)` produces a finite Hermitian quotient operator `T_r`. Its matrix-valued Gauss spectral measure reproduces the supplied Krylov moments.

An exact finite closure certificate can be obtained from rank stabilization,

\[
\boxed{
\operatorname{rank}\mathcal H_{r+1}
=\operatorname{rank}\mathcal H_r,
}
\]

provided an upstream numerical error certificate is supplied. A visually small singular value is not sufficient.

The preferred production certificate remains the direct block-Lanczos residual norm because it tests termination without reconstructing high powers.

---

## 5. Current Euclidean history frontier

For the orthonormal 32-state q=2 boundary seed `V0`, the Euclidean normal master has

\[
\mu_0=I_{32},
\qquad
\mu_1=V_0^\dagger M_EV_0=M_{EE}.
\]

The history-targeted producer computes

\[
Y=M_EV_0=\sum_v H_v^E H_v^E V_0
\]

in the frozen Hermitian Euclidean convention, after which

\[
\mu_2=Y^\dagger Y.
\]

The first genuine master-Krylov residual Gram is

\[
\boxed{
R_1
=\mu_2-\mu_1^\dagger\mu_1
=Y^\dagger(I-V_0V_0^\dagger)Y
\succeq0.
}
\]

If `R1=0` within a certified error bound, the Euclidean boundary cyclic space closes at depth zero and

\[
V_0^\dagger e^{-\sigma M_E}V_0=e^{-\sigma\mu_1}
\]

is exact **for that finite Euclidean regulated sector**.

If `R1` has positive rank, factor it to form `Q1`, evaluate `M_E Q1`, obtain `A1` and the next residual `B2`, and continue. No dense ambient Peter--Weyl matrix is required.

The historical production run did not finish all 32 `M_E b_i` columns, so the actual numerical `rank(R1)`, `Q1`, `B1`, `A1` and `B2` are not claimed yet. The calculation has been restarted on the recovered research branch.

---

## 6. Source-complete spectral graph

Closing only

\[
V_0^\dagger P_0V_0
\]

is not sufficient to claim a physical two-point kernel or `Gamma^(2)`.

For a declared source family `O_A`, the seed must contain every state needed by the desired source derivative. For a quadratic source claim a safe finite seed is generated from

\[
\boxed{
\mathcal V_{\Gamma^{(2)}}
=\operatorname{span}\{V_0,\ O_AV_0,\ O_AO_BV_0\}_{A,B\in\mathcal S},
}
\]

with the exact relational/time-ordered convention used by the existing source machinery.

A block spectral graph built from this source-enriched seed can supply all finite projector matrix elements needed by the zero-source amplitude, one-source derivative and quadratic source Hessian. Connected correlators still require the legal `W=log Z` construction and a source-dressed refinement certificate.

---

## 7. Spectral diagnostics imported without phenomenological fitting

Once the graph is derived rather than guessed, several diagnostics are legitimate.

### Seed return probability

\[
P_V(\sigma)
=\frac1{d_V}\operatorname{Tr}
\left(V^\dagger e^{-\sigma M}V\right).
\]

### Constraint spectral dimension

\[
\boxed{
d_{s,M}^{(V)}(\sigma)
=-2\frac{d\ln P_V}{d\ln\sigma}
=2\sigma\,
\frac{\sum_a\lambda_a e^{-\sigma\lambda_a}\operatorname{tr}W_a}
{\sum_a e^{-\sigma\lambda_a}\operatorname{tr}W_a}.
}
\]

This is a **master/constraint spectral-flow diagnostic**. It is not automatically the dimension of physical spacetime. Any identification with continuum spacetime spectral dimension requires an independent refinement/geometry theorem.

### Seed-weighted zeta function

On the positive finite spectral support,

\[
\boxed{
\zeta_V(s)
=\sum_{\lambda_a>0}
\operatorname{tr}(W_a)\lambda_a^{-s}.
}
\]

At finite regulator this is an ordinary finite sum; no analytic continuation is needed. It can diagnose refinement stability and scaling, but is not used to manufacture particle masses or fundamental constants.

---

## 8. What is explicitly rejected

The following moves are forbidden:

```text
choose N,K,p,xi to obtain desired physics;
fit graph-channel weights to observed forces;
select a spectral resonance because it matches a known particle mass;
identify sqrt(master eigenvalue) with particle mass;
interpret master heat sigma as physical time;
interpret a master resolvent parameter z as physical omega;
tune beta/lambda to create a zero mode or a dark sector;
add Krylov depth until a desired near-zero appears without a predeclared convergence test;
set inline booleans in a moment packet and call that a physical certificate.
```

Positive eigenvalues of `M_BQG` measure violation of the declared constraints at the finite regulator. Physical particle poles, masses and propagation laws may only be read later from the source-dressed physical effective kernel `Gamma_BQG^(2)(omega,k)`.

---

## 9. Hardened fail-closed emission rule

The spectral graph layer **never emits a physical BQG history/projector by itself**.

It may certify only

```text
finite_spectral_history_closed = true
```

and expose a

```text
candidate_zero_spectral_weight
```

when the finite cyclic quotient is closed.

Inline fields such as

```text
domain_complete
master_constraint_certified
quantum_hda_or_explicit_dtarget_certified
source_seed_complete_for_claim
```

are retained only as unverified declarations/backward-compatible diagnostics. Even if all are `true`, the strongest spectral-layer status is

```text
FINITE_SPECTRAL_HISTORY_CLOSED_PHYSICAL_CERTIFICATE_REQUIRED
```

and

```text
physical_history_closed = false
physical_projector_emitted = false.
```

Physical promotion is delegated exclusively to

```text
scripts/bqg_physical_spectral_history_certificate_gate.py
```

which independently binds the spectral result to the complete full-master/HDA artifact through matching

```text
habitat_hash
domain_hash
convention_hash
master_pencil_hash
```

and then additionally requires the refinement/rigging-map certificate defined in `BQG_PHYSICAL_SPECTRAL_HISTORY_CERTIFICATE.md`.

Therefore a successful finite spectral calculation cannot, by construction, close the repository-level `PHYSICAL_PROJECTOR_HISTORY` gate without independent evidence.

---

## 10. Status

The finite spectral-history **mathematical architecture is closed**:

\[
\boxed{
M\to\mathcal K(M,V)\to J_V\to d\Sigma_V
\to V^\dagger e^{-\sigma M}V
\to V^\dagger\mathbf 1_{\{0\}}(M)V.
}
\]

The actual physical BQG history is **not** closed. Its evidence ladder is now explicit:

```text
S0  finite spectral cyclic closure
S1  complete full-master + matching HDA + spectral hash linkage
S2  refinement/rigging projector and boundary-history convergence
S3  source-dressed connected-history convergence
```

`PHYSICAL_PROJECTOR_HISTORY` may close only at S2, and connected source/correlator claims require S3.

The current actual Euclidean BQG calculation is still before completed S0: `mu0` and `mu1` exist, the direct `mu2 -> R1 -> Q1,B1` producer exists, one historical real master-image shard was recovered, and the full 32-column calculation has been restarted. The production Lorentzian/full-master/HDA and refinement evidence remain additional independent requirements.

This is the intended scientific boundary: **spectral closure is now a strong exact computational theorem, not a shortcut around physicalization.**
