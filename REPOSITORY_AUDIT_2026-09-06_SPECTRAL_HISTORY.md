# Binary Quantum Gravity — spectral-history and claim audit

**Audit date:** 2026-09-06  
**Base main at audit start:** `0f048599f9d31d8b5aff33fa22288cfd6742abe0`  
**Audit branch:** `spectral-history-closure`

This note audits the active BQG theory after introducing an operator-derived spectral-history programme.  It is deliberately stricter than a project narrative: exact mathematics, finite numerical evidence, model-selection assumptions and open physical arrows are separated.

## 1. Executive scientific status

The repository contains a substantial internally consistent finite mathematical architecture.  Its strongest pieces are exact finite algebra/representation theory, explicit negative controls, finite Peter–Weyl constraint dynamics, Regge held-out controls, master-constraint projector mathematics and an unusually good fail-closed physicalization ledger.

The central unresolved issue is no longer whether useful finite structures can be built.  It is whether one **single microscopic constraint/history construction**, under controlled refinement, produces the physical connected metric effective action and its Lorentzian poles without inserting the desired continuum answer.

The strongest corrected status is therefore

```text
selected q=2 route-family mathematics          = exact in stated assumptions
selected spatial PL completion                 = exact/finite existence + tested stability
finite quantum-geometry carrier                = exact/finite
finite Peter-Weyl constraint spectral data     = exact/finite
finite master-projector theorem                 = exact
operator-derived finite spectral history       = now executable with exact stopping rule
unique microscopic selection of route rewrite  = NOT proved
continuum/refinement physical rigging history  = OPEN
connected physical W[J]                        = OPEN
physical Gamma^(2)(omega,k)                    = OPEN
physical six-Wilson vector                      = OPEN
external experimental confirmation             = NO
```

## 2. Correction of audit provenance

At the start of this audit the actual `main` HEAD was

```text
0f048599f9d31d8b5aff33fa22288cfd6742abe0
```

not a previously quoted `e6fc...` spectral-history commit.  The files `BQG_SPECTRAL_HISTORY_GRAPH.md` and `scripts/bqg_spectral_history_graph_gate.py` did not exist on `main` at that point.

The spectral-history work was therefore rebuilt explicitly on branch `spectral-history-closure` rather than treated as pre-existing evidence.

This correction is scientifically important: repository state, CI artifacts and mathematical claims must be independently verifiable rather than inherited from conversational history.

## 3. External spectral-graph theory: what survives audit

The phenomenological small-world spectral-graph framework considered during this audit contains useful mathematical technology but cannot be imported as physical evidence.

### Useful machinery

- spectral measures;
- heat kernels;
- spectral gaps;
- zeta functions when their domains/continuations are defined;
- spectral dimension as a diagnostic;
- refinement/RG comparison;
- Krylov/Jacobi operator reduction.

### Rejected inference pattern

The BQG theory must not infer physical constants by choosing coordination number, small-world probability, channel scales or weights and then interpreting agreement with target numbers as independent prediction.

In particular, the external derivation that identifies a three-dimensional lattice sum

\[
\sum_{n\in\mathbb Z^3\setminus\{0\}}|n|^{-2}
\]

with a finite multiple of `zeta(2)` requires an explicit regularization/analytic continuation.  As an ordinary lattice sum it diverges.  A later switch from a `zeta(4)` contribution to `zeta(2)` also requires a derivation rather than a numerical substitution.

BQG therefore retains the spectral apparatus but replaces the phenomenological graph by a graph **derived from its operator**.

## 4. New exact spectral-history principle

For a finite Hermitian/positive constraint operator `M` and seed block `Q0`, block Lanczos generates

\[
M Q_n=Q_{n-1}B_n^\dagger+Q_nA_n+Q_{n+1}B_{n+1}.
\]

The graph vertices, ranks and edge blocks are outputs of `M`, not chosen network parameters.

The exact finite stopping rule is

\[
\boxed{B_{r+1}=0.}
\]

If this occurs, the cyclic Krylov space is invariant and for every finite spectral function `f`,

\[
\boxed{Q_0^\dagger f(M)Q_0=E_0^\dagger f(J_r)E_0.}
\]

Consequently the reduced Jacobi graph reproduces, on the seed sector, the full finite heat history and zero spectral projector, not merely a few moments.

The deterministic implementation regression closes a synthetic reachable space and reproduces heat kernels/projector to floating-point errors of order `1e-15`.

## 5. Two spectral routes must not be conflated

### Route A — positive master constraint

Use

\[
\mathbb M_E=\sum_v(H^{sine}_{E,v})^\dagger H^{sine}_{E,v}\ge0.
\]

This is the route naturally aligned with the finite common-kernel theorem

\[
\ker\mathbb M_E=\bigcap_v\ker H_{E,v}
\]

for the stated finite constraint family.

The new first-edge gate computes `A0`, `mu2`, `B1^dag B1` and the actual residual rank on the complete 32-state logical K5 seed sector.

### Route B — parity Krylov chain of `H_E0+H_E1`

The repository already contains much deeper exact data for

\[
H=H^{sine}_{E,0}+H^{sine}_{E,1}.
\]

Exact parity gives

\[
P H P=0.
\]

The certified higher-shell artifact already establishes

\[
B_1^\dagger B_1=K=P H^2P,
\]

\[
B_2^\dagger B_2=\Lambda
=K^{-1/2}(P H^4P-K^2)K^{-1/2}
\]

with reconstruction errors near `1.6e-13`.

Thus BQG already possessed the first two edges of an actual operator-derived spectral graph; they were not previously organized as a complete history-closure programme.

The new B3 calculation extends the exact source columns from `H^2|i>` to `H^3|i>` and tests whether the next hopping block vanishes.

These two routes answer different questions.  A closed `H_E0+H_E1` parity chain is not automatically the common-kernel master projector of the complete constraint family.

## 6. q=2 audit: exact result versus selection assumption

The active theory uses a local route family with `q` binary route choices.  It imposes local valence homogeneity:

\[
\boxed{q+2=2^q.}
\]

For integers `q>=1`, the unique solution is indeed

\[
q=2.
\]

This is exact **given the route-family definition and the valence-homogeneity requirement**.

The requirement itself is an architectural/model-selection axiom.  It is not currently derived from the graph-changing Hamiltonian, a variational principle, an information-theoretic uniqueness theorem or external observation.

Therefore the strongest valid statement is

> `q=2` is uniquely selected inside the declared locally valence-homogeneous binary route family.

The stronger wording “binary microphysics uniquely implies q=2” would overstate the result.

## 7. Dimension-three audit: mathematically exact but not statistically independent

For the frozen route rewrite,

\[
B=2^q,
\qquad
E_{g+1}^{active}=2B E_g^{active},
\qquad
L_{g+1}=2L_g.
\]

Hence generally

\[
\boxed{d_*^{causal-volume}=\log_2(2B)=q+1.}
\]

With the selected `q=2`,

\[
\boxed{d_*=3.}
\]

This theorem is exact for the frozen rewrite.  However, the value `3` is mathematically downstream of the same route-family assumptions that selected `q=2`.

Likewise:

- `q=2 -> four labels -> C4`;
- suspension of `C4` gives an octahedral `S2` link;
- the selected global 16-cell boundary gives an `S3` PL completion;
- the route growth gives `d*=3`.

These are valuable mutually consistent consequences, but they are **not four statistically independent measurements of an unknown dimension**.  They share common microscopic assumptions.

Independent evidence would require a dimension-blind dynamical ensemble generated by the frozen microscopic operator/rule and a blind measurement of its large-scale topology/diffusion/volume exponents.

This qualification is consistent with the repository's own `BINARY_TO_GEOMETRY_GATE.md` and `HODGE_DIMENSION_SELECTOR.md`, both of which retain the stronger microscopic geometrogenesis arrow as open/conditional.

## 8. Global topology audit

The boundary of the 16-cell provides an exact economical PL `S3` completion with correct local link structure and tested refinement stability.

What is proved:

- existence of the selected completion;
- exact combinatorial/topological properties of that completion;
- finite refinement stability of the tested sequence.

What is not proved:

- uniqueness of `S3` under all admissible microscopic histories;
- dynamical dominance of the 16-cell completion in a physical path/history measure;
- suppression of alternative topology sectors under refinement.

Therefore `selected global PL 3-manifold` is an accurate label; `microscopic uniqueness of global spatial topology` is not.

## 9. Walsh/tetrahedral carrier audit

The `Z2^2` Walsh construction is one of the cleanest exact steps in the theory.

Three nontrivial real characters give four unit vectors with

\[
\sum_an_a=0,
\qquad
n_a\cdot n_b=-1/3\quad(a\ne b),
\]

which is exactly the Gram structure of regular tetrahedral normals.

Likewise the four spin-1/2 face carriers have an exact two-dimensional SU(2) Gauss-singlet sector.

The representation-theory statements are exact.  The remaining physical question is not whether the carrier exists, but whether the frozen microscopic history dynamically selects and propagates the corresponding semiclassical geometric sector with the required correlations.

## 10. Plebanski/Urbantke/Regge/ADM audit

The repository has strong finite bridges demonstrating that the selected geometric carrier can be mapped into familiar continuum gravitational structures and that negative controls fail as intended.

These controls establish compatibility and internal consistency.  They do not yet establish that the actual microscopic physical history has Einstein gravity as its unique IR effective action.

The important distinction is

```text
candidate carrier can realize GR structures
```

versus

```text
microscopic physical measure dynamically flows to GR with the correct residue and interactions.
```

Only the first is currently structurally closed.

The internal Regge `L=9,10` held-out test is a good genuine held-out test of the fixed-scaffold continuum correction law.  It remains an **internal numerical validation**, not a held-out observation of nature and not an upstream test of binary geometrogenesis.

## 11. HDA audit

Finite graph-changing HDA calculations are meaningful nontrivial consistency tests.  They show that the implementation is not merely a static geometry dictionary.

They do not yet prove:

- an anomaly-free constraint algebra uniformly over arbitrary refinements/habitats;
- existence of the continuum physical inner product;
- equivalence of the master-projector limit to a physical Lorentzian history;
- physical time.

The current fail-closed wording in `THEORY_STATUS.md` is correct on these points.

## 12. Master-projector audit

The finite theorem

\[
\mathbb M_G=C_A^\dagger G^{AB}C_B\ge0,
\qquad
\ker\mathbb M_G=\bigcap_A\ker C_A
\]

for positive definite `G` is exact.

When zero is isolated, the finite projector and heat convergence are exact as well.

The major continuum danger is gap closing.  If

\[
\Delta_M(\epsilon)\to0
\]

under regulator removal, the raw finite projector/heat convergence estimate is not enough.  The theory then needs a normalized spectral-window/rigging limit with convergence of physical matrix elements.

The new spectral graph helps because it exposes the seed spectral measure and gap flow directly, but it cannot remove this requirement.

## 13. Constraint spectrum is not physical frequency

This is a major strength of the present repository after the physicalization correction.

Neither

- Feshbach spectral `z`,
- Lanczos eigenvalues,
- C8 character angles,
- master gaps,
- higher-shell `Lambda` eigenvalues

may be renamed particle masses or physical `omega` without the history/effective-action bridge.

The required chain remains

```text
actual constraints
 -> physical projector / relational or boundary history
 -> source-dressed Z[J]
 -> connected W[J]
 -> gauge-reduced Gamma[g]
 -> physical Gamma^(2)(omega,k)
 -> TT/scalar/Maxwell response.
```

This guard should remain non-negotiable.

## 14. Scalar/cosmology/dark-sector audit

The recent q=2 scalar source/Legendre calculations are legitimate finite effective-action positive controls.

They are not yet cosmological scalar response because the theory still lacks, in one derived physical construction:

- physical gauge reduction;
- conserved matter-source coupling;
- connected spatial history;
- frequency/momentum-dependent physical kernel;
- homogeneous physical background effective action.

Therefore present `mu_BQG`, `Sigma_BQG`, dark-matter and dark-energy outputs correctly remain `OPEN_PHYSICAL` / `NOT_DERIVABLE` in the active ledgers.

Any future claim should be rejected if it uses a constraint eigenvalue or a local static source Hessian as a dark component by relabelling.

## 15. Experimental-status audit

No current internal calculation establishes experimental quantum-gravity confirmation.

The proper first external test must occur only after:

1. the physical history construction is frozen;
2. the physical TT/scalar/photon kernels are frozen;
3. all dimensionless outputs are fixed;
4. one common scale rule is fixed;
5. the external dataset/likelihood is preregistered;
6. no sector-by-sector retuning occurs after opening the data.

The existing repository truth ledger already encodes this philosophy and should be preserved.

## 16. Reproducibility risk discovered

The exact higher-shell B1/B2 matrices and sparse source columns currently depend on GitHub Actions artifacts from August 2026.  The repaired final artifact has a finite retention lifetime.

For the present B3 computation the artifacts are still available and their provenance is frozen by workflow run IDs plus per-column SHA256 recorded by the new gate.

Long-term repository hardening should preserve at minimum:

- final small-matrix certificates (`K`, `Lambda`, `B1`, `B2`, future `B3`);
- SHA256 manifest for every large sparse source artifact;
- generating commit SHA;
- exact script SHA;
- regulator declaration;
- a reproducible regeneration workflow.

Large raw sparse columns need not all live in Git, but their digest/provenance must survive artifact expiry.

## 17. Highest-value next calculations

The recommended order is now:

### A. Finish B3 for existing exact H01 parity chain

This determines whether the already-proven two-edge finite spectral history closes or grows another exact shell.

### B. Continue the positive-master chain

The positive `M_E` route is the one tied directly to the finite common constraint kernel/projector.  Continue `A1/R2/B2...` until closure or a clearly bounded truncation/refinement programme exists.

### C. Repeat across refinement/patch size

For each regulator/patch record:

- Jacobi block ranks;
- spectral measure;
- zero-sector weight;
- first positive gap;
- heat trace;
- constraint spectral dimension;
- terminator dependence when not closed.

Then ask whether normalized physical/projector matrix elements converge as the regulator changes.

### D. Only then build source-dressed physical history

Insert derived gauge-invariant metric sources into the projector/boundary history, compute connected `W[J]`, perform correct gauge/constraint reduction and only then read physical poles.

## 18. Final audit verdict

BQG is strongest when treated as a **candidate operator/geometry programme with exact finite certificates and explicit falsification gates**, not as a numerological theory of constants.

The new spectral-history direction improves it because it turns the vague question “how many history steps should be summed?” into a precise operator question:

\[
\boxed{\text{does the next block-Lanczos residual vanish?}}
\]

If yes, the finite seed spectral history is exactly closed.  If no, the operator itself demands another shell.

The remaining scientific leap is larger and separate: proving that the refinement-compatible physical projector/history and connected effective action generated from the actual BQG constraints have the observed Lorentzian gravitational IR behaviour.

That leap is **not yet closed**, and keeping it open is a strength of the audit rather than a weakness of the mathematics already established.
