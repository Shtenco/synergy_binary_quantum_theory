# BQG unified system equation and proof ledger — 2026-11-11

## Purpose and status discipline

This document joins the principal finite constructions into one operator chain. It is a dependency map and a proof ledger, not a claim that the full physical theory has been derived.

Status labels:
- **Exact finite identity:** follows algebraically on explicitly declared finite-dimensional spaces.
- **Finite numerical evidence:** computed or tested only for a stated carrier, cutoff, and floating-point representation.
- **Conditional theorem:** valid only under listed analytic/physical hypotheses.
- **Open:** the necessary operator, limit, artifact, or independent verification is absent.
- **No-go on a declared carrier:** a negative result for that carrier only; not automatically a theorem about the full/refined theory.

## 1. Kinematics and finite regulator

For an oriented graph \(\Gamma=(V,E)\), a finite Peter–Weyl link cutoff \(J\) gives

\[
\mathcal H^{\mathrm{kin}}_{\Gamma,J}
=
\bigotimes_{e\in E}
\left(
\bigoplus_{j\le J} V_j\otimes V_j^*
\right).
\]

The gauge-invariant carrier is

\[
\mathcal H^{\mathrm{GI}}_{\Gamma,J}
=
\bigoplus_{\{j_e\le J\}}
\bigotimes_{v\in V}
\operatorname{Inv}_{SU(2)}
\left(
\bigotimes_{e\ni v}V_{j_e}^{\,\varepsilon(v,e)}
\right),
\]

with orientation-dependent dual representations encoded by \(\varepsilon(v,e)\). In the K5 implementation, basis keys carry the ten link spins and five local recoupling labels; local intertwiners are assembled from Clebsch–Gordan coefficients and the declared dimension normalization.

This is a finite carrier. It is not the unbounded Peter–Weyl Hilbert space and not a continuum limit.

## 2. Local volume and node operators

For a declared ordered triple of incident links, the raw oriented triple product is

\[
Q_v^{\mathrm{raw}}
=
\epsilon_{abc}J_{e_1}^{a}J_{e_2}^{b}J_{e_3}^{c}.
\]

The production code records the raw anti-Hermitian defect before constructing

\[
Q_v^{H}
=
\frac12(Q_v^{\mathrm{raw}}+Q_v^{\mathrm{raw}\dagger}),
\qquad
V_v=\sqrt{|Q_v^{H}|}.
\]

The finite spectral definition is

\[
Q_v^H=U\operatorname{diag}(\lambda_i)U^\dagger,
\qquad
V_v=U\operatorname{diag}(\sqrt{|\lambda_i|})U^\dagger.
\]

The symmetrization correction obeys the exact finite-matrix estimate

\[
\|Q_v^{\mathrm{raw}}-Q_v^H\|_2
\le
\|Q_v^{\mathrm{raw}}-Q_v^H\|_F
=
\frac12\|Q_v^{\mathrm{raw}}-Q_v^{\mathrm{raw}\dagger}\|_F.
\]

This bounds the correction to \(Q\), not the nonlinear map \(Q\mapsto\sqrt{|Q|}\). The recorded scalar derivative \(1/(2\sqrt{|\lambda|})\) is not an operator-Lipschitz or eigensolver backward-error bound.

The K5 node operator is assembled from oriented path-hit monomials and the local volume action. In the reviewed implementation the Hermitian node combination has the structure

\[
H_v^{(J,\epsilon)}
=
\frac12
\sum_{\alpha\in\mathcal S_v}
s_\alpha
\left(
T_{v,\alpha}^{(J,\epsilon)}
+
T_{v,\alpha}^{(J,\epsilon)\dagger}
\right),
\]

where each \(T\) is a finite signed sum of path-hit/volume monomials. The finite combinatorics, selected primitive adjoints, and selected projected node blocks have independent tests. Full-domain self-adjointness and full quantum HDA closure do not follow from those selected blocks.

## 3. Pruning is part of the represented operator

Write the implemented constraint family as

\[
H_v^{(J,\epsilon,\mathrm{mode})}.
\]

The legacy path has multiple amplitude thresholds at different accumulation/composition layers. The explicit no_prune path bypasses the reviewed amplitude-threshold deletions and uses coherent aggregation, but it still uses finite \(J\), floating-point CG coefficients, tensor contractions, and finite-precision summation. Therefore

\[
\texttt{no\_prune}
\ne
\text{symbolic exact arithmetic}
\ne
J=\infty.
\]

The independent pruning helper defines a policy and a discarded-vector ledger. For each stage \(s\), it can record the discarded vector \(\delta_s\), its \(\ell^1/\ell^2\) norms, count, maximum amplitude and digest. A valid downstream error bound has the form

\[
\|\Delta y\|
\le
\sum_s
\|A_{n-1}\cdots A_{s+1}\|\,\|\delta_s\|.
\]

The helper is not yet wired through every production K5 cache/composition layer, so the actual production path does not yet emit a complete stage-by-stage pruning error certificate.

## 4. One global constraint map and the master operator

For \(m=|V|\) node constraints, define the direct-sum map

\[
\mathcal C_{\Gamma,J,\epsilon}:
\mathcal H^{\mathrm{GI}}_{\Gamma,J}
\longrightarrow
\mathbb C^m\otimes\mathcal K,
\qquad
\mathcal C\psi
=
\sum_{v=0}^{m-1}|v\rangle\otimes H_v\psi.
\]

Equivalently, \(\mathcal C\psi=(H_0\psi,\ldots,H_{m-1}\psi)\). For a Hermitian strictly positive node metric \(G\succ0\),

\[
\boxed{
M_G
=
\mathcal C^\dagger(G\otimes I)\mathcal C
}
\]

and

\[
\langle\psi,M_G\psi\rangle
=
\|(G^{1/2}\otimes I)\mathcal C\psi\|^2.
\]

Since \(G^{1/2}\otimes I\) is invertible,

\[
\boxed{
\ker M_G
=
\ker\mathcal C
=
\bigcap_v\ker H_v.
}
\]

This is an exact finite-dimensional theorem for the represented matrices. It fails in this form for singular positive-semidefinite \(G\), because a singular metric can erase constraint combinations.

The Loewner bounds are

\[
\boxed{
\lambda_{\min}(G)\,\sigma_i(\mathcal C)^2
\le
\lambda_i(M_G)
\le
\lambda_{\max}(G)\,\sigma_i(\mathcal C)^2.
}
\]

Thus changing a strictly positive \(G\) can reshape the nonzero spectrum and conditioning, but cannot change the exact common kernel.

## 5. Shell compression must retain coherent leakage

For an orthogonal shell projector \(P\), set

\[
C_P=PHP,
\qquad
L_P=(I-P)HP.
\]

Then the exact identity is

\[
\boxed{
PH^\dagger HP
=
C_P^\dagger C_P+L_P^\dagger L_P.
}
\]

The leakage Gram contains off-diagonal interference between different input columns that reach the same outside basis state. Replacing it by separate column norms or a diagonal leakage penalty changes the operator.

The bqg_sparse_leakage_gram.py helper is now wired through the downstream bqg_k5_graph_sector_master_shell.py wrapper. For each node, the wrapper obtains complete no-prune K5 output columns, partitions them into shell action C_v and outside-shell leakage, and assembles the coherent G_v^leak=L_v^dagger L_v. The sparse-pruning workflow now includes the wrapper and a seven-method synthetic dense-oracle test file. This is a real integration step, but the tests patch the production action with synthetic columns; no full production shell calculation has been execution-confirmed. The helper and wrapper still need scaled-overflow protection and a checked orthonormal-basis precondition.

## 6. Direct singular values, rank, and error certification

For a tall sparse constraint map \(\mathcal C\), the preferred numerical object is its direct singular spectrum, not a rank conclusion inferred only from \(\mathcal C^\dagger\mathcal C\):

\[
\mathcal C=U\Sigma V^\dagger,
\qquad
\kappa(\mathcal C^\dagger\mathcal C)=\kappa(\mathcal C)^2.
\]

The current streaming implementation uses Householder QR in a binary-counter TSQR tree and an SVD of the final compact \(R\). It exposes the complete right singular basis and flags common-scale separation/underflow. It explicitly does not certify exact rank because it has no upstream construction-error bound.

If \(\widehat{\mathcal C}=\mathcal C+E\), Weyl's singular-value inequality gives

\[
|\sigma_i(\widehat{\mathcal C})-\sigma_i(\mathcal C)|
\le\|E\|_2.
\]

Therefore, only with a justified \(\eta\ge\|E\|_2\) and

\[
\boxed{\sigma_{\min}(\widehat{\mathcal C})>\eta}
\]

can full column rank of the exact target map be certified. The missing upstream ledger must cover CG/coefficient evaluation, intertwiner construction, volume spectral reconstruction, branch aggregation, sparse assembly, and all pruning stages. TSQR backward error belongs in a separate numerical-solver term.

For adjacent cutoffs, the code compares

\[
\Delta\mathcal C_{J,J'}
=
\mathcal C_{J'}-\mathcal C_J
\]

in a common labelled row embedding, zero-extending missing canonical output-state labels. The direct TSQR/SVD norm is a finite-float64 comparison. The sufficient Weyl condition

\[
\|\Delta\mathcal C_{J,J'}\|_2
<
\sigma_{\min}(\mathcal C_J)
\]

protects full column rank of the represented maps, but does not certify the exact underlying operators without the upstream error budget.

## 7. Physical projector and refinement are not yet one implemented map

On a finite carrier the exact spectral projector is

\[
\Pi_{0,\Gamma,J,\epsilon}
=
\mathbf1_{\{0\}}(M_G).
\]

A finite tolerance window

\[
\Pi_{\delta}
=
\mathbf1_{[0,\delta]}(M_G)
\]

is a different object. It is not automatically an exact kernel projector, a rigging map, or a physical inner product.

A continuum/refinement construction requires graph embeddings \(I_{\Gamma\to\Gamma'}\), compatible output embeddings \(\mathcal J_{\Gamma\to\Gamma'}\), and a controlled intertwining defect, for example

\[
\mathcal C_{\Gamma',J',\epsilon'}I_{\Gamma\to\Gamma'}
-
\mathcal J_{\Gamma\to\Gamma'}\mathcal C_{\Gamma,J,\epsilon}
\longrightarrow0
\]

in a specified operator topology and with a uniform error budget. This is a required target relation, not a theorem currently supplied by the code.

The actual 32-dimensional all-\(j=1/2\), \(K\in\{0,2\}^5\), five-node, \(J_{\max}=5/2\) pilot reported no numerical zero sector and a large positive identity-metric minimum eigenvalue. This is strong finite-carrier numerical evidence for an empty common kernel on that declared represented carrier. It is not a no-go theorem for larger/refined Hilbert spaces, and the finite floating-point/pruning construction error is not yet bounded rigorously. The exact positive-metric kernel theorem means that, for the exact represented finite map, no other strictly positive \(G\) can manufacture a zero mode.

## 8. Lorentzian dynamics and the HDA gate

The full target is an off-shell same-domain identity of the form

\[
\boxed{
[H[N],H[M]]
\stackrel{?}{=}
D\!\left[q^{ab}(N\partial_bM-M\partial_bN)\right].
}
\]

The old factorized Hamiltonian \(H_{\mathrm{geom}}[N]\otimes I_{\mathrm{path}}\) has an exact route-channel no-go: its commutator has no path-derivative component. The route-normal principal-symbol construction is a finite/semiclassical constructive result, but the full quantum gravitational HDA remains open.

For the Lorentzian sine scalar channel, the representation-theoretic selection rule from three rank-0/1 legs returning to scalar is exact:

\[
\boxed{000,\;011,\;101,\;110,\;111},
\qquad
J_2\in\{0,1\}.
\]

This authorizes deleting the forbidden intermediate \(J_2=2\) sector after the middle action. It does not authorize thresholding retained amplitudes. The accelerated 24-term full \(H_L\) still requires an unrestricted amplitude-by-amplitude reference, an independently checked Hermitian ordering, mixed commutators at their derived spin walls, and the full HDA comparison.

## 9. Physical TT and scalar observables

The finite S4 parity-even quartic TT response has a six-dimensional invariant space. Its exact extraction map is nonsingular. Under additional physical assumptions—local spatial rotation invariance, a locally Lorentz-invariant vacuum, a single gapless metric graviton and a local analytic derivative expansion—the six finite-regulator coefficients collapse to one isotropic TT coefficient, and a nonzero physical massless-pole \(k^4\) shift is excluded within that family.

This is a conditional EFT theorem. The physical BQG six-vector, physical projector, history measure, physical residue and common physical scale are not yet computed.

The scalar physical object must be a connected relational-history kernel,

\[
C^{\mathrm{phys}}_{AB}(\Delta t;b,c)
=
\left\langle
\delta\mathcal O_A(t,b)
\delta\mathcal O_B(t+\Delta t,c)
\right\rangle_{\mathrm{phys},c},
\qquad
\Gamma^{(2)}_{\mathrm{scalar}}
=
(C^{(2)}_{\mathrm{phys}})^{-1},
\]

with a derived relational clock and physical inner product. The exact q=2 local shape Legendre transform is a positive finite control; it does not supply the missing conformal/lapse direction or connected inter-block history dynamics.

## 10. The single end-to-end target equation

All pieces must ultimately be joined into one regulator-labelled object:

\[
\mathfrak B_{\Gamma,J,\epsilon}
=
\left(
\mathcal H^{\mathrm{GI}}_{\Gamma,J},
\{H_v^{(\Gamma,J,\epsilon)}\}_{v\in V},
G_\Gamma,
\mathcal C_{\Gamma,J,\epsilon},
M_{G,\Gamma,J,\epsilon},
\Pi_{\mathrm{phys}},
I_{\Gamma\to\Gamma'},
\mathcal R_{\mathrm{phys}},
\mathcal O_{\mathrm{phys}}
\right).
\]

The physical predictions are not the finite spectra by themselves. They are the limits of normalized physical observables after the physical map, graph refinement, spin removal, pruning control, and relational history construction have been made compatible:

\[
\boxed{
\mathcal A_{\mathrm{BQG}}[\mathcal O_1,\ldots,\mathcal O_n]
=
\lim_{\substack{
\Gamma\nearrow\Gamma_\infty\\
J\to\infty\\
\epsilon\downarrow0
}}
\mathcal A^{\mathrm{phys}}_{\Gamma,J,\epsilon}
[\mathcal O_1,\ldots,\mathcal O_n].
}
\]

The existence, path independence, normalization and physical interpretation of this limit are open. The formula is the system-level target, not a result already computed by the repository.

## 11. Current proof ledger

| Layer | What is established | What is not established |
|---|---|---|
| Finite SU(2)/Peter–Weyl representation | Exact finite coupling and selection identities; selected basis/adjoint oracles | Uniform all-sector/refinement theorem |
| K5 volume | Raw-defect diagnostics, finite permutation covariance, independent finite SVD/eigenspectral oracles | Operator-norm propagation through \(\sqrt{|Q|}\); full-domain error bound |
| K5 node constraints | Selected finite primitive and projected adjointness tests | Full same-domain quantum self-adjointness and HDA closure |
| Sparse pruning | Deterministic coherent aggregation helper; legacy/no-prune modes | Production-wide stage ledger and certified accumulated error |
| Leakage | Coherent Gram integrated through graph-shell wrapper; seven synthetic dense-oracle tests wired into workflow | Production-shell execution, overflow-safe accumulation, verified orthonormal basis and cutoff completeness |
| Master constraint | Exact \(G\succ0\) kernel theorem and Loewner bounds | Physical projector/rigging map on refined carrier |
| Direct spectrum | TSQR/direct-SVD and dense synthetic oracles | Upstream operator-construction error certificate |
| Actual 32D carrier | Strong finite numerical empty-kernel evidence | General no-go for enlarged/refined spaces |
| S4-sign rank | Exact selector/multiplicity reductions and shard smoke | Complete 128-shard/global rank certificate |
| Lorentzian sine/HDA | Exact scalar-channel selection rule; route-factorization no-go | Full Hermitian \(H_L\), mixed commutators, full quantum HDA |
| TT | Exact finite six-dimensional invariant quotient; conditional EFT collapse theorem | Physical BQG six-vector and physical TT pole |
| Scalar/history | Exact local positive-control source | Connected physical scalar kernel and relational time |
| Continuum/experiment | Explicit target and dependency map | Joint regulator limit, calibrated observables and experimental confirmation |

**Scientific rule:** no finite numerical gate, successful workflow, representation identity, or conditional EFT theorem may be promoted into a claim of complete physical BQG closure without the missing maps and limits above.


## 12. Graph-shell integration, cutoff completeness, and numerical hardening — source recheck

The graph-shell wrapper now assembles the finite pullback from complete supplied K5 node columns:

\[
M_{\mathrm{projected}}=\sum_v C_v^\dagger C_v,
\qquad
M_{\mathrm{leak}}=\sum_v L_v^\dagger L_v,
\qquad
M_{\mathrm{full}}=M_{\mathrm{projected}}+M_{\mathrm{leak}}.
\]

Here \(C_v=P H_v P\) and \(L_v=(I-P)H_vP\). The leakage Gram is assembled coherently across columns sharing the same outside-shell key, so its off-diagonal interference terms are retained. The workflow runs the graph-shell wrapper's 12-method synthetic test suite; this is source configuration, not proof of a fresh Actions PASS.

### Spin-wall enforcement — repaired in source

The earlier audit text said that the wrapper accepted depth and cutoff independently. That statement is stale. The current source defines

\[
J_{\max}^{(2)}\ge 1+2(d+1)=2d+3
\]

and enforces it at the start of `run()`, before shell generation. The bound follows from the T-monomial edge-hit multiplicities \((1,1,1,2)\): one node-H can raise doubled spin by at most two, and the pullback applies one additional H after the generated depth. The default depth 1 / jmax2 5 meets the bound; depth 2 / jmax2 5 is rejected, while depth 2 requires jmax2 at least 7.

The regression `test_depth_two_requires_complete_pullback_spin_wall` checks that the invalid depth/cutoff pair fails before `generate_shell` is called. This closes the specific fail-open regulator-wall issue for the current shell-generation rule. It does not prove the hit-count model for arbitrary future operator words; changes to T-monomial combinatorics must update the bound and its independent reachability oracle together.

### Basis metric prerequisite — finite shell check now present

The earlier statement that the wrapper did not check shell-basis normalization is also stale. The current source calls `validate_shell_basis_orthonormality(basis)` on the generated basis. It groups every local \((v,\text{spins})\) sector actually used, constructs the corresponding oriented intertwiner tensors, measures each local Gram defect, and fails closed if the declared tolerance is exceeded. The result is included in `structural_ok`. A dedicated test checks the initial Peter–Weyl carrier.

This is a finite numerical check over all local sectors represented in the particular generated shell. It is not a general proof for all possible graph bases, all spin cutoffs, or alternate normalization conventions. The exact global product-basis interpretation still relies on the declared Peter–Weyl orthogonality and tensor-product conventions.

### Overflow and PSD diagnostics — partial hardening present

The current leakage helper validates finite input amplitudes and checks the assembled shell action, leakage Gram, shell Gram and pullback for non-finite values. It raises `FloatingPointError` if the Gram cannot be represented in float64; a regression deliberately supplies an amplitude whose square overflows. The PSD classifier now uses half-before-addition Hermitianization,

\[
A_H=\tfrac12 A+\tfrac12 A^\dagger,
\]

and validates finite, nonnegative, non-boolean tolerance controls. A near-maximal finite Hermitian one-by-one matrix is covered by a regression.

The helper still uses ordinary float64 outer products and accumulation rather than a scaled/compensated Gram algorithm. Its present safety property is **fail-closed on detected non-finite results**, not scale-free representability for every possible finite input. Any future attempt to support more extreme dynamic range should add scaled accumulation and compare against an independent dense/high-precision oracle.

### What remains open

- No production-generated multi-state shell run has been freshly confirmed in this audit turn; the current workflow's tests are synthetic dense-oracle fixtures plus selected finite-sector production primitive tests.
- The shell pullback is exact only relative to the supplied finite float64 columns and the checked basis conventions. Upstream CG/coefficient, volume square-root, aggregation, cutoff and solver errors still lack one propagated global bound.
- A full physical projector/rigging map, cutoff-independent graph embedding, refinement intertwining estimate and continuum limit remain unimplemented.
- The K5 no-prune label means removal of amplitude-threshold pruning on the reviewed path; it does not mean exact arithmetic, infinite spin, or a completed physical HDA.

## 13. Positive-metric bound implementation: dimension and underflow contract

The finite Loewner theorem in §4 is now reflected in the constraint-stability implementation with explicit dimensions:
- the node metric dimension is checked against the number of node constraints;
- the singular-spectrum length is checked against the input carrier dimension;
- TSQR's full right-nullspace convention pads the singular spectrum with exact zeros when the output row count is smaller than the input dimension;
- underflow in \(\sigma_i^2\), \(\lambda_{\min}(G)\sigma_i^2\), and \(\lambda_{\max}(G)\sigma_i^2\) is reported separately and makes the numerical bound set unresolved.

This is a correctness/diagnostic hardening of the finite theorem implementation. It does not turn the singular values into exact symbolic values, nor does it supply the missing upstream construction-error bound \(\eta\) required for an exact-operator Weyl certificate. The regression file currently contains 40 test methods by static source count. The metric SPD guard uses a scale-free relative minimum-eigenvalue gap, so a tiny but well-conditioned positive metric is not rejected merely because an absolute threshold underflows; no fresh CI PASS is claimed.
