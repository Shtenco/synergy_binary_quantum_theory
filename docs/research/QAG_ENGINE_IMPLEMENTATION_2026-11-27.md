# QAG-1 finite postselection engine — implementation checkpoint

**Date:** 2026-11-27  
**Repository:** Shtenco/synergy_binary_quantum_theory  
**Scope:** finite-dimensional quantum-channel benchmark only.

## Implemented files

- `qag_postselection.py` — NumPy implementation of a two-branch controlled source–probe interaction.
- `qag_periodic_probe.py` — periodic-grid position/momentum operators, Gaussian probe state, exact branch-dependent momentum-kick unitary, and weak-value comparison benchmark.
- `tests/test_qag_periodic_probe.py` — tests for Hermiticity, kick sign convention, weak-value sign reversal, and dephased convex-branch response.
- `tests/test_qag_postselection.py` — unit tests for pure/mixed source equivalence, coherence suppression, normalization, positivity, complete measurement probabilities, and zero-probability outcomes.
- `.github/workflows/qag-postselection-regression.yml` — isolated GitHub Actions workflow that compiles the module and runs the QAG-specific unittest suite.
- `docs/research/QAG_POSTSELECTION_ANTIGRAVITY_AUDIT_2026-11-27.md` — theoretical model and physical interpretation boundaries.

## Implemented channel

For source density matrix \(\rho_s\), rank-one postselection vector \(|\phi\rangle\), probe density matrix \(\rho_p\), and controlled branch unitaries \(U_L,U_R\), the code evaluates

\[
\widetilde\rho_p^\phi=
\sum_{i,j\in\{L,R\}}
(\rho_s)_{ij}\langle\phi|i\rangle
\langle j|\phi\rangle U_i\rho_pU_j^\dagger.
\]

It then calculates \(P_\phi=\operatorname{Tr}\widetilde\rho_p^\phi\) and normalizes only when the probability is numerically resolved. The pure-source route independently constructs

\[
K_\phi=\sum_i\langle\phi|i\rangle a_i U_i,
\qquad \widetilde\rho_p^\phi=K_\phi\rho_pK_\phi^\dagger.
\]

The source helper uses an off-diagonal coherence parameter \(\nu\), with \(|\nu|\le1\). The fully coherent and fully dephased limits are \(|\nu|=1\) and \(\nu=0\), respectively, for the specified phase convention.

## Regression coverage

| Test | Property |
|---|---|
| Pure-source Kraus equivalence | Density-matrix channel equals the single-Kraus construction |
| Dephased source | Cross terms vanish exactly |
| Partial coherence | Unnormalised output is affine in real \(\nu\) |
| Weak value | Exact value for the chosen pre/postselection parameterisation |
| Success probability | Exact finite-\(\epsilon\) expression, not its asymptotic approximation |
| Complete measurement | Two orthogonal outcomes sum to probability one |
| Zero-probability branch | Complete measurement handles a valid zero-probability outcome |
| State validity | Conditional state is normalized and positive semidefinite |
| Periodic-grid momentum | Spectral momentum matrix is Hermitian |
| Kick convention | (U=\exp(i\Delta p X/\hbar)) shifts probe momentum by (+\Delta p) |
| Weak-value kick | Exact finite-grid conditional shift approaches (\operatorname{Re}(A_w)\Delta p_R) |
| Dephased kick control | Conditional shift has the ordinary convex-branch sign |

## Exact-versus-asymptotic correction

For the parameterisation

\[
|\psi\rangle=\frac{\sqrt{1+\epsilon}|L\rangle+
\sqrt{1-\epsilon}|R\rangle}{\sqrt2},
\qquad
|\phi\rangle=\frac{|L\rangle-|R\rangle}{\sqrt2},
\]

in the zero-coupling limit (when the two branches do not distinguish the probe), the exact source-overlap probability is

\[
P_\phi=\frac{(\sqrt{1+\epsilon}-\sqrt{1-\epsilon})^2}{4},
\]

while \(P_\phi\sim\epsilon^2/4\) is only the leading small-\(\epsilon\) asymptotic. Likewise,

\[
(P_R)_w=-\frac{\sqrt{1-\epsilon^2}+1-\epsilon}{2\epsilon}
\sim-\frac1\epsilon.
\]

With branch-distinguishing interaction enabled, the exact success probability is instead computed from Tr(K rho_probe K^dagger) and depends on the probe state and branch unitaries. The audit note and tests preserve this distinction. This avoids treating a source-only overlap or a leading-order expression as an exact interacting probability.

## Validation status and limitations

- Source files and tests have been committed to the default branch.
- A dedicated CI workflow has been added to compile and run the QAG tests on relevant pushes and pull requests.
- An independent numerical smoke calculation of the same periodic-grid equations was run in the analysis environment: with epsilon=0.1 and a right-branch kick of -0.001 (dimensionless units), the coherent postselected shift was +0.0094739473 versus weak-value prediction +0.0094749372 (relative deviation 1.0448e-4); the fully dephased control shifted -0.00045. The source-only overlap probability was 0.00250628144669, matching its analytic expression.
- GitHub Actions QAG regression run **passed**: run ID 46027910833, commit 0d53d8eb51ba1e5b0ec7983e02dab92485a0a04d. The job completed both "Compile QAG module and tests" and "Run QAG regression suite" successfully. Run URL: https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/46027910833.
- The CI result verifies the committed module/test syntax and all 11 QAG unit tests at that commit. It does not validate a gravity-specific Hamiltonian or establish physical antigravity.
- No gravity-specific Hamiltonian, Newtonian force calibration, source-position observable, measurement-apparatus model, or BQG physical projector is implemented here.
- The toy probe rotation is a generic unitary benchmark. It is not a gravitational interaction model.
- This module cannot establish antigravity, a propulsion effect, a violation of conservation laws, or physicalization of BQG.

## Next implementation gate

1. Construct branch unitaries from a declared toy interaction Hamiltonian and compare exact finite evolution with the weak expansion.
2. Add a Stinespring measurement-apparatus register so all outcomes and momentum exchange are represented in one unitary model.
3. Add parameter sweeps over \(\epsilon,\Delta p,\nu\), grid size and packet width, with convergence checks.
4. Keep the BQG embedding blocked until physical relational position operators and a physical projector are constructed.
