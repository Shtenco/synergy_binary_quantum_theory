# W5 production audit — 2026-11-13

## Scope and provenance

This audit targets the actual V09/V10/V11 W5 Peter–Weyl production path, not a replacement toy operator.

The source snapshot was obtained from the W5 source-staging artifact associated with PR #81:

- source commit: `e7f6b983e2990b1c268d9164f7f2e438f54e9e6e`
- source tree: `60e32745e0fb630964714876775e8cde20e1adab`
- 760 manifest-listed files passed byte-count and SHA-256 verification.
- The staged runner recorded Python 3.12.16, NumPy 2.5.4, SymPy 1.14.0, mpmath 1.3.0. The production audit workflows use Python 3.11 and currently install unpinned NumPy/SymPy; pinning the numerical environment is still required for bit-for-bit reproducibility.

## 1. L1 sparse-accumulator audit

GitHub Actions run `43444754151` completed successfully. The local reproduction of `tests.test_bqg_production_linearity` reports 9 tests, 6 passing and 3 expected failures.

- The strict reference `SparseLinearAccumulator` passes exact linearity/additivity, complex homogeneity, and order-independence tests.
- Legacy production helpers `COMP.add` and `CK.add_state` still prune a running partial sum. The expected failures deliberately preserve this known sub-threshold nonlinearity.
- The strict reference accumulator is **not yet wired** into the active V09/V10/V11 production call graph.

A passing L1 job therefore does not certify production linearity.

## 2. L2 real-production-primitive audit

GitHub Actions run `43444729954` completed successfully in 92.8 seconds: 12 tests, 9 ordinary passes and 3 expected failures.

The safe-domain complex homogeneity and independent-basis additivity checks pass for the real production primitives. The three expected failures characterize threshold-boundary additivity violations in the volume, `H_E`, and full covariant `C(V)` paths.

These are diagnostic expected failures, not acceptance of full linearity. A fixed matrix-level sparsification policy or a proven error-controlled accumulator must replace state-dependent intermediate pruning before a production linearity claim is made.

## 3. Full V10 W5 closure: reproducible blocker

GitHub Actions run `43444229133` executed the five W5 seeds and all five K5 source nodes at `Jmax2=10,12,14`, `max_states=500`, `max_depth=8`, `amp_tol=1e-12`, production pruning, cache on.

**All three cutoff jobs failed at the same guard before producing closure JSON.**

The stack ends in `peter_weyl_covariant_K_composition_cached_gate.restore_scalar_charge_M`, called by the sine-ordered H_E path. The reported transition is:

```text
node 0
input J2=0, canonical M2=0
output J2=2, M2=0
RuntimeError: scalar H_E violated canonical J,M preservation
```

This invalidates the current highest-weight-only scalar Wigner–Eckart reduction for this charged-sector path unless the underlying representation/selection rule is re-derived. Do not suppress the guard to make CI green.

It is not yet a no-go theorem for BQG or proof that the full physical Hamiltonian is non-scalar. The charged label may describe only the local four-leg sector while an external holonomy index carries additional SU(2) charge. The correct total representation and tensor rank must be derived.

## 4. Discriminating diagnostics launched

Two diagnostics were launched on this isolated branch:

- Run `43445958758`: one-node W5 closure in `shadow_no_prune` mode, with all currently instrumented pruning sites bypassed. It was still executing at the latest status check; no artifact/result is available. The ledger itself still reports incomplete threshold-site coverage, so a passing shadow run would not yet prove zero loss.
- Run `43446220696`: paired `production` and `shadow_no_prune` executions of the existing `peter_weyl_w5_scalar_failure_context_probe.py`, restricted to the K=0 seed and source node 0. The `shadow_no_prune` matrix job was cancelled after the runner received a shutdown signal before producing an artifact; this is not a numerical result. The `production` matrix job was still executing at the latest status check. The probe preserves the guard and captures the exact input key, operator stage, per-sequence output J,M weights, support, and largest complex amplitudes.

The production probe is the next artifact to inspect. If the wrong-sector weight appears before pruning or persists in a completed shadow run, the representation/canonicalization or operator construction must be corrected. If it appears only after thresholded accumulation, the exact sequence trace will identify the first asymmetric pruning site. No conclusion is drawn from the cancelled shadow job.

## 5. Required repair sequence

1. Inspect both artifacts from run `43446220696` and the shadow result from `43445958758`.
2. Derive the correct SU(2) action on the full charged state, including external holonomy/matrix indices. Decide whether the operator is rank-0 on the total coupled charge or a nontrivial tensor on the local charge label.
3. Build a finite full-M reference for one fixed W5 seed and node. Compare every magnetic block against the proposed reduced Wigner–Eckart implementation. Keep the existing guard active.
4. Only after the representation is correct, address the threshold-induced nonlinearity. Inventory every pruning site transitively; `lower_level_pruning_coverage_complete=false` remains a hard blocker.
5. Integrate coherent accumulation before a frozen, globally declared output-basis mask. Prove convergence and error bounds; do not replace intermediate pruning with an unbounded dictionary and call it solved.
6. Rerun V10 at Jmax2 10/12/14 with cache on/off and production/shadow comparisons. Require complete closure, no spin-wall omissions, no cap hits, zero unexplained scalar-sector leakage, direct-adjoint residual convergence, and full complex support comparison.
7. Only then rebuild V11 common carrier and the W5 joint-Gram candidate. The current W5 cone certificate is a reduced spectral lower bound, not a measured Lorentzian cancellation or physical kernel.

## Verdict

- L1 strict reference accumulator: **PASS**.
- L2 safe-domain production primitive controls: **PASS**, with three expected threshold-boundary failures.
- Legacy production linearity: **NOT CERTIFIED**.
- V10 full W5 closure at Jmax2 10/12/14: **BLOCKED** by charged-sector J,M preservation assertion.
- Full no-pruning coverage: **INCOMPLETE**.
- Physical Lorentzian W5 operator, joint kernel, and physical causal cone: **OPEN**.

No physical claim is promoted by these finite computational controls.
