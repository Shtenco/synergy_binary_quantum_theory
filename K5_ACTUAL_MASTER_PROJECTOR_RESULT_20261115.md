# Actual K5 master-constraint projector — finite result, 2026-11-15

Status: **PASS_K5_ACTUAL_FINITE_MASTER_PROJECTOR** for the declared finite K5 model under the strict numerical-omission contract.

- Strict workflow run: 43698505888
- Strict source HEAD: 74a89b47e0a692509f5df4304e8192eb11bbde86
- Independent successful rerun: 43698498561 (source commit 0306ec29c4671906d7dfa59b8515f8fa67b9f9b8)
- Workflow: .github/workflows/bqg-k5-actual-master-projector.yml
- Implementation: scripts/bqg_k5_actual_master_projector_gate.py
- Domain: 32D fully-active \(j=1/2\) K5 sector
- Codomain: 140D Gauss-invariant K5 Hilbert
- Constraint family: all 60 local graph-changing Euclidean Thiemann-type maps

## Numerical result

| Metric | Zero-sector rank | Positive spectral gap |
|---|---:|---:|
| Identity | 1 | 2712.094489653526 |
| Dense positive metric 1 | 1 | 5190.006355422270 |
| Dense positive metric 2 | 1 | 5866.276172309323 |
| Dense positive metric 3 | 1 | 6225.848254204168 |

Additional checks:

- Maximum pairwise difference between the four zero-sector projectors: \(1.8526\times10^{-15}\).
- Maximum all-constraint residual on the projectors: below \(2.3\times10^{-14}\).
- Independent 4-simplex target constraint residual: \(1.2638\times10^{-13}\).
- Maximum meaningful Gauss-projection loss fraction: \(7.33\times10^{-15}\).
- 720 roundoff-scale raw outputs were below the fixed \(10^{-13}\) column floor; maximum such raw norm \(9.16\times10^{-16}\); aggregate numerical column-mask Frobenius bound \(9.98\times10^{-15}\).
- The strict master-operator perturbation bound uses only this numerical column-mask omission, because the declared constraint maps are \(C_a=P_G T_a\): Gauss projection defines the codomain map and is not an omitted term in \(C_a\). The separately reported Gauss-projection omission ledger is \(1.3227\times10^{-5}\) in the raw-to-projected comparison, while the actual maximum relative projection loss on non-roundoff columns is only \(7.33\times10^{-15}\). This distinction is explicit in code and regression tests.
- Strict relative master-operator error bounds for the four metrics are \(9.30\times10^{-16}, 3.24\times10^{-15}, 2.96\times10^{-15}, 2.88\times10^{-15}\), all below the \(10^{-10}\) acceptance threshold.
- The fixed \(10^{-12}\) link-matrix mask omitted zero nonzero entries.
- The finite induced physical Gram matrix was positive semidefinite and had rank 1 after quotienting the kinematical null space.
- Heat-kernel convergence matched \(\|e^{-TM_G}-P_0\|=e^{-T\Delta_M}\) with maximum observed absolute discrepancy below \(7\times10^{-16}\).

## Strict error-bound interpretation

The first strict attempt (workflow 43697222676) failed because it counted the difference between raw output \(T_a\) and projected output \(P_GT_a\) as a numerical omission from \(C_a\). That is not the declared operator: the actual finite constraint maps are defined as \(C_a=P_GT_a\) with codomain the 140D Gauss-invariant Hilbert space. The corrected contract therefore keeps the Gauss-projection loss as an independent diagnostic and bounds only numerical omissions from the assembled \(C_a\) columns. The correction was tested twice successfully (runs 43698498561 and 43698505888), including the full actual 60-constraint assembly and strict error-bound gate.

## Claim boundary

This is an actual finite K5 constraint-family construction, not the earlier synthetic random-matrix selftest. It establishes a one-dimensional zero sector and a positive induced inner product for the declared 32D domain and 60 implemented Euclidean local constraints.

It does **not** establish the full Lorentzian constraint family, arbitrary-graph/refinement compatibility, anomaly-free continuum limit, a continuum rigging map, or the physical Born measure of the complete BQG theory. The physical gate PHYSICAL_PROJECTOR_HISTORY remains open.
