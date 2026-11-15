# Actual K5 master-constraint projector — finite result, 2026-11-15

Status: **PASS_K5_ACTUAL_FINITE_MASTER_PROJECTOR** for the declared finite K5 model and the source snapshot below.

- Workflow run: 43697119767
- Source HEAD: 6516ad1e567b712b63c13b4b68bda45d689f93f3
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
- 720 roundoff-scale raw outputs were below the fixed \(10^{-13}\) column floor; maximum such raw norm \(9.16\times10^{-16}\); aggregate discarded operator Frobenius bound \(9.98\times10^{-15}\).
- The fixed \(10^{-12}\) link-matrix mask omitted zero nonzero entries.
- The finite induced physical Gram matrix was positive semidefinite and had rank 1 after quotienting the kinematical null space.
- Heat-kernel convergence matched \(\|e^{-TM_G}-P_0\|=e^{-T\Delta_M}\) with maximum observed absolute discrepancy below \(7\times10^{-16}\).

## Strict follow-up

A stricter follow-up now includes the Gauss-projection omission itself in the master-operator perturbation bound and requires the fixed link-matrix mask to omit no nonzero entries. That exact implementation is being rerun in workflow 43697222676. The result above is a valid finite-regulator PASS under the earlier finite gate contract, but the stricter error-bound follow-up must pass before the expanded result is promoted into the central physicalization registry.

## Claim boundary

This is an actual finite K5 constraint-family construction, not the earlier synthetic random-matrix selftest. It establishes a one-dimensional zero sector and a positive induced inner product for the declared 32D domain and 60 implemented Euclidean local constraints.

It does **not** establish the full Lorentzian constraint family, arbitrary-graph/refinement compatibility, anomaly-free continuum limit, a continuum rigging map, or the physical Born measure of the complete BQG theory. The physical gate PHYSICAL_PROJECTOR_HISTORY remains open.
