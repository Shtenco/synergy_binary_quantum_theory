# W5 scalar covariance: pilot status and numerical residual audit

Date: 2026-11-10

## Scope

This note records the current W5 seed-column blocker and the diagnostic pipeline status. It does not alter the Hamiltonian, project away any representation sector, select a physical ordering, or infer a physical kernel.

## Observed v10 pilot residual

GitHub Actions run [42702071445](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42702071445), at W5 node 0 and doubled cutoff Jmax2=10, failed in the scalar H_E preservation guard.

The log reported:
- expected source sector: J2=0, M2=0;
- expected-sector weight: 6.641081550043974;
- wrong-sector weight: 5.390512613173637e-20;
- first wrong component: J2=2, M2=0;
- number of wrong components: 15.

Thus the wrong-sector weight fraction is approximately 8.117e-21, and its relative amplitude scale is approximately 9.009e-11. This is a small but nonzero computed residual. It is not yet classified as floating-point noise: cutoff dependence, per-oriented-spec cancellation, and precision sensitivity remain to be established.

The pilot did not produce a seed-column JSON artifact. No real W5 columns, five-node closure, physical Gram matrix, or joint physical kernel have therefore been certified.

## Diagnostic-run serialization defect

The first three-cutoff scalar diagnostic run [42702605869](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42702605869) reached final JSON serialization after approximately 18 minutes, then failed because a NumPy int64 value was not JSON serializable. This is an output-format defect, not a numerical-operator failure.

The diagnostic script was updated in commit 941123865ba48b8860006b546c69db5ce63d3720 to convert NumPy scalar values through .item() during JSON serialization. A fresh run is active as [42781449616](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42781449616).

## Threshold caveat

The diagnostic call to project_branch_complete_charges(..., tol=0.0) does not remove the lower-level hard-coded coefficient test abs(c) > 1e-13 in the shared projector. Consequently, the current diagnostic is not a literal zero-threshold decomposition. This must be reported and controlled before interpreting tiny residuals as converged symmetry leakage or numerical noise.

Do not change the shared physical operator or weaken its scalarity guard to make the pilot pass. The next diagnostic should compare the current projector against an isolated audit projector with an explicit, recorded coefficient threshold, without changing the frozen backend.

## Independent frozen-replay status

The numerical-continuation workflow [42705095822](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42705095822) completed successfully. This confirms that the frozen second-stage replay passed at that commit; it does not establish W5 scalar covariance or physical closure.

The core-regression workflow [42705095807](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42705095807) had failures in higher-shell and L1/Q4 jobs due to cross-repository artifact download responses returning HTTP 301. The directional-Regge job passed. Those infrastructure failures are separate from the W5 scalar-covariance calculation.

## Acceptance sequence

1. Obtain the scalar diagnostic JSON for cutoffs 10, 12, and 14.
2. Separate the residual by oriented spec, T, T-dagger, and their complete sums.
3. Compare a threshold-controlled independent decomposition and quantify discarded coefficient weight.
4. Repeat at increased numerical precision if the residual is not explained by cancellation and cutoff behavior.
5. Only after scalar covariance is established, rerun all five seed columns and verify norms, support, direct-adjoint agreement, and cutoff stability.
6. Then attempt full five-node closure.
7. Only after closure and physical coefficients/order/regulator are frozen, construct the physical Gram and test the joint kernel.

Until these gates pass, the physical joint-kernel question remains open.
