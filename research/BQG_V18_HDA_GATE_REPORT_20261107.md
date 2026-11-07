# BQG V18 HDA progress — verified postprocessing + isolated infrastructure check

Date: 2026-11-07. No full C1/C2/C3 proof and no new true Lorentzian SU(2) amplitudes.

## Key numeric postprocessing of genuine prior H_E data

12 true orientation norm-squares from `HE_SINE_12_ORIENTATIONS_TRUE_POSTSUM_CHARGE_20261106.json` sum to `2.0697524883561496`; previously measured coherent full H_E sine norm-squared `1.8234484519418555`; aggregate interference contribution `-0.2463040364142941` (~-11.90017% of orientation quadrature). This proves summed *norm squares* are not full coherent operator amplitudes. Orientation 6 has raw2 `7.718124510711135e-19` despite incoherent2 `0.45531606300320676`. Near-cancellation requires full complex amplitudes and careful arithmetic; do not infer numerical error from their ratio.

For the previously frozen K5 parity theorem with `a=-2/3,c=-32/9`, `||EE w||/|w|=138.2137831058379`, the **necessary but not sufficient** `||Z||` thresholds for exact compensation of odd mixed HDA term `X` on eps=(1/4,1/8,1/16,1/32,1/64) are approximately `351.7506705,705.4094243,1412.7897540,2827.5826975,5657.1849504`; Z=(ES+SE)psi. This depends on parity, same fixed family, bounded Z and independent parity-even D. Real Z amplitude has NOT been measured. Therefore no off-shell full HDA assertion.

W17 reminder: true `E0|W17` rank 17, Galerkin `PE0P` rank 2, so 15 compressed zero vectors not physical kernels. W43 rank 43 E-only, cannot be used as an invariant full H habitat without leakage test.

## CI isolation experiment

- Main v16 first-leg GitHub Actions run `42367951178` failed with steps=null, no artifacts, after removal of three PowerShell statements from bash.
- A NEW minimal runner smoke workflow `.github/workflows/bqg-v18-runner-smoke.yml` contains only echo and python3, with no checkout, pip or BQG computation. Commit `23d78313602a51ea2706597a7f56d1cdff31182f`.
- Minimal run `42371971150` still failed, steps=null, artifacts=0. A repository/account/Actions infrastructure/configuration issue is now strongly indicated, NOT a numerical SU(2) result. Exact cause unknown: check account Actions permissions, billing/quota and runner activity/visibility before further computational restarts.

## Local package

`BQG_V18_AUDIT_AND_VALIDATORS.zip` contains:

- `bqg_v18_numeric.py`: published-data replay and conditional thresholds.
- `bqg_v18_contracts.py`: 32→12 first-leg identity, strict stage-completion metadata gate, leakage vs compression, complete nine-channel HDA algebra for input matrices.
- `bqg_v18_c1_compare.py`: file-backed exact complex sparse coefficient comparison of 32 *completed* Gauss V4 paths with independent four V2 primitives on same SHA; requires genuine data, none supplied this session.
- `bqg_v18_intermediate_overlap.py`: complex first-stage R01 overlap; forbids concluding pruning from orthogonality alone.
- `test_v18.py`: 16/16 local unittest PASS using **synthetic fixtures only for code tests**, not scientific amplitudes.

Critical: metadata flag alone is NOT evidence of provenance; use independently SHA-authenticated, original final V2/V4 artifacts. No physical Lorentzian pair or off-shell HDA is claimed.

Next physical gate: unblock Actions; run 32/32 true stage-3 original V4 paths plus four fully independent V2, compare coefficients 1e-10, then evaluate B on common W5 and nine HDA channels with separately constructed D.
