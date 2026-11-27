# W5 scalar-covariance failure trace — 2026-11-27

Recovered from GitHub Actions run #11 (ID `42787724099`) and the paired first-failure context probe run #2 (ID `42787724085`).

## Reproducible failing path

The W5 adjoint-closure pilot run #14 (ID `42787724159`) failed in:

`bqg_v10_w5_lorentzian_adjoint_closure.py`
→ `peter_weyl_lorentzian_sine_ordered_triple_gate.py:he_sine_reduced`
→ `peter_weyl_covariant_K_composition_cached_gate.py:restore_scalar_charge_M`.

Exact exception:

`scalar H_E violated canonical J,M preservation at node 0: input J2=0, canonical M2=0, output=(2, 0)`.

The failure context pins:
- source node 0, target node 3, (i=j=0), Jmax2=10;
- first (H_E) call before volume;
- scalar seed (J2=0,M2=0);
- 159 output states;
- final allowed-sector weight 6.641081550043974;
- final wrong-sector (J2=2,M2=0) weight 5.3905126131736365e-20, corresponding to amplitude norm about 2.32175e-10.

## Cancellation profile

The observer captured individual signed sequence outputs before complete oriented-spec aggregation. Many individual sequence rows contain wrong-sector weights around (3.7\times10^{-2}) to (6.7\times10^{-2}). These contributions cancel across the signed sequence sum; after all 12 oriented specs, the residual wrong-sector weight is approximately (5.39\times10^{-20}).

This is a high-cancellation path. It is exactly the type of path where intermediate state-dependent pruning and summation order need to be tested at the actual operator boundary, rather than inferred from a small random accumulator test.

The failure context also reports:
- `operator_changed=false`;
- `projection_changed=false`;
- `guard_disabled=false`;
- max volume leakage (5.38\times10^{-16});
- max basis leakage 1.0.

The W5 scalar-covariance diagnostic run #11 separately failed while serializing a NumPy `int64` value to JSON. That is a reporting bug in the diagnostic job, not a physical covariance PASS or FAIL.

## Interpretation and next experiment

1. Re-run this exact pinned seed and capture the full per-key amplitude map at each oriented-spec boundary.
2. Compare (a) production accumulator, (b) Neumaier accumulation without intermediate magnitude pruning, and (c) high-precision or exact-CG reference where feasible.
3. Preserve the same branch generation, projection masks, operator sequence, and basis cutoffs; vary only the accumulator/pruning layer.
4. Check scalar (J,M) leakage after each oriented-spec sum and after the full 12-spec sum.
5. Do not treat the current (2.32\times10^{-10}) residual as proof of a physical scalar-sector violation; it is a finite numerical implementation failure that requires root-cause isolation.

**C2 remains OPEN.** This trace provides a concrete high-cancellation regression target, not a completed linearity certificate.
