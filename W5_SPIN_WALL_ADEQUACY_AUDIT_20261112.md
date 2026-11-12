# W5 spin-wall adequacy audit — 2026-11-12

## Inputs from source

- W5 seed doubled edge spin: `s_e=4` (`w5_input_key` in V10).
- `apply_hit_branch` changes a doubled edge spin by `s_out=s_in±1` and silently skips `s_out>Jmax2`.
- Conservative link-hit counts from the production operator sequence audit:
  - one `C(V)`: at most 2 hits;
  - one `C(K)`: at most 7 hits;
  - one ordered `C(K) C(K) C(V)` path: at most 16 hits.
- V10 `run_closure` default: `Jmax2=10`; it rejects values below 10 but does not require 20.
- The V10 call forwards this `Jmax2` to the production C(V)/C(K)/KKV composition.

## Derivation

A single link hit changes doubled spin by at most one. If all 16 hits were concentrated on the same edge, the conservative bound after one raw KKV step is

[
s_e^{max}(1)le 4+16=20,
qquad j_e^{max}(1)le10.
]

At repeated raw depth `d`:

[
s_e^{max}(d)le4+16d,
qquad
j_e^{max}(d)le2+8d.
]

These are conservative path bounds, not a claim that a production amplitude actually reaches the bound.

## Finding

The V10 default `Jmax2=10` corresponds to maximum `j=5`, while the conservative one-step support bound is `Jmax2=20` (`j=10`). Therefore the default does not establish that one complete raw KKV action fits inside the spin carrier.

Because the shared `apply_hit_branch` silently continues when `s_out>Jmax2`, any nonzero branch reaching that wall is dropped without its amplitude being recorded. The current maximum-spin diagnostic only reports reached keys and cannot establish that no nonzero key was rejected at the wall.

This is a **proven insufficiency of the default bound as a completeness guarantee**, not proof that a nonzero production amplitude is actually lost in a particular run. Actual wall-crossing amplitudes remain unmeasured.

## Required correction

1. Instrument `apply_hit_branch` before the `continue`, recording the complete input branch, edge, input/output doubled spin, auxiliary path, and complex amplitude.
2. Separate exact forbidden negative-spin transitions from upper-wall truncation:
   - `s_out<0`: exact representation selection rule;
   - `s_out>Jmax2`: numerical carrier/spin-wall omission unless a separate exact bound proves it unreachable.
3. For a one-step conservative no-wall pilot, use `Jmax2>=20`; for depth `d`, use `Jmax2>=4+16d`, subject to memory feasibility.
4. Better: do not use a guessed finite wall as a proof. Retain and report every crossing amplitude, then establish a convergent support bound or exact invariant carrier.
5. Reject `ZERO_LOSS` whenever any nonzero wall event exists or any wall event is uninstrumented.

## Status

- Conservative spin bound: mathematically derived from the audited path hit count.
- Default-wall completeness: **NOT ESTABLISHED**; default 10 is below conservative 20.
- Actual nonzero wall-crossing amplitude: **NOT MEASURED**.
- Production ZERO_LOSS: **NOT CERTIFIED**.
