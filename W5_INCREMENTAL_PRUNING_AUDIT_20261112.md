# Incremental pruning counterexample — 2026-11-12

## Source control flow

The production helpers use the pattern:

```python
z = dst.get(key, 0j) + scale * amp
if abs(z) > tol:
    dst[key] = z
elif key in dst:
    del dst[key]
```

This is not equivalent to coherent accumulation followed by one final tolerance check.

## Reproduction

Take `tol=1e-11` and repeatedly add the same output-key contribution `a=0.6*tol`. Each individual contribution is below the threshold. If the key is absent, every update computes `z=a`, declines to insert it, and loses it permanently.

| Contributions N | Exact coherent sum | Incremental-pruned result | Error / tol |
|---:|---:|---:|---:|
| 2 | 1.2e-11 | 0 | 1.2 |
| 10 | 6.0e-11 | 0 | 6 |
| 100 | 6.0e-10 | 0 | 60 |
| 1,000 | 6.0e-9 | 0 | 600 |
| 10,000 | 6.0e-8 | 0 | 6,000 |

A sum-then-prune implementation retains the full coherent sum in every row.

## Mathematical consequence

There is no error bound of the form `global_error <= tol` for this incremental algorithm, independent of the number of contributions. For N same-key increments of magnitude `0.6 tol`, the lost coherent amplitude is `0.6 N tol`, which grows linearly with N.

This is more severe than merely dropping one final small residual. It means that per-update pruning can erase a sequence whose eventual coherent sum is arbitrarily larger than the tolerance.

The same pattern appears in the audited W5 `add_dict`, `add_state`, composition `add`, and `add_sparse` helpers (with different tolerances). This establishes a source-level algorithmic risk and a synthetic control counterexample. It does **not** measure how often or how strongly this occurs in a real W5 production run.

## Required correction

1. Accumulate all contributions for a full output key without tolerance pruning.
2. Use compensated or pairwise summation if cancellation/conditioning warrants it.
3. Apply a declared threshold only after the coherent group is complete.
4. Preserve all pre-threshold contributions in a path ledger.
5. If memory forces streaming, use a bounded accumulator that preserves subthreshold partial sums and emits a mathematically justified error bound; do not reset a key to zero after each small increment.
6. Report pathwise and coherent omitted norms separately.

## Status

- Counterexample: locally reproduced.
- Production frequency/magnitude of this effect: not measured.
- Production ZERO_LOSS: not certified.
