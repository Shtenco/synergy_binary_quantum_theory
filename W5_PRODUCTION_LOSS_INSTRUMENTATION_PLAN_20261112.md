# W5 production loss instrumentation plan — 2026-11-12

## Scope and evidence boundary

This is a source-level instrumentation map, not a production run result. The files below were fetched from branch `research/lemaitre-bqg-causal-cone-20261110` at source commit `b69b105756dd7d9be093b9bf217ce048bd4a2554` (the source-staging workflow commit). No production loss norm, ZERO_LOSS verdict, or physical closure claim is made by this document.

The source directory contains 176 files. The actual KKV/sine path has a dependency chain; testing only the final script without its imported modules would be incomplete.

## Production call graph

```
k5_peter_weyl_safe_hda_column
  ├─ charged_intertwiner_recoupling_gate
  ├─ peter_weyl_lorentzian_K_block_gate
  └─ peter_weyl_covariant_volume_leg_gate
       ├─ charged_intertwiner_recoupling_gate
       └─ peter_weyl_lorentzian_K_block_gate
  └─ peter_weyl_covariant_K_leg_gate
       ├─ volume_leg_gate
       └─ K_block_gate
  └─ peter_weyl_covariant_composition_gate
  └─ peter_weyl_covariant_K_composition_gate
       └─ peter_weyl_covariant_K_composition_cached_gate
            └─ peter_weyl_covariant_K_sine_composition_gate
                 ├─ peter_weyl_euclidean_sine_ordering_gate
                 └─ peter_weyl_lorentzian_ordered_triple_gate
                      └─ peter_weyl_lorentzian_sine_ordered_triple_gate
```

## Numerical pruning sites to instrument

| Module | Source locations | Existing rule | Required loss accounting |
|---|---:|---|---|
| `k5_peter_weyl_safe_hda_column.py` | 194, 254, 261, 264, 280, 282, 292–293, 331, 337–341, 353 | incremental dictionary threshold 1e-13; intertwiner coefficient 1e-12; output coefficient 1e-12 then 1e-11; cached state 1e-10; final state 1e-8 | record every dropped complex coefficient, keyed by stage and input column; measure L1 and L2 discarded mass; separately record spin-wall omissions at 136 |
| `peter_weyl_lorentzian_K_block_gate.py` | 51, 62, 66, 76 | local projection 1e-12; V state 1e-10; H_E and K state 1e-9 | per-stage discarded norm and propagation through subsequent V/H_E actions |
| `peter_weyl_covariant_volume_leg_gate.py` | 98–151, 175, 239 | charged projection 1e-13; branch merge/output tolerance 1e-11; closed state 1e-10 | per-branch and post-coherent-sum loss; never conflate branch cancellation with threshold discard |
| `peter_weyl_covariant_K_leg_gate.py` | 53–56, 89–134, 164, 181–185, 207–236 | 1e-13, 1e-11, 1e-10 thresholds in projection, charge decomposition and completed H_E | report primitive-branch losses, complete-basis reconstruction leakage, and post-sum pruning separately |
| `peter_weyl_covariant_composition_gate.py` | 40–88 and 130–152 | TOL=1e-11; output projection and state merges | per composition leg, before/after each holonomy close |
| `peter_weyl_covariant_K_composition_gate.py` | 40–225 | inherits composition thresholds and all-J charge projection | preserve exact component comparison against frozen Gauss reference |
| `peter_weyl_covariant_K_composition_cached_gate.py` | 1–256 | claims no additional pruning; caches already-pruned underlying outputs; thresholds 1e-9, 1e-10, 1e-18 in acceptance | prove cache transparency by comparing cached and uncached full state dictionaries and cache key identity |
| `peter_weyl_covariant_K_sine_composition_gate.py` | 41–68 | explicitly preserves historical Gauss K cutoff 1e-9 | report sine path loss separately from inherited cutoff |
| `peter_weyl_lorentzian_ordered_triple_gate.py` | 47–54, 122–164, 201 | TOL=1e-11; exact scalar-channel rank selection after middle K | classify rank-forbidden J=2 removal as exact representation selection, not numerical loss; verify with selection-rule residual |
| `peter_weyl_lorentzian_sine_ordered_triple_gate.py` | 59–154 | cache wrappers around sine path | cached/uncached equality and complete per-leg loss ledger |

## Three distinct truncation classes

1. **Numerical amplitude pruning**: `abs(amplitude) <= tolerance`. Must record the actual discarded complex amplitude, stage, state key, input column, and threshold.
2. **Representation/regulator wall**: for example `so > Jmax2` in `apply_hit_branch`. This is not a floating-point tolerance. Record omitted branch count and exact pre-wall amplitudes, then repeat at a larger wall to test convergence.
3. **Exact selection-rule elimination**: after the second K leg, a rank-(0+1) operator cannot map source J=2 to final J=0. This is a mathematically justified structural zero only if the relevant reduced matrix elements and SU(2) selection rule are independently verified. It must not be added to numerical-loss totals.

## Required audit ledger

For every numerical prune site, record:
- site identifier, source SHA, line/function, threshold;
- number of candidate coefficients and number retained/dropped;
- `sum |a_drop|`, `sqrt(sum |a_drop|^2)`, maximum dropped magnitude;
- pre/post state norm and relative change;
- whether the operation occurs before or after coherent summation;
- downstream propagated error bound (at minimum a conservative triangle-inequality bound);
- deterministic rerun at tolerances `1e-8, 1e-10, 1e-12, 1e-14` where numerically meaningful.

For regulator walls, rerun at successively larger `Jmax2`, preserving the same input, branch ordering, and amplitude tolerance. A stable output norm alone is insufficient; compare the full sparse output vector and its support.

## ZERO_LOSS acceptance

Do not label the production path ZERO_LOSS merely because reported losses are zero after thresholding or because final vectors agree after coherent cancellation. The claim requires:
- all numerical prune sites instrumented and zero discarded coefficients (or a separately approved nonzero error budget);
- all regulator-wall omissions absent or bounded by a convergent wall study;
- cache-on/cache-off component equality;
- branch-resolved and coherently summed diagnostics kept distinct;
- exact selection-rule eliminations proved independently;
- deterministic repeated execution with identical source SHA and parameters.

Until those tests run on the staged source archive, production loss and ZERO_LOSS remain OPEN.
