# Deterministic sparse-amplitude policy — design and verification record (2026-11-11)

> **Status:** reusable numerical primitive and regression suite committed. The K5 producer has not yet been wired to this API; existing K5 HH outputs remain legacy-thresholded until that integration and a fresh computation are completed.

## 1. Artifacts

- Implementation: [`scripts/bqg_sparse_pruning.py`](scripts/bqg_sparse_pruning.py)
- Tests: [`scripts/test_bqg_sparse_pruning.py`](scripts/test_bqg_sparse_pruning.py)
- Existing producer audited separately: [`scripts/k5_peter_weyl_safe_hda_column.py`](scripts/k5_peter_weyl_safe_hda_column.py)
- Related audit: [`K5_PRUNING_AND_NUMERICAL_INTEGRITY_AUDIT_2026-11-11.md`](K5_PRUNING_AND_NUMERICAL_INTEGRITY_AUDIT_2026-11-11.md)

## 2. API and semantics

The new utility deliberately separates two operations:

1. `coherent_aggregate(contributions)`: group all complex contributions by canonical output key, sort deterministically, and use `math.fsum` separately on real and imaginary parts. No threshold is applied.
2. `apply_pruning(aggregated, policy)`: only accepts an already coherently aggregated mapping. It returns the kept mapping and a ledger entry for the exact removed vector.

Supported policies:
- `no_prune`: retain every nonzero representable amplitude.
- `aggregate_then_prune`: prune only after canonical-key aggregation at an explicitly named stage and threshold.

The old producer's legacy thresholds are intentionally not represented as a policy mode in this helper: its historical behavior includes intermediate deletion and cannot honestly be reproduced by a post-aggregation policy. The old producer remains the legacy baseline until a separate, explicit compatibility mode is integrated.

## 3. Ledger contract

Each pruning record contains:
- stage and policy;
- threshold;
- input, kept and discarded key counts;
- maximum discarded amplitude;
- discarded (\ell_1) and (\ell_2) norms;
- SHA-256 digest of sorted canonical-key/real/imaginary records;
- an explicit coherent-aggregation-complete flag.

`verify_ledger_identity` checks the difference between full and pruned maps against the ledger's count, norms and digest.

The helper also implements the telescoping bound
[
\|x_n-\widetilde x_n\|_2
\leq\sum_s
\|A_{n-1}\cdots A_{s+1}\|_2\,\|\delta_s\|_2.
]
The caller must supply valid downstream operator-norm products. The helper does not infer or fabricate them.

## 4. Regression coverage authored

The committed unit suite covers:
1. deterministic aggregation under contribution permutations;
2. cancellation before thresholding;
3. constructive accumulation of individually small contributions;
4. comparison with an independent dense-vector oracle;
5. pruning-ledger identity;
6. no-prune retention of very small nonzero amplitudes;
7. a finite-dimensional telescoping-bound check;
8. refusal to prune when coherent aggregation is not certified complete.

## 5. Verification status

- The implementation and test files were fetched back from GitHub after commit and inspected.
- A local independent Python harness passed the core deterministic-sum, cancellation/constructive-accumulation, ledger-norm, dense-vector, and telescoping-bound checks.
- **The committed `unittest` file itself has not yet been executed against a local checkout**, and the K5 producer has not yet been integrated with the new API. No fresh K5 HH numerical result is claimed.
- The local harness is a core-logic check, not a substitute for running `python -m unittest scripts/test_bqg_sparse_pruning.py` in the repository environment.

## 6. Integration requirements before promoting K5 results

The API must be threaded through:
- primitive path evaluation and local intertwiner-channel extraction;
- `add_dict` and all coherent accumulation sites;
- T-cache and H-cache keys/results;
- state pruning and HH composition;
- commutator construction;
- graph-shell column generation.

The integration must preserve an explicit legacy mode for regression comparison while the new no-prune mode is used for small reference cases and threshold-ladder mode is used for feasible production-size cases. No code path may retain an unlabelled hard-coded threshold.

A no-prune request must not silently fall back to pruning. If memory/time limits prevent completion, the run must fail with a resource-limit status and retain partial diagnostics separately.

## 7. Acceptance gate

Do not call the K5 operator pruning-certified until:
- all hidden cutoffs are removed or surfaced in the policy;
- no-prune aggregation is permutation-invariant;
- tiny dense-oracle action agrees;
- primitive and node adjoints pass bilinear tests;
- raw volume Hermiticity and reconstruction residuals pass;
- threshold ladder stabilizes output and relevant spectral data;
- shell leakage is built coherently as (L^\dagger L);
- all reported results include the policy, threshold ledger, source SHA and run manifest.


## 8. K5 producer integration status

The producer-specific integration contract is maintained in [K5_NO_PRUNE_INTEGRATION_CONTRACT_AND_TEST_PLAN_2026-11-11.md](K5_NO_PRUNE_INTEGRATION_CONTRACT_AND_TEST_PLAN_2026-11-11.md). It identifies all current K5 thresholds and defines the required policy threading and validation gates.

The standalone helper is not yet wired into `k5_peter_weyl_safe_hda_column.py`. In particular, no-prune mode must bypass the local coefficient cutoff, primitive output cutoff, cached T cutoff, state/HH pruning, and top-level H/commutator pruning. The policy identity must be part of every affected cache key. The historical path must remain explicitly reproducible as `legacy`; post-aggregation pruning must not be described as reproducing the historical intermediate-pruning semantics.

Until the producer-specific dense oracle, bilinear adjoint tests, raw volume integrity checks, and threshold/error-propagation gates pass, all current HH values remain legacy-thresholded finite regression values. This document and helper do not claim a production fix.


## 9. Input-integrity hardening — 2026-11-11

The helper now rejects non-finite amplitudes supplied directly to `apply_pruning`, rejects non-finite scales even when the corresponding sparse mapping is empty, validates that ledger-verification `atol` is finite and nonnegative, and rejects overflow in the computed full-minus-pruned difference. Four dedicated negative-path tests were added across these cases (including overflow), bringing the committed regression suite to twelve tests. This is reference-helper hardening only; K5 producer integration remains open.
