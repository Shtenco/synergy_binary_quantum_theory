# BQG C1/C2 continuation audit — 2026-11-27

## Scope and evidence discipline

This is a local finite audit and a candidate numerical-accumulation implementation. It is **not** a proof that C1 or C2 is closed. It does not claim a full production dependency-graph run. Results are classified as finite checks, source-level findings, or open obligations; none is promoted to a theorem without a complete proof.

Repository baseline inspected: `main` at `049f84697523b0e056b06cc7cd8da7baabfc3bb1` (2026-11-26 16:09:19 UTC).

## C1 — K5 depth-6 sector census

An independent local recomputation of the S5 irreducible dimensions using the hook-length formula gives:

| Partition | Irrep dimension |
|---|---:|
| `[5]` | 1 |
| `[4,1]` | 4 |
| `[3,2]` | 5 |
| `[3,1,1]` | 6 |
| `[2,2,1]` | 5 |
| `[2,1,1,1]` | 4 |
| `[1^5]` | 1 |

The weighted multiplicity sum is `3,111,637`, matching the stated K5 depth-6 Hilbert dimension exactly. The unweighted multiplicity sum is `676,346` and must not be confused with the Hilbert dimension. The S4/S3 and S5/S4 coset partitions were independently checked as disjoint exhaustive covers (4 and 5 cosets respectively).

Structural block/column census was complete for the following sectors:

| Sector | Blocks | Columns |
|---|---:|---:|
| `[3,2]` | 2755/2755 | 130903/130903 |
| `[3,1,1]` | 2719/2719 | 153455/153455 |
| `[2,2,1]` | 2749/2749 | 130503/130503 |
| `[2,1,1,1]` | 2712/2712 | 103318/103318 |

**Interpretation:** these are arithmetic-consistency and structural-coverage checks only. They do not calculate the rank of the master map inside each block. In particular, `[5]` vacuum inclusion does not prove nonvacuum injectivity; `[4,1]` still lacks a verified rank witness; and complete census in `[3,2]`, `[3,1,1]`, `[2,2,1]` is not an injectivity proof.

**C1 status: OPEN.** Required next artifact: reproducible exact or certified numerical per-block rank witnesses covering every nonvacuum block in all seven S5 sectors, with explicit operator construction, row/column ordering, rank threshold or exact field, and residual/certificate verification.

## C2 — linearity risk in sparse accumulation

The inspected production sources still contain magnitude-thresholded accumulation patterns. Examples include local `add`/`add_state` helpers that delete a key when an intermediate partial sum falls below a tolerance, plus downstream complete-state filtering. Such state-dependent deletion is not a linear operation and can depend on grouping/order. Exact SU(2) selection rules and fixed, predeclared structural masks are distinct and should remain available; they must not be implemented by amplitude thresholding on intermediate partial sums.

For `tau = 1e-11`, the minimal counterexample `a=1.2 tau`, `b=-0.9 tau` gives `F(a)=1.2 tau`, `F(b)=0`, but `F(a+b)=0`; additivity fails. Reordering cancellation triples can likewise change the result. A Neumaier sparse accumulator candidate is supplied in `scripts/bqg_neumaier_sparse_accumulator.py`; it retains all keys during accumulation and removes only exact zeros at finalization.

### Finite candidate tests

A 50,000-trial local stress run (complex values, 4–32 terms, magnitudes spanning approximately `1e-300` to `1e250`) found:

- Neumaier matched `math.fsum` exactly in 99.74% of trials (naive sum: 95.66%).
- Mean normalized Neumaier error: `2.4033e-19`; maximum: `2.0463e-16`.
- Measured permutation defect: zero in this sample; partition-merge maximum normalized error: `2.0463e-16`.
- Complex linearity defect mean/max: `3.7098e-17` / `1.2631e-15`.
- Adversarial sparse cancellation: 128 keys × 1,000 cancellation triples each (384,000 contributions total). Thresholded accumulation retained 0 keys; Neumaier retained all 128. Maximum absolute error against the target was `1.65436e-24`.

A separate checked-in pytest suite executes seeded single-pass, partition-merge, permutation, adversarial cancellation, and complex-scaling checks against `math.fsum`; it passed locally (4 tests).

**Limit:** compensated floating-point summation is not symbolic exact arithmetic. These tests validate the primitive candidate, not the full W5 production operator.

### Integration obligations

1. Replace every state-dependent magnitude-pruning accumulation site in the production call graph, not merely the final output filter.
2. Audit all direct dictionary writes, merges, and post-operator filters for hidden amplitude pruning.
3. Add basis-vector, pairwise-superposition, complex-scaling, and permutation tests at each public operator boundary.
4. Re-run the complete paired W5 V09/V10 matrix: `Jmax2 ∈ {10,12,14}`, production/shadow, `K ∈ {0,2,4,6,8}` — 30 cells. The currently available paired aggregate had 0/30 reports; no production PASS is available.
5. Preserve the full dependency graph and exact run configuration in the generated evidence bundle.

**C2 status: OPEN.** The accumulator is a research candidate, not an integrated production fix.

## Gate summary

- C1: `OPEN` — finite multiplicity/coset/census consistency checks pass; nonvacuum rank witnesses missing.
- C2: `OPEN` — nonlinear threshold-pruning risk demonstrated; candidate accumulator passes primitive tests; full integration and 30-cell W5 audit missing.
- C3/C4 and all physicalization gates: unchanged and `OPEN`.
- `core_theory_closed_declared=false` and `experimentally_confirmed=false` remain the required declarations.

## Reproduction

Run from repository root:

```bash
python -m pytest -q tests/test_bqg_neumaier_sparse_accumulator.py
```

The tests are deterministic. Passing them does not certify the full operator.
