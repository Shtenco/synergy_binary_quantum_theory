# W5 production source cutoff map and runtime staging — 2026-11-12

## Repository / revision

Repository: `Shtenco/synergy_binary_quantum_theory`
Production source ref: `research/lemaitre-bqg-causal-cone-20261110`
Research PR branch: `research/physical-history-causality-20261112`
PR: #81 (draft, open)
Source modules were fetched through the connected GitHub repository interface and inspected by exact line number. This is a static source audit; the files have not yet been staged into the local runtime or executed there.

## 1. Exact source identities

| File | Blob SHA | Lines |
|---|---|---:|
| `scripts/k5_peter_weyl_safe_hda_column.py` | `e7e1e0705c4b833dcee167a2036ef5f8929a4388` | 389 |
| `scripts/peter_weyl_covariant_K_leg_gate.py` | `60b9b3ec878cd41383a0d0ac7862db852c9d9552` | 343 |
| `scripts/peter_weyl_covariant_volume_leg_gate.py` | `7d1fb7fe364e8e05431dcc109a954365c3274b80` | 299 |
| `scripts/peter_weyl_covariant_composition_gate.py` | `1ef229b1dd4b117c890fdce5f7b0025199c5f534` | 311 |
| `scripts/peter_weyl_covariant_K_composition_gate.py` | `4f837adb14e5051df85044196e7bd14a9b670d94` | 243 |
| `scripts/peter_weyl_lorentzian_ordered_triple_gate.py` | `0f9106200ce15499a01d4da9896c33c4a1fec1b3` | 228 |
| `scripts/bqg_v10_w5_lorentzian_adjoint_closure.py` | `54bcf101d1fab06b1c5df1ef0d8aa5fdb9704aac` | 589 |

## 2. Concrete truncation sites

### 2.1 Base K5 Peter–Weyl column

File: `k5_peter_weyl_safe_hda_column.py`

- Lines 130–137: `apply_hit_branch` silently omits `s_out<0` and `s_out>Jmax2` (finite spin wall).
- Lines 191–195: `add_dict` drops accumulated coefficients unless `abs(v)>1e-13`.
- Lines 254–264: local intertwiner coefficients below `1e-12` are omitted; product amplitudes below `1e-12` are omitted; final primitive outputs below `1e-11` are omitted.
- Lines 277–282: cached state output below `1e-10` is omitted.
- Lines 292–293: `prune_state` defaults to `1e-8`.

These are multiple serial truncation boundaries. The source does not preserve each dropped complex amplitude with its complete state key and reason.

### 2.2 Covariant volume leg

File: `peter_weyl_covariant_volume_leg_gate.py`

- Lines 98–151: projection to covariant sectors has local coefficient floors `1e-13`, a configurable `tol=1e-11`, and a final output filter using the same tolerance.
- Lines 154–175: inverse/forward holonomy path; output is filtered at `1e-10`.
- The charged-volume routine also reports a basis leakage residual. A small reconstruction residual is not a complete amplitude ledger: it does not identify every dropped coefficient or distinguish exact zero from numerical omission.

### 2.3 Covariant K leg

File: `peter_weyl_covariant_K_leg_gate.py`

- Lines 53–57: `add_state` filters at `tol=1e-11`.
- Lines 89–134: complete charge projection discards local coefficients below `1e-13`, then filters product amplitudes and output at `1e-11`.
- Lines 181–185: additional local coefficient and projected output filters.
- Lines 230–236: `apply_HE_complete_state` returns only amplitudes above `1e-10`.
- Lines 253 onward: `K=[V,H_E]` is assembled from these already-filtered components.

The all-J charge-basis reconstruction leakage is diagnostic but does not recover the complex amplitudes omitted at the local coefficient filters.

### 2.4 First C(V) composition

File: `peter_weyl_covariant_composition_gate.py`

- Line 45: `TOL=1e-11`.
- Lines 48–56: every `add` may delete the merged output key at or below `TOL`.
- Lines 71–88: covariant projection also filters at `TOL`.
- `C_volume_component` forms direct and `h V h` contributions using these helpers. There is no event-level record of a term removed from either side before their final subtraction.

### 2.5 C(K) composition

File: `peter_weyl_covariant_K_composition_gate.py`

- Line 40: `TOL=1e-11`.
- It delegates aggregation to `COMP.add`, so it inherits the first composition cutoff.
- Lines 54–60: `K=[V,H_E]` combines independently computed `V H_E` and `H_E V` terms.
- Lines 75–94: the outer covariant K-leg composition combines direct and holonomy-dressed K terms. Per-path complex amplitudes are not retained by the returned diagnostic structure.

### 2.6 Ordered K-K-V triple

File: `peter_weyl_lorentzian_ordered_triple_gate.py`

- Lines 47–54: `TOL=1e-11`, used by `add`.
- Lines 99–151: the path is evaluated in auxiliary indices `(i,j,k)` and stages `C(V)`, first `C(K)`, second `C(K)`.
- After the first K action, lines 155–157 retain only doubled source spins `J2 in (0,2)`.

The `J=0,1` restriction is an exact representation-theoretic selection rule for a final rank-0/1 leg to return to source `J=0`. It must be tagged `exact_selection_rule`, not `numerical_threshold`. The same numerical `TOL` used by the sparse accumulator remains a distinct source of possible loss.

### 2.7 V10 closure

File: `bqg_v10_w5_lorentzian_adjoint_closure.py`

- Line 26: `TOL=0.0`; `add_sparse` itself can preserve all nonzero binary64 values when used with the default.
- Lines 53–57: the default sparse accumulator only drops exact zero at `TOL=0`.
- Lines 66–69: scalar-channel filter retains `J2 in (0,2)`; this is the same exact final-rank selection rule.
- Lines 292–298: `run_closure` requires `amp_tol>0`, finite `Jmax2`, `max_states`, and `max_depth`.
- Lines 350–351: forward and direct-adjoint columns are filtered by `abs(z)>amp_tol`.
- Lines 357–361: state/depth caps count omitted support keys but do not retain their full complex amplitudes.

The V10 layer is not lossless merely because its own `TOL=0`: its imported lower-level modules already apply nonzero filters.

## 3. Dependency closure required for a real local run

The seven requested modules are not a self-contained executable bundle. Imports reveal additional required production dependencies:

- `charged_intertwiner_recoupling_gate.py`
- `peter_weyl_lorentzian_K_block_gate.py`
- `peter_weyl_zeroaware_volume_migration_experiment.py`
- `peter_weyl_lorentzian_sine_ordered_triple_gate.py`

The last module is imported by V10 and patches the raw ordered-triple engine; the zero-aware volume experiment is imported and applied by the ordered-triple run. The correct staging unit is therefore the complete repository `scripts/` tree, not just the seven target files. This is why a read-only bundle workflow packages all of `scripts/`.

## 4. Loss accounting semantics

For each truncation point record:
`stage, node, oriented_branch, auxiliary_indices, input_key, output_key, J_in, M_in, J_out, M_out, amplitude_real, amplitude_imag, cutoff, reason`.

Reason codes:
- `exact_selection_rule`
- `numerical_threshold`
- `spin_wall`
- `carrier_projection`
- `state_cap`
- `depth_cap`
- `underflow_or_nonfinite`

For every stage and reason, calculate:
1. pathwise L2 before coherent grouping;
2. coherent L2 after grouping by full output key;
3. number of nonzero omitted contributions;
4. largest and smallest omitted magnitude;
5. per-key complex sum;
6. downstream amplification bound.

Do not record only the residual after accumulation: pathwise events that cancel before thresholding would otherwise disappear from the audit trail. Conversely, do not call the pathwise norm a physical output error; report both norms.

## 5. Exact versus numerical pruning

The middle source-spin projection `J2 in (0,2)` is exact only under the documented final-rank-0/1-to-singlet selection argument and only after the complete preceding K action. It is safe to retain as an explicit algebraic selection rule, with a dedicated regression test that injects J2=4 and proves the final scalar channel is zero.

All `abs(amplitude)>tol` filters are numerical truncations unless separately proved to act only on exact algebraic zeros. A small reconstruction leakage, a successful adjoint check, or a zero coherent residual does not prove zero pathwise loss.

## 6. Local staging workflow

A read-only GitHub Actions workflow was added to the research PR branch:
`.github/workflows/w5-production-source-staging.yml`

Commit:
`625f16eaa50e47feabd4e9cb35b6124d3b71ab90`

Workflow run:
`43203239147`
https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/43203239147

Purpose: package the complete `scripts/` tree with a SHA-256 manifest into a seven-day artifact, solely to transfer the exact branch sources into the local runtime. It does not run scientific validation or modify production source code.

At the time of this report the workflow was still `queued`; no artifact was available to materialize. Local `git clone` also failed because the container has no DNS/network access to github.com. Therefore the source transfer and production execution remain pending, not silently treated as completed.

## 7. Production ZERO_LOSS decision

**NOT CERTIFIED.**

This decision is based on concrete source-level numeric filters, finite spin/state/depth boundaries, and absence of a complete runtime loss ledger. It does not assert that the actual omitted production amplitude is nonzero; the magnitude has not been measured.

A valid certification requires:
- complete source/dependency staging;
- runtime instrumentation before every numerical deletion;
- exact selection-rule events separately proven;
- all spin-wall and state/depth cap amplitudes retained or bounded;
- both pathwise and coherent loss norms;
- downstream operator-norm amplification;
- threshold convergence and independent matrix/adjoint/SU(2) checks;
- proof that the carrier is complete for the claimed domain.

Until all are satisfied, `ZERO_LOSS` remains fail-closed.
