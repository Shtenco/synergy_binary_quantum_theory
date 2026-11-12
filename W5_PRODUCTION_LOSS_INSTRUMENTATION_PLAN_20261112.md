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


## Additional source findings from the exact W5 call path

The following modules are part of the executable path and must be included in the dependency closure, not treated as optional examples:

- `peter_weyl_zeroaware_volume_migration_experiment.py`
- `peter_weyl_covariant_K_projection_audit_gate.py`
- `peter_weyl_volume_cross_rep_audit.py`
- `peter_weyl_q_block_gauge_audit.py`
- `peter_weyl_covariant_K_sine_composition_gate.py`
- `peter_weyl_lorentzian_sine_ordered_triple_gate.py`
- `bqg_v10_w5_lorentzian_adjoint_closure.py`

The latter is the W5 closure wrapper. It calls the sine-ordered engine while intercepting the raw ordered-triple executor, then constructs forward and direct-adjoint columns across the selected source nodes. A standalone run of the one-triple gate is therefore not a substitute for the V10 closure path.

### Spectral zero-aware cutoff — separate from amplitude pruning

The volume path applies a spectral convention before state amplitudes are generated. In `peter_weyl_zeroaware_volume_migration_experiment.py`, `ZERO_FACTOR=1000.0` and

```
tau = 1000 * machine_epsilon * matrix_dimension * max(1, spectral_scale)
```

Eigenvalues with `abs(ev) <= tau` are assigned zero before `sqrt(abs(ev))`. The same scale-dependent zero-aware convention is independently present in `peter_weyl_volume_cross_rep_audit.py` and `peter_weyl_q_block_gauge_audit.py`.

This is not the same operation as dropping a sparse-state amplitude. It is a spectral/nullspace regularization choice in the definition of the volume block. The first baseline-vs-no-amplitude-pruning comparison must keep this convention fixed on both sides. A separate spectral sensitivity study must vary the factor and compare eigenvalue masks, volume matrices, and downstream complex state vectors. Setting this cutoff to zero is not a valid shortcut: numerical near-zero eigenvalues may then be square-rooted and change the operator itself.

### V10 closure controls that the reference evaluator must expose

At source SHA `54bcf101d1fab06b1c5df1ef0d8aa5fdb9704aac`, `bqg_v10_w5_lorentzian_adjoint_closure.py` has these defaults:

| Control | Default | Consequence |
|---|---:|---|
| `Jmax2` | 10 | finite doubled-spin wall |
| `max_states` | 500 | reachable-support cap |
| `max_depth` | 8 | closure-depth cap |
| `amp_tol` | 1e-12 | forward/adjoint output support pruning |

The wrapper rejects `amp_tol <= 0`, filters the coherently aggregated forward and adjoint maps at `amp_tol`, and marks closure incomplete when the state/depth caps are hit. Its `TOL=0.0` accumulator does not make the full path lossless: lower-level operators still prune, the final column filter is positive, and `apply_hit_branch` still silently rejects upper-wall transitions.

The existing conservative one-step KKV support estimate is doubled spin `Jmax2=20` from a seed doubled spin 4 if all 16 link hits land on one edge. This is a worst-case path bound, not evidence that an amplitude reaches the bound. The V10 default 10 therefore cannot certify one-step carrier completeness. Record every attempted `s_out > Jmax2` branch before rejecting it and repeat with larger walls. For repeated depth, the conservative bound grows as `4 + 16d`; a finite default wall is not a proof of closure.

### Independent reference evaluator — implementation contract

The evaluator must be a shadow execution of the exact pinned source tree, not a reimplementation that silently changes operator definitions.

1. **Pin and verify.** Stage all 176 `scripts/` files plus project metadata and record the repository ref, commit SHA, per-file SHA-256, Python version, and dependency versions. Refuse execution if any source hash differs.
2. **Isolate modes.** Run production and reference modes in separate Python processes/import namespaces. Clear all `lru_cache` instances between modes. Never reuse a cached production result in the reference path.
3. **Reference amplitude accumulation.** In the reference mode, preserve every nonzero binary64 contribution through all numerical amplitude filters and incremental sparse merges. Accumulate by the complete output key first; apply any requested comparison threshold only after the coherent group is complete. Keep exact SU(2) selection rules and negative-spin forbidden transitions as structural exclusions, with explicit reason codes.
4. **Do not conflate the spectral convention.** Hold the zero-aware volume eigenvalue mask fixed for the first comparison. Run the `ZERO_FACTOR` sensitivity ladder separately and report changes to the volume matrix before interpreting downstream state differences.
5. **Instrument the spin wall.** Capture the full branch identity and complex amplitude immediately before the `s_out > Jmax2` rejection. A wall-crossing counter without amplitudes is insufficient.
6. **Instrument closure caps.** Preserve outputs rejected by `amp_tol`, state-cap and depth-cap events. A capped run must return `INCOMPLETE`, never a partial `PASS`. The reference mode must have an explicit no-prune setting; the current V10 positive-only `amp_tol` guard must be bypassed only in the shadow harness, never by editing the pinned production source.
7. **Compare vectors, not only norms.** At each stage compare full complex sparse dictionaries over the union of keys: absolute and relative L2 residual, support symmetric difference, max coefficient residual, pathwise omitted L1/L2, and coherent omitted L2. Propagate signed residual vectors through the remaining operator chain where feasible; use triangle-inequality bounds only as conservative supplements.
8. **Separate cache and operator checks.** Compare cached/uncached full dictionaries under the same mode. Separately compare Lorentzian and sine paths, forward and direct-adjoint columns, branch-resolved and coherent SU(2) commutators, and the complete five-node epsilon sum.

A double-precision no-pruning run is the first control, not an exact-arithmetic proof. If near-cancellation or conditioning makes its result unstable, repeat bounded low-spin subproblems with higher precision or exact SymPy coefficients. Memory/time caps must be reported as incomplete rather than silently converted into pruning.

## Updated evidence status

- Additional source-level findings: confirmed by reading the pinned production files.
- Spectral cutoff sensitivity: not measured.
- Full dependency archive in local runtime: not available.
- Independent no-pruning reference evaluator: design specified, not implemented or executed.
- V10 production closure with instrumented wall/cap events: not run.
- Production ZERO_LOSS: OPEN.


## 2026-11-12 staging completeness correction

A second inspection of the repository root and both staging workflows found that the earlier W5 staging bundles were not complete repository snapshots:

- `.github/workflows/stage-w5-production-source.yml` packages `scripts/` plus a short allowlist of root dependency/configuration files.
- The initial `.github/workflows/w5-production-source-staging.yml` on this research branch packages only `scripts/`; it omits the root-level Python modules, `proof_tools/`, `tools/`, `tests/`, and the root dependency manifest.
- The repository root at the audited revision contains additional tracked source locations (including root Python modules, `proof_tools/`, `tools/`, and `tests/`). Therefore a scripts-only artifact cannot be described as the complete production source tree, even if all 176 scripts are present.

The research-branch staging workflow has now been corrected to package `git archive HEAD` for the exact event commit, rather than copying only `scripts/`. It records the source commit and tree IDs, event ref, Git tree manifest, Python/pip/runtime metadata, pip freeze, per-file SHA-256 manifest, and a SHA256SUMS file. A guard requires `git rev-parse HEAD == GITHUB_SHA`, so a workflow_dispatch or push run cannot silently package a different branch tip. The source ref is taken from `GITHUB_REF` (not the detached checkout's symbolic branch name).

This is a packaging/provenance correction only. It does not alter the physical operators or establish numerical loss values. The corrected workflow's run must complete and its artifact must be downloaded and independently verified before the source-transfer gate can be marked closed. The older scripts-only run is not an acceptable substitute.

Current validation boundary: workflow source was committed and reread from GitHub; no runner has yet executed the corrected packaging steps, so archive contents, hashes, and artifact availability remain unverified. Production ZERO_LOSS remains OPEN.


### Tracked-path completeness hardening

A final packaging review found a subtle completeness issue in the first correction: `git archive` honors `export-ignore` attributes and can therefore omit tracked paths. The workflow was hardened again to enumerate the Git index with `git ls-files -z` and stream those exact paths through tar, preserving symlinks and avoiding shell word-splitting. Before packaging it explicitly fails if any mode-160000 submodule gitlink exists, because a shallow checkout without recursive submodule materialization would otherwise package an empty directory and falsely claim completeness. This repository revision has no `.gitmodules` file; the runtime guard remains as a fail-closed safeguard.

The SHA-256 manifest is generated from the resulting staged tree. The Git tree manifest remains the authoritative path/mode/object inventory, and the archive must be checked against it after download. This is still workflow-source verification only until the runner executes and the artifact is independently inspected.


### 2026-11-12 intermediate spin-wall observability correction

A focused review of the exact V10 execution path found a false-completeness risk that the earlier final-state spin diagnostic did not close:

- `k5_peter_weyl_safe_hda_column.apply_hit_branch()` silently skips `s_out > Jmax2` before constructing the branch.
- V10's `max_spin_reached` is calculated from the final Gauss closure basis after the three-leg KKV maps and scalar projection. It does not describe intermediate link spins inside the KKV legs. An intermediate branch can reach the configured wall, be rejected on its next raising hit, and later legs can otherwise have returned the surviving path to a lower final spin. Therefore the final-basis maximum alone cannot certify that no intermediate wall branch was omitted.
- The V10 wrapper now probes only a primitive input edge already at the configured wall. It calls the original hit routine with a diagnostic wall of `Jmax2+2` to expose the one-step `Jmax2+1` candidate, records its full spin tuple, edge, auxiliary indices and complex amplitude, then returns the original production result computed at the original `Jmax2`. The probe does not add the omitted branch to the production state.
- The closure acceptance predicate now requires `spin_wall_omitted_branch_count == 0`. The result reports pathwise omitted L1/L2, maximum amplitude, per-edge counts and bounded branch samples. These pathwise aggregates are not coherent output norms and must not be interpreted as such.
- Regression tests were added for the fail-closed acceptance condition and for capturing a complex branch immediately before wall rejection.

This is a source-level instrumentation correction, not a measured production result. The test suite has not yet executed on a runner, the W5 staging artifact is still unavailable, and the new diagnostic has not yet been exercised against the real production KKV closure. Numerical amplitude pruning, state/depth cap accounting, and the independent no-pruning shadow evaluator remain open. Production ZERO_LOSS remains OPEN.


### Intermediate-spin support observability extension

The spin probe was extended beyond upper-wall events: it now observes every retained primitive hit branch and records the maximum intermediate link spin before later KKV legs can lower it. It separately counts and reports pathwise amplitudes for intermediate branches above the declared `J=7/2` acceptance ceiling. The closure's `max_spin_reached` now takes the maximum of the final Gauss basis and all observed intermediate hit outputs, and acceptance requires this maximum to remain at or below `3.5`. This closes the related blind spot where a high intermediate spin could disappear from the final scalar basis without ever touching the configured upper regulator wall. The new regression tests reject such a hidden intermediate excursion and verify its branch capture. This remains unexecuted source-level hardening until Actions runs; it is not a measured production result.


### Branch-amplitude representation in the spin probe

The spin-wall probe does not treat the scalar branch prefactor as the entire state amplitude. Each hit branch is factorized as `(spins, local_tensors, scalar_factor)`, so the probe now computes branch L1 and L2 norms as the scalar magnitude multiplied by the corresponding product of local tensor norms, and computes the maximum full product-basis coefficient from the product of local tensor maxima. It records the scalar complex factor, local tensor shapes and a SHA-256 fingerprint for every sampled branch; the first five samples at each diagnostic class also carry the exact flattened real/imaginary tensor payloads. This keeps the aggregate pathwise norms meaningful while bounding JSON output size. The test suite includes a direct factorized-norm check. The probe still records pathwise branch mass, not coherent sums between branches.
