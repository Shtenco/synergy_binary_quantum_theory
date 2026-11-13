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


## 2026-11-12 coefficient-ledger and isolated-threshold-comparison increment

### V10 coefficient ledger added (source-level; execution pending)

The V10 closure wrapper now preserves every coefficient it removes at its own final coherent output threshold and every coefficient attached to a new support key rejected by the state/depth cap. Records are split into forward and direct-adjoint channels and include input key, output key, exact binary64 real/imaginary components, absolute magnitude, threshold (for `amp_tol`), input depth and active cap (for closure-cap rejection). The result contains full record arrays plus per-channel count, pathwise L1, pathwise L2 and maximum magnitude.

The zero output threshold is now accepted by V10 for a controlled sensitivity run. It retains every nonzero output coefficient at this final wrapper threshold; it does **not** disable lower-level pruning in imported operator modules.

### Isolated process comparison harness added

`scripts/bqg_v10_w5_amp_tol_shadow_compare.py` runs the same V10 source in two separate Python processes with fresh import namespaces/caches: the requested positive `amp_tol` and `amp_tol=0`. For completed closures it compares the complete complex sparse forward columns over the union of input/output keys, reporting support symmetric difference, absolute/relative L2 residual, maximum coefficient residual and per-input-column diagnostics. If either closure hits a state/depth cap, it returns `INCOMPLETE_CLOSURE_LIMIT_REACHED` and does not compare partial matrices.

The harness deliberately returns `INCOMPLETE_NO_PRUNING_REFERENCE` even after a successful threshold comparison. Lower-level amplitude pruning remains active, so this is a final-output-threshold sensitivity control, not the requested complete no-pruning shadow evaluator and not a ZERO_LOSS proof. Unit tests were added for complex sparse-vector decoding and comparison.

### Lower-level pruning-site inventory revalidated

The exact source files still contain the following active numerical elimination sites requiring independent before/after coherent-sum instrumentation:

- `k5_peter_weyl_safe_hda_column.py`: `add_dict` threshold 1e-13; intertwiner coefficient 1e-12; primitive output 1e-12/1e-11; cached output 1e-10; state threshold default 1e-8.
- `peter_weyl_lorentzian_K_block_gate.py`: local volume 1e-12; V-state 1e-10; H_E/K-state 1e-9.
- `peter_weyl_covariant_volume_leg_gate.py`: charged projection 1e-13; branch projection/merge `tol` default 1e-11; completed state 1e-10.
- `peter_weyl_covariant_K_leg_gate.py`: state merge 1e-11; charged projection 1e-13; completed branch 1e-11; completed C(K) state 1e-10.
- `peter_weyl_covariant_composition_gate.py`: `TOL=1e-11` in state merges and covariant completion.
- `peter_weyl_covariant_K_composition_gate.py`: inherits the composition/branch thresholds; exact scalar-relevance projection must remain classified separately from numerical pruning.
- `peter_weyl_covariant_K_composition_cached_gate.py`: cache wrappers must be tested for full dictionary equality, not treated as an independent numerical cutoff.
- `peter_weyl_covariant_K_sine_composition_gate.py`: historical Gauss H_E cutoff 1e-9.
- `peter_weyl_lorentzian_ordered_triple_gate.py`: local `TOL=1e-11` branch accumulation; exact J=2 scalar-channel elimination is a separate selection-rule class.
- `peter_weyl_lorentzian_sine_ordered_triple_gate.py`: cached wrappers around the same sine path; cache-on/cache-off equality remains required.

The lower-level ledger must retain both (a) branch coefficients removed before coherent summation and (b) coefficients removed after coherent summation. These are not interchangeable: pre-sum pathwise L1/L2 are not the coherent residual norm. Every site needs a stable site identifier, complete state/branch key, input-column identity, exact complex value and threshold. Spectral zero-aware eigenvalue regularization remains a separate operator-definition sensitivity axis, not amplitude pruning.

### Current completion boundary

- V10 final `amp_tol` coefficient ledger: implemented in source, not yet executed.
- V10 state/depth cap coefficient ledger: implemented in source, not yet executed.
- Isolated two-process full complex-vector comparison at final output threshold: implemented, unit tests added, not yet executed.
- Complete lower-level pruning instrumentation: OPEN.
- True no-pruning shadow evaluator: OPEN; the current harness explicitly refuses to claim this.
- Jmax2/amp_tol convergence ladder: OPEN; must run only after staging archive and complete coefficient instrumentation are verified.
- All cap-hit results remain `INCOMPLETE`; no partial matrix can be accepted as PASS.
- Production ZERO_LOSS remains OPEN.


## 2026-11-12 lower-level audit collector implementation increment

A context-local collector now exists in `scripts/bqg_pruning_audit.py`. It is disabled outside an explicit capture context and does not change threshold decisions. At each wired site it records the exact complex discarded aggregate, complete key representation, threshold, site identifier, and (for incremental merges) prior value plus incoming contribution. Exact-zero cancellation/deletion events are retained with zero norm so they are distinguishable from missing instrumentation.

The collector is currently wired into:
- `k5_peter_weyl_safe_hda_column.py`: incremental `add_dict`, primitive intertwiner/product projection, primitive output, cached T-state output, and `prune_state`;
- `peter_weyl_covariant_volume_leg_gate.py`: Gauss/charged projections, branch projection/output, and closed-state output;
- `peter_weyl_covariant_K_leg_gate.py`: incremental `add_state`, charged/Gauss projections, complete-branch projection/output, volume matrix-element/state output, and completed H_E output;
- `peter_weyl_covariant_composition_gate.py`: incremental `add`;
- `peter_weyl_lorentzian_ordered_triple_gate.py`: incremental `add`.

V10 now captures and serializes the lower-level ledger for both complete and cap-incomplete closure results under `lower_level_pruning_ledger`. The original production numerical comparisons remain unchanged. The collector reports only events at instrumented sites; an empty ledger is not a ZERO_LOSS certificate until the coverage manifest is complete.

Regression tests:
- `tests/test_bqg_pruning_audit.py`: context isolation, nested context restoration, threshold accounting, incremental merge accounting, and exact cancellation;
- `tests/test_bqg_pruning_audit_wiring.py`: verifies the four central merge hooks retain the original threshold behavior while emitting records.

The W5 staging workflow now runs these tests, the V10 shadow-comparison tests, and the V10 closure algebra tests before packaging the exact tracked source tree.

### Remaining uninstrumented sites / explicit open boundary

This increment is partial instrumentation, not complete pruning coverage. In particular:
- `peter_weyl_lorentzian_K_block_gate.py` still has its local-volume, V-state, H_E-state and K-state thresholds;
- `peter_weyl_covariant_K_composition_gate.py`, its cached wrapper, and the sine-composition adapters need their own explicit site inventory and cached/uncached event equivalence;
- remaining direct dictionary comprehensions and coefficient filters in the volume/K/composition dependency chain need stable site IDs;
- the spectral zero-aware volume eigenvalue convention remains a separate operator-definition sensitivity axis and must not be counted as amplitude pruning;
- the spin-wall probe remains a separate regulator-omission ledger.

Therefore `complete_lower_level_pruning_coverage=false`, `true_no_pruning_reference=false`, and `Production ZERO_LOSS=OPEN`. The new workflow must pass on the exact final commit and its artifact must be downloaded and verified before these source-level changes are treated as executed evidence.


## 2026-11-12 bounded no-pruning shadow mode — current source state

The earlier sections above are chronological implementation notes. Their statements that the shadow evaluator and lower-level bypass were entirely open have been superseded by the following source changes.

### Implemented in the current branch

- `scripts/bqg_pruning_audit.py` now provides a context-local `bypass()` mode and a single `keep(value, threshold)` policy. Outside the explicit shadow context, `keep` is exactly the original strict `abs(value) > threshold) decision. In shadow mode it retains every exact nonzero coefficient while still excluding exact zeros.
- `record_threshold` and `record_merge` distinguish a genuine exact cancellation (nonzero prior plus nonzero incoming value producing exact zero) from a structural zero or an ordinary subthreshold nonzero value. The latter is recorded as `would_prune_bypassed` in shadow mode.
- The reachable W5 operator path now routes its amplitude cutoffs through this policy: K5 primitive and cached state construction; covariant volume projection and closure; covariant K projection and volume/K state construction; covariant conversion; Lorentzian K-block volume and state pruning; Euclidean sine ordering; ordered-triple accumulation; and the K-composition/cached/sine adapters through those shared functions.
- V10 accepts `pruning_mode=production|shadow_no_prune`. The shadow mode is bounded by the same explicit `Jmax2`, `max_states`, and `max_depth` limits; hitting any cap remains an incomplete result, never a partial PASS.
- `bqg_v10_w5_amp_tol_shadow_compare.py` now launches isolated production and bounded `shadow_no_prune` subprocesses and compares the complete complex sparse forward columns when both closures finish. It also retains both lower-level event ledgers and explicitly refuses to claim `ZERO_LOSS`.
- The independent `peter_weyl_covariant_K_projection_audit_gate.py` has also been connected to the same bypass/ledger policy for its merge and projection stages.
- Regression coverage now includes collector isolation, nested contexts, genuine cancellation, production threshold behavior, merge-site wiring, state pruning, covariant projection, and preservation of subthreshold nonzero coefficients in shadow mode.

### Evidence status and remaining closure conditions

This is now a **bounded shadow evaluator for the explicitly instrumented amplitude-threshold sites**, not yet a certified complete no-pruning reference.

1. The latest staging and PR-gate workflows for the current source head are queued/pending; no successful run or artifact for that exact head has yet been inspected. Source edits and unit-test definitions are not execution evidence.
2. The full call graph still needs a line-by-line coverage manifest proving that every reachable approximate amplitude threshold is routed through `PRUNE.keep`. The independent audit scripts and any optional diagnostic paths must be classified separately from the production V10 call graph.
3. The zero-aware spectral eigenvalue convention is an operator-definition sensitivity axis, not an amplitude threshold; it remains separately controlled and must not be silently bypassed as though it were coefficient pruning.
4. The spin-wall probe and state/depth caps remain independent finite-regulator boundaries. A shadow run that touches either is incomplete for the corresponding support claim.
5. Only after the exact-source staging artifact is verified and both production/shadow closures finish without caps can the full complex support, coherent L1/L2, pathwise L1/L2, adjoint residual, and `Jmax2 × amp_tol` ladder be interpreted.
6. Continuum/refinement convergence and physical regulator removal remain separate obligations beyond this finite W5 calculation.

Current gates therefore remain: `complete_lower_level_pruning_coverage=false`, `certified_true_no_pruning_reference=false`, and `Production ZERO_LOSS=OPEN` until the above evidence is produced. The earlier finite causal-order theorem, reduced analytic cone, and finite interference/decoherence results remain unchanged; this W5 instrumentation work does not upgrade those mathematical results into a full physical BQG proof.


## 2026-11-13 exact-head staging trigger and acceptance guard

This commit intentionally updates the plan file because the W5 staging workflow watches this path. The push must therefore produce a new staging run for this exact commit. The staging job is required to run, before packaging:

1. `test_bqg_pruning_audit.py`;
2. `test_bqg_pruning_audit_wiring.py`;
3. `test_bqg_v10_w5_amp_tol_shadow_compare.py`;
4. `test_bqg_v10_w5_lorentzian_adjoint_closure.py`;
5. `py_compile` over the complete W5 production dependency chain.

The artifact is acceptable only if its embedded `SOURCE_COMMIT.txt` equals the workflow event SHA, the source-tree manifest matches the exact tracked Git tree, every per-file SHA-256 and `SHA256SUMS.txt` entry verifies, and the submodule guard passes. A passing artifact upload without these checks is not sufficient.

The staging run is not itself the production closure experiment. It validates source transfer, regression wiring, and syntax only. Complete transitive threshold coverage, unbounded production/shadow execution, cache equivalence, regulator-wall convergence, and production ZERO_LOSS remain separate gates. No gate is upgraded by this trigger commit alone.


## 2026-11-13 transitive threshold-site static guard

A new regression in `tests/test_bqg_pruning_audit_wiring.py` now walks the repository-local Python import closure starting from the V10 production entry point. At the current source head the closure is 18 modules, including the sine engine, cached K-composition adapter, charged-intertwiner rules, and the zero-aware volume audits.

The guard checks that every threshold-looking `abs(...)` / `np.abs(...)` comparison in that closure is either:
- in a function that routes the actual amplitude decision through `PRUNE.keep`;
- delegated to the instrumented `PW.prune_state`; or
- explicitly classified as a non-amplitude operation.

The explicit non-amplitude classifications are narrow and documented in the test:
- V10 `add_sparse` uses `tol=0` at every call site and removes exact zeros only;
- spin-wall probes test exact nonzero branch existence, while wall omissions are recorded separately;
- source-J and charged-label comparisons are exact SU(2) selection rules;
- zero-aware eigenvalue masks are a spectral operator-definition convention, not sparse-amplitude pruning;
- the audit collector's own threshold comparisons implement the policy being tested.

Local execution of the exact same scanner against the staged source tree passed: 18 local modules discovered, zero unclassified threshold comparisons. The new test is committed at `6017327569c87f2cd77dc02475dea61efe71c968`; CI must still execute on that commit before the guard is accepted as runner evidence.

**Important boundary:** this closes a static source-inventory blind spot, not the full pruning proof. AST scanning cannot prove dynamic dispatch, all runtime branches, cache replay semantics, or that every branch amplitude is propagated correctly. It also does not close the spin-wall, state/depth cap, or spectral sensitivity studies. The acceptance flag `complete_lower_level_pruning_coverage` remains false until runtime branch coverage and ledger consistency are demonstrated.

## 2026-11-13 execution and cache status

A fresh V10 production pilot was attempted with one source node, the minimum five W5 seed states, `Jmax2=10`, `max_depth=0`, and `amp_tol=1e-12`; it did not complete within 240 seconds and produced no result JSON. A standalone sine-ordered triple pilot also did not complete within 120 seconds. These are timeouts, not PASS or numerical FAIL results.

No cache-on/cache-off full-matrix comparison has yet completed. The existing sine engine installs local memoization wrappers and restores monkey-patched functions in `finally`, while V10 clears module-level LRU caches before its audited walk. This is good isolation hygiene, but it is not an equivalence proof. The required next test must execute the same pinned operator columns with caches enabled and disabled, compare complete complex dictionaries and diagnostics, and compare threshold events with explicit accounting for cache-hit replay versus unique operator-column construction.

Production V10 closure, independent no-pruning reference, cache equivalence, `Jmax2 × amp_tol` convergence, and Production ZERO_LOSS remain OPEN.


## 2026-11-13 follow-up: exact staging evidence, fail-closed comparison, and runtime bottleneck

The W5 staging and physical-history/reference workflows both completed successfully for branch head `8f79f75c64d216de387b19bd731eaeb14046860a`. The W5 staging artifact was downloaded and unpacked independently. Its embedded source commit is the PR merge-source commit `e7f6b983e2990b1c268d9164f7f2e438f54e9e6e`; the embedded per-file SHA-256 manifest verified. Five W5-focused test files were rerun directly from that artifact: **59 passed, 8 subtests passed**. This is source/regression evidence only, not a production closure result.

The cache comparison helper now has fail-closed result semantics:
- `PASS_CACHE_EQUIVALENCE` requires exact equality of forward and direct-adjoint complex sparse matrices (zero support symmetric difference and zero coefficient residual) and exact equality of every substantive ledger/diagnostic group;
- `FAIL_CACHE_EQUIVALENCE` is returned for either a matrix mismatch or any ledger mismatch;
- `COMPLETED_CACHE_EQUIVALENCE` is no longer used as a status that could be mistaken for a pass when the comparison itself found a difference.

Regression tests cover successful equality, spin-wall/spectral-ledger mismatches, and a numerical matrix mismatch. The W5 runner now has an explicit context-local cache-off execution mode; it is implemented but has not yet completed a real full-matrix comparison.

A fresh local execution of the staged V10 source was attempted with `Jmax2=10`, one source node `[0]`, the five W5 seed states, `max_depth=0`, and `amp_tol=1e-12`. It did not complete within 260 seconds. A separate local-only timing probe (not committed and not used as scientific output) sampled the hot path: `project_branch_complete_charges` consumed approximately 23.5 seconds of a 25-second interval over 4,074 calls; the projection path made approximately 2.58 million charged-tensor/projection calls in that interval. This identifies a computational scaling bottleneck, not a mathematical failure or a loss bound. A vectorized projection prototype was explored locally, but its first end-to-end timing did not establish a material speedup; it was not promoted into the production branch.

Current execution status:
- Production V10 full matrix: **OPEN**; minimal real run timed out.
- Independent full no-pruning reference: **OPEN**; current shadow shares the production operator implementation and remains partial.
- Cache-on/off full forward and direct-adjoint matrix equality: **OPEN**; harness and fail-closed comparison are implemented, but no full pair completed.
- Cache-on/off complete ordered ledger equality: **OPEN**; same reason.
- Full `Jmax2 × amp_tol` grid and sparse-support convergence: **OPEN**; runner compares forward and direct-adjoint support, but no full grid completed.
- Spin-wall omissions: instrumentation exists; complete production measurement remains **OPEN**.
- State/depth caps: independently classified and reported; cap-free run remains **OPEN**.
- Spectral sensitivity: factor ladder and diagnostics exist; full matrix ladder remains **OPEN**.
- Production ZERO_LOSS: **NOT PROVED / OPEN**.

Do not upgrade any of these scientific gates from CI success, static inventory, unit fixtures, or a partial timeout profile. A future closure claim still requires completed real matrices, a genuinely independent reference, exact cache-on/off ledger/matrix comparison, and complete pathwise/coherent forward-and-adjoint loss accounting.


## 2026-11-13 latest exact-source replay and 300-second production pilot

The latest branch head at this checkpoint is `43970dd2d457fd7516a40b6bc4a4af6a94c7d85e`. Both required GitHub workflows completed successfully for that exact PR head:
- W5 source staging: run [43445461224](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/43445461224), success;
- physical-history/reference gates: run [43445461071](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/43445461071), success.

The W5 artifact was downloaded and independently verified. It contains source commit `dd87ab8dea68cc800b1acbc336d298835d19c26b`, source tree `2a2ecd27044b001646c3d2784bdf017947af4ce1`, and 760 tracked source files. All 767 entries in the outer `SHA256SUMS.txt` verified; all 760 per-file manifest hashes and byte counts verified. The five W5-focused suites were rerun from that exact artifact: **61 passed, 8 subtests passed**.

The cache comparison acceptance logic was tightened after review. It now refuses a full-ledger PASS unless both runs contain every required ledger payload, the observed substantive ledgers match, and both lower-level coverage flags are explicitly `coverage.complete=true`. Missing ledgers produce an INCOMPLETE status; a matrix or observed-ledger mismatch produces FAIL; incomplete coverage cannot be promoted to PASS. The static transitive-import audit also now checks that every sparse threshold/filter decision is connected to an audit record, and that the shared sparse-state filter records each nonzero coefficient before applying its threshold. This is stronger static wiring evidence, not runtime path coverage.

A further production attempt on the exact staged source was run locally with `Jmax2=10`, source node `[0]`, `max_depth=0`, `max_states=500`, `amp_tol=1e-12`, and cache enabled. The process did not complete within 300 seconds and was terminated by the execution timeout; it produced no complete matrix or ledger. This is a compute timeout, not a numerical FAIL and not a bound on omitted loss.

Updated gate status remains:
- Production V10 full matrix: **OPEN — 300-second minimal pilot timed out**.
- Independent full no-pruning reference: **OPEN — shadow uses the same operator implementation and is not an independent derivation**.
- Cache-on/off full forward and direct-adjoint matrix equality: **OPEN — harness exists; no complete real pair**.
- Cache-on/off full ordered-ledger equality: **OPEN — fail-closed coverage gate now prevents false PASS; no complete real pair**.
- Full `Jmax2 × amp_tol` grid and sparse-support convergence: **OPEN — forward and direct-adjoint grid comparisons are implemented; no complete grid**.
- Spin-wall omissions: **OPEN — instrumentation and separate diagnostics exist; no completed unbounded measurement**.
- State/depth caps: **OPEN — independently recorded and a cap-free mode exists; no completed cap-free run**.
- Spectral sensitivity: **OPEN — eigenvalue diagnostics and factor-ladder runner exist; no completed full-matrix ladder**.
- Production ZERO_LOSS: **NOT PROVED / OPEN**.

The next engineering priority is now computational scaling in `project_branch_complete_charges` and its charged-basis contractions. Optimizations must preserve the exact coefficient/event ordering and be checked against the current implementation before promotion. Do not weaken the thresholds, spin wall, or closure criteria to obtain a nominal completion.

## 2026-11-13 rank-four K5 leg contraction optimization

A targeted profile of one complete apply_HE_complete_key pass at Jmax2=1 identified the dominant repeated kernel: apply_axis_np was called 27,438 times and consumed 0.778 s of a 1.436 s profiled run. The old implementation performed two moveaxis operations around tensordot; NumPy reported 55,180 moveaxis and 27,676 tensordot calls in that profile.

The K5 node tensors in this path are rank four. The production helper now dispatches rank-four tensors to four fixed einsum equations (optimize=False) and retains the old generic moveaxis/tensordot fallback for other ranks. No thresholds, spin cutoffs, branch ordering, path enumeration, or physics assumptions were changed.

### Differential evidence

The old contraction was retained as an independent test reference. For four full apply_HE_complete_key comparisons with cold operator caches:

| Input basis index | Jmax2 | Old (s) | New (s) | Speedup | Exact state / ordered ledger |
|---:|---:|---:|---:|---:|---|
| 0 | 1 | 0.5582 | 0.3256 | 1.71x | exact |
| 1 | 1 | 0.5539 | 0.3161 | 1.75x | exact |
| 2 | 1 | 0.5447 | 0.3098 | 1.76x | exact |
| 0 | 2 | 11.7047 | 9.0147 | 1.30x | exact |

For all four comparisons, sparse output keys and complex amplitudes were exactly equal, diagnostics were exactly equal, and every ordered audit event (including its value and key representation) matched exactly. The respective event counts were 1,034; 1,052; 1,047; and 102,305.

Added regressions check every tensor axis on a non-cubic complex rank-four tensor, negative-axis handling, and a full end-to-end H_E operator plus exact ordered-ledger equivalence against the old contraction. The five W5-focused suites now pass 63 tests after this addition; the full instrumented dependency chain also passes py_compile.

### Interpretation and remaining gate

This is a measured local-kernel improvement, not evidence that the production V10 matrix is now tractable. The gain decreases at Jmax2=2, and no Jmax2=10 full matrix has yet completed with the new kernel. The next required experiment is to rerun the same preregistered production pilot without changing amp_tol, max_depth, max_states, or any physical cutoffs, then compare wall time, state counts, support, diagnostics, and complete audit traces. All ten scientific closure gates remain OPEN until their full required evidence exists; ZERO_LOSS remains NOT PROVED.

## 2026-11-13 bounded-cache memory control

The unbounded cache-on production pilot was killed by the 4-GiB cgroup OOM guard. A cache-off attempt stayed memory-light but produced no complete JSON within 260 seconds. This separates two resource problems: unbounded retention causes memory exhaustion, while bypassing all caches makes the operator path too slow for the same bounded wall-time window.

To preserve the numerical operator while limiting retained intermediates, the large state/tensor caches are now bounded LRU caches. The audit-aware wrapper still stores and replays the exact ordered pruning events on hits; eviction only causes the original function body to be recomputed. No threshold, branch rule, spin wall, state/depth cap, spectral convention, or physical coefficient was changed.

The current limits are deliberately concentrated on large objects:
- K-composition and sine-composition state-result wrappers: 2 entries each.
- Complete H_E state cache and local volume-column cache: 2 entries.
- Primitive T branch cache: 16 entries.
- Charged tensor / charged basis / all-M basis / Q123 matrix caches: 32 / 8 / 4 / 2 entries.
- Large K5 intertwiner, spin-matrix, CG-leg, orientation, and volume caches have explicit bounded capacities; the largest volume-matrix caches retain at most 2 entries.

Regressions now verify that an actual complete H_E cache eviction recomputes to the identical complex sparse output and replays the identical ordered event ledger. Separate tests assert the configured capacities and exercise generic bounded-LRU eviction. The five W5-focused suites pass **66 tests**, and the instrumented dependency chain compiles.

### Current full-pilot status

A fresh cache-on pilot is running with the exact previous production settings (nodes=[0], Jmax2=10, max_depth=0, max_states=500, amp_tol=1e-12) and one BLAS/OpenMP thread. At the latest checkpoint it had run for about 4.5 minutes, used about 360 MiB RSS for the Python process, and had not triggered another OOM; it had not yet written a complete result JSON. This is a resource-control improvement only, not a completed matrix. The run must finish and its matrix/ledger be inspected before any gate status changes.

All ten scientific closure gates remain OPEN. Production ZERO_LOSS remains NOT PROVED.

## 2026-11-13 KKV epsilon-orbit repeated C(V) memoization

The bounded-cache production pilot from the previous checkpoint remained CPU-bound for approximately 40 minutes (one worker thread, about 1.8 GiB RSS, no cgroup OOM events) and had not written its result JSON. It was stopped deliberately to avoid leaving a non-terminating resource experiment running; it produced no matrix, ledger, or numerical verdict. This is recorded as an incomplete runtime attempt, not as a numerical failure.

Inspection of the actual KKV contraction loop found a deterministic repeated subproblem: the 24-term epsilon orbit revisits identical C_volume_component requests for the same exact sparse input, source/target leg and auxiliary indices. In each node_kkv_maps call the source covariant state is fixed, and the same volume component is requested repeatedly as the first operator in different epsilon permutations.

The V10 closure runner now wraps this exact component with PRUNE.audit_lru_cache(16). The key contains the complete sorted sparse state (including exact complex amplitudes) and all operator indices. The cached value retains the original output insertion order and leakage diagnostic. On hits, the audit wrapper replays the exact nested pruning-event rows in their original call order. cache_mode=off still bypasses this cache through the existing context-local policy. The original component function is restored in the runner's finally block. A cache hit/miss/capacity report is attached to the lower-level audit metadata.

No operator formula, coefficient, threshold, spin wall, state/depth cap, or physical selection rule changed. The five W5-focused suites pass **66 tests** and the complete instrumented dependency chain passes py_compile on the updated local source snapshot. This is an implementation-level memoization change; the full component-level cache-on/off equality and full ledger comparison still require a real production run.

The previous 40-minute pilot did not complete, so no scientific gate status changes. All ten W5 scientific gates remain OPEN and production ZERO_LOSS remains NOT PROVED. The next run must use the exact production settings and verify that the new C(V) cache actually records hits, reduces wall time, and preserves the complete ordered lower-level ledger.
