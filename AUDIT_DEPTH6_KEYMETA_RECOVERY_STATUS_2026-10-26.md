# BQG depth-6 [3,2] KEYMETA recovery status — 2026-10-26

## Executive status

The [3,2] KEYMETA recovery is **not complete**. Do not infer rank, kernel, or numerical closure from the available artifacts.

Latest audited durable coverage remains:
- Original source snapshot: 20 orbit blocks, 11,169 columns.
- Earlier recovery artifacts: 196 distinct orbit blocks, 23,531 columns.
- Union: 216 / 2,755 orbit blocks (7.84%) and 34,700 / 130,903 columns (26.51%).
- The previous archive audit found no CRC, manifest, or payload SHA-256 integrity errors in the artifacts it inspected.
- No new KEYMETA recovery artifacts were produced by the 2026-10-25/26 attempt.

## Latest execution evidence

Workflow run [40245635062](https://github.com/Shtenco/binary_quantum_theory/actions/runs/40245635062) was triggered by push at 2026-10-25 23:37:54 UTC. The `prepare` job failed in approximately three seconds; the REST jobs record shows no runner ID and no executed steps. A retry (run attempt 2) failed similarly. The accompanying test workflow run [40245635061](https://github.com/Shtenco/binary_quantum_theory/actions/runs/40245635061) also failed immediately without recorded steps.

The run/job logs endpoint did not return a usable log. Therefore, the exact failure cause is **unconfirmed**. The observed pattern is consistent with a runner/billing/Actions availability problem, but this report does not claim that as proven.

## Previous long-run failure

Recovery run [40193824201](https://github.com/Shtenco/binary_quantum_theory/actions/runs/40193824201) ended in failure. Available job records show runner shutdown/cancellation during the numerical replay, and its current run produced only the plan and source snapshot artifacts. The prior durable recovery run [39743112374](https://github.com/Shtenco/binary_quantum_theory/actions/runs/39743112374) still contains the previously audited recovery artifacts.

## Code safeguards prepared

On branch `bqg/depth6-32-resumable-keymeta`:
- Orbit-chunk planning supports a maximum block count as well as a proxy budget.
- A checkpoint manifest is written after each successfully completed orbit block.
- A test simulates failure after the first orbit and checks that the first block's KEYMETA and partial manifest survive.
- Recovery validation verifies prior KEYMETA SHA-256 and rejects duplicate orbit IDs.
- The recovery workflow now accepts a `runner_label` and `max_parallel` input so it can be dispatched to a configured Linux self-hosted runner if hosted runner execution is unavailable.

The checkpoint is only written after the existing frozen engine finishes computing a whole orbit block. It does **not** checkpoint inside the expensive `master_map` computation, so a shutdown during that operation can still lose the in-progress block.

## Required next actions

1. Inspect repository Actions/billing availability and runner status, because the latest jobs failed before a runner started and the logs were unavailable.
2. If hosted runners cannot start, configure an authorized Linux self-hosted runner and manually dispatch the workflow with its exact label. Do not switch to a self-hosted label before that runner is registered and online.
3. Resume from the audited durable artifacts in run `39743112374`; do not recompute already verified blocks.
4. Re-audit every newly produced ZIP, manifest, orbit ID and SHA-256. Only claim complete KEYMETA coverage after the union contains exactly all 2,755 expected blocks and 130,903 columns with no overlap or integrity errors.
5. Treat full numerical closure, rank, kernel, and physical interpretation as separate claims requiring their own checks.

**Scientific claim boundary:** `rank_certified=false`; `numerical_closure_claimed=false`.
