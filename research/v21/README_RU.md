# BQG V21 — fail-closed genuine Lorentzian C1 stage checkpoints

Date: 2026-11-07. Scope: scientific software / local tests, **no new Lorentzian SU(2) amplitudes**; C1/C2/C3 remain OPEN.

## Corrected V20 flaws

1. Prior stage worker reported `passed=true` for finite output even if volume/outer covariance-basis leakage exceeded production threshold. `v21_checkpoint_guard.py` now refuses this status for leak >=1e-8. **Charge leakage is separate:** individual `complete_charge_basis_leakage=1` cannot be interpreted as failure of coherently summed full HE.
2. Previously stages 1->2 were incorrectly forced to match all `aux` and `order` indices. Exact first operators can coincide for different path indices; `v21_stage_logic.py` validates actual first-leg signature, source, graph and backend. 32 complete auxiliary paths for one pair require only 12 distinct first operations. Stage 2->3 keeps strict `aux`/`order` checks and scalar filter.
3. `pl_v21_pair_collector.py` requires all 32 original stage-3 outputs, physical basis leakage guard, SHA-verified covariant NPZ, consistent operator signatures, complete 8-channel complex sum per each of four primitive orderings, independent original Gauss projector and direct complex comparison against FOUR independently obtained V2 primitive NPZ. V2 reference metadata must also be SHA-checked and match source/node/orientation/mode/cutoff; NPZ SHA alone is not evidence of independence.

## Costs

Naive three-stage operator calls per physical pair: 32*3=96. Frozen exact first-stage signatures: 12 distinct, so 12+32+32=76 calls without any higher-stage reuse (20 fewer, a 20.833% decrease in call *count*, not runtime). This **does not mean 32 physical paths have been computed**. Real completed paths from this local run: **0/32**.

## Scope and reproducibility

Local standalone `BQG_V21_LOCAL_C1_REPAIR.zip` includes 18 unchanged V20 unit tests, 26 V21 unit tests, and SHA manifest. They pass locally (44/44), using mock/small test vectors only for programming contracts; NO mocked state counts as physical evidence. The full original private repository cannot be cloned from this local offline execution container. GitHub branch contains `scripts/v21_checkpoint_guard.py`, `scripts/v21_stage_logic.py`, amended `scripts/pl_v20_stage_worker.py`, `scripts/pl_v21_pair_collector.py` and `tests/test_v21_c1_checkpoint_repair.py` (18 focused regression tests). GitHub tests have not been shown to pass Actions.

From full repository checkout:

```bash
python -m unittest discover -s tests -p 'test_v21_c1_checkpoint_repair.py' -v
python scripts/pl_v20_stage_worker.py --help
python scripts/pl_v21_pair_collector.py --help
```

C1 collector manifest requires `source_node, omit, cycle, source_commit, Jmax2=7, volume_backend_sha256, input_state_sha256, paths` with one entry for each of 32 paths, each containing `stage1,stage2,stage3` checkpoint prefixes, and `independent_v2` entries `fp,fq,ap,aq` with `npz,sha256,metadata_json,metadata_sha256`. All files must be from the SAME frozen original production operator definition, with independent V2 execution and consistent normalization and basis. A self-declared manifest cannot by itself prove independent provenance; inspect original job/commit artifacts too.

**Research state remains** `C1=false`, `C2=false`, `C3=false`, `HDA_proved=false`, `physical_TT=false`.

## Infrastructure

2026-11-07 new minimal no-dependency GitHub Actions runner failed with zero job steps. Do not interpret this as a physical Hamiltonian result. Resolve runner admission separately before scheduling heavy 32-path execution.
