# W5 SU(2) transitive pruning inventory — 2026-11-15

Status: **PASS_TRANSITIVE_THRESHOLD_INVENTORY**

- Workflow run: 43697302001
- Audited source HEAD: 145fddb5f9fcfa60530bfde43fe156b2df6233e4
- Workflow: .github/workflows/bqg-w5-su2-pruning-coverage.yml
- Inventory implementation: scripts/bqg_w5_su2_pruning_coverage_inventory.py

## Results

| Measure | Result |
|---|---:|
| Modules in direct all-M W5 SU(2) transitive local import closure | 17 |
| Threshold-like comparisons inventoried | 41 |
| Shared PRUNE decision calls inventoried | 26 |
| Unclassified threshold sites | 0 |
| Unclassified PRUNE calls | 0 |
| Plain, unaudited LRU caches | 0 |

The 17 modules were:

- bqg_pruning_audit
- bqg_sparse_linear_accumulator
- charged_intertwiner_recoupling_gate
- k5_peter_weyl_safe_hda_column
- peter_weyl_covariant_K_composition_cached_gate
- peter_weyl_covariant_K_composition_gate
- peter_weyl_covariant_K_leg_gate
- peter_weyl_covariant_K_projection_audit_gate
- peter_weyl_covariant_K_sine_composition_gate
- peter_weyl_covariant_composition_gate
- peter_weyl_covariant_volume_leg_gate
- peter_weyl_euclidean_sine_ordering_gate
- peter_weyl_lorentzian_K_block_gate
- peter_weyl_q_block_gauge_audit
- peter_weyl_volume_cross_rep_audit
- peter_weyl_w5_direct_all_M_su2_gate
- peter_weyl_zeroaware_volume_migration_experiment

## Interpretation

This PASS establishes static coverage of the threshold-like comparisons and pruning decisions in the direct W5 all-M SU(2) operator's transitive local Python import closure. Every identified numerical pruning decision is routed through the shared PRUNE audit/bypass mechanism or is explicitly classified as a non-pruning representation, operator, acceptance, or domain rule.

It therefore supports a genuine instrumented shadow-no-prune run for this exact import closure.

It does **not** establish that production and shadow outputs are numerically equivalent, that the truncation error is negligible, or that the full W5 operator is self-adjoint on its global domain. Those are separate numerical and mathematical gates. Any later source change that expands this import closure requires rerunning the inventory.
