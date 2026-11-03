# Full Lorentzian node1 KA: independent 13-Gram and next genuine micro test (2026-11-03)

**Scientific status:** 13 real full-outgoing columns of ONE auxiliary triple channel calculated and SHA-verified; full corrected Lorentzian HDA / TT are OPEN.

GitHub run [41321196422](https://github.com/Shtenco/binary_quantum_theory/actions/runs/41321196422) **success**: 13/120 real KB input basis keys, Gram rank13, 449 unique full outgoing covariant states, coherent norm 0.03939790876960724, incoherent quadrature norm 0.039092802112961594, squared-norm interference +2.3948038375120337e-05, max normalized overlap 0.5353919106912871 between keys12/13. Input squared weight of 13 computed keys = 0.14883574093624263; the uncomputed 0.8511642590637574 cannot be advertised as output norm/error bound.

Frozen parent KB SHA256: `f486e2250d4a5678a34438467b9e338d730ebf6e1d1ec956d2de7f8c2521b293`.

**New independent micro work**: [run 41619692055](https://github.com/Shtenco/binary_quantum_theory/actions/runs/41619692055) for input keys 105,59,113. Their archived weights are respectively 0.08575943012830507, 0.08249218242433651, 0.048759819051150434. If all succeed, 16 of 120 input keys will be resolved and squared-input coverage rises from 0.14883574093624263 to **0.36584717254003457**. Check full complex amplitude, original source/orientation/initial index and SHA. None of three outputs counted until finished.

**Important engineering correction:** repository `LP.install_sine_cached_stack()` ALREADY has per-process `functools.lru_cache` for the HE, K, inverse, direct, close, CK and CV kernels. Do not reimplement identical in-process caching. Valuable next speedup requires verified **persistent cross-run** basis-element cache keyed by primitive physical operator, normalization/cutoff, full covariant source, charge nodes, target, operator source hash and exact recoupling convention; benchmark and independent recomputation before using it.

**Open blockers:** max branch diagnostic `CK_complete_charge_basis_leakage=1.0` is still not a certified *post-summation* charge gate; need true global leakage. The third step (C_a(K)) on all 120 source keys and all eight auxiliary channels is not complete. The 24 oriented terms for both actual source nodes and Hermitian adjoint not obtained. Do not promote HDA, physical projector, TT, time crystal or Bell as physical.
