# W5 scalar-covariance pruning audit — 2026-11-27

## Executive result

The W5 scalar-covariance blocker is **not resolved**. The historical unpruned diagnostic reported wrong-sector weight (4.016710228941321\times10^{-30}) for the K=0, node-0, Jmax2=10 seed, compared with the production failure weight (5.3905126131736365\times10^{-20}). However, “unpruned” in that comparison means only that two outer tolerances were set to zero. The operator path still contains several internal amplitude cutoffs. The result is therefore evidence that the production branch/state pruning strongly affects the reported residual, not an exact-arithmetic proof of scalar covariance.

No production threshold is changed by this audit. No residual is projected away, and no W5 PASS is claimed.

## Source-level findings

Inspected on branch `research/null-cone-physical-bridge-v2-20261127`:

1. `scripts/peter_weyl_w5_scalar_covariance_diagnostic.py`
   - The standalone diagnostic calls complete-charge projection and accumulation with `tol=0.0`.
   - It decomposes the output into oriented specifications and T/T† contributions.
   - Its support counts at (10^{-13},10^{-11},10^{-9}) are reporting thresholds, not operator thresholds.

2. `scripts/k5_peter_weyl_safe_hda_column.py`
   - `apply_primitive_oriented` drops local intertwiner coefficients below (10^{-12}) (line 254).
   - It drops complete output amplitudes below (10^{-12}) (line 261), then filters the returned map below (10^{-11}) (line 264).
   - `apply_T_cached_state` filters state amplitudes below (10^{-10}) (line 282).
   - `add_dict` drops accumulated entries at or below (10^{-13}) (lines 191–195).
   - `prune_state` has a default (10^{-8}) threshold (lines 292–293).
   Consequently, setting the outer projection/accumulator tolerance to zero does not make this underlying primitive path threshold-free.

3. `scripts/w5_scalar_tolerance_ablation_onepass.py`
   - It correctly evaluates production and zero outer projection tolerances on identical generated branches, and separately evaluates production and zero outer accumulator tolerances.
   - It compares flat, per-oriented-spec, reversed-spec, and `math.fsum` accumulation.
   - It explicitly retains the internal (10^{-13}) basis coefficient cutoff (lines 128–130 and 234–238).
   - Its output is diagnostic only and compares full-map L2 residual, support, allowed-sector weight and forbidden-sector weight.

4. `scripts/peter_weyl_lorentzian_sine_ordered_triple_gate.py`
   - The final scalarity test uses a hard (10^{-8}) wrong-sector fraction, while its underlying composition helpers have their own cutoffs. That gate must not be interpreted as an exact covariance certificate without an independent cutoff-convergence study.

## Historical evidence and workflow status

- Failure-context probe run 42787724085 completed successfully and recorded production wrong-sector weight (5.3905126131736365\times10^{-20}).
- Standalone diagnostic run 42787724099 failed as a workflow, despite the reported artifact/result summarized in `BQG_W5_SCALAR_SYMMETRY_BLOCKER.md); the numerical value is not an independent exact-arithmetic certificate.
- The original tolerance-ablation run 42795399991 is **cancelled**, not passed.
- The one-pass ablation script and workflow are present in the repository. The workflow is configured for `workflow_dispatch`, but no run of the current one-pass workflow has been verified in this audit.
- The current research branch head remains `6f9b3857aecdcd588011c768d50e43ad6230b28b` at the time this note was prepared; main is `a0d3142b1d1d4113987890bc51657e77985dc516`. The branch is 10 commits ahead and 0 behind main.

## Required next experiments

### A. Run the one-pass outer-tolerance ablation

Run `w5-scalar-tolerance-ablation.yml` on this exact branch/commit and preserve its artifact. Read all eight variants, not only the forbidden-sector weight:
- production projection + production accumulation;
- zero projection + production accumulation;
- production projection + zero accumulation;
- zero projection + zero accumulation;
- grouped production;
- grouped zero;
- reversed grouped zero;
- flat `math.fsum` reference.

Require stable allowed-sector norm and full-map L2 agreement as well as reduced forbidden-sector weight. A smaller forbidden-sector residual alone is not a pass.

### B. Remove ambiguity from the internal cutoff ladder in an isolated audit implementation

Do **not** alter the frozen production operator. Create an audit-only variant with each internal cutoff independently exposed as a named parameter:
- local intertwiner coefficient cutoff;
- branch output coefficient cutoff;
- primitive-map output cutoff;
- T-state output cutoff;
- dictionary accumulation cutoff;
- outer complete-charge projection cutoff;
- outer state accumulation cutoff.

Run a one-factor-at-a-time ablation, then a joint zero-cutoff floating-point reference. Record discarded (L^2) amplitude at each layer. Do not use a pseudothreshold or silently restore the production cutoff.

### C. Establish numerical precision and summation sensitivity

For the same fixed generated branch terms, compare:
- standard complex-double summation;
- compensated real/imaginary summation;
- deterministic reversed and permuted term orders;
- extended precision where the backend genuinely supports it.

Do not assume NumPy `longdouble` provides extended complex precision on the runner. Record dtype mantissa bits and platform. If the primitive CG/volume backend remains binary64, higher-precision accumulation alone is not a full high-precision operator evaluation.

### D. Cutoff dependence

Repeat at Jmax2=10, 12, 14 only after preserving identical conventions and seed labels. Compare the forbidden-sector amplitude norm, allowed-sector norm, full-map L2 differences and internal discarded-weight ledger. Do not compare only relative fractions.

### E. Oriented-spec and adjoint decomposition

For every oriented specification, preserve raw T, raw T†, their local Hermitian combination, and the coherent 12-spec sum as separate maps. Report pairwise cancellation norms and the residual after each addition in fixed and reversed order. A single aggregate scalar number is insufficient.

## Acceptance conditions

The W5 scalar-covariance gate remains **OPEN** until:
1. internal cutoff effects are explicitly controlled and their discarded weights reported;
2. full-map outputs are stable under summation order and precision escalation;
3. the forbidden-sector amplitude converges with cutoff without projection;
4. allowed-sector amplitude and support remain stable;
5. the complete five-seed W5 calculation and adjoint-closure checks are rerun on the same declared basis.

Even passing these finite-regulator numerical checks would establish only the declared finite Peter–Weyl operator's covariance/adjoint properties. It would not establish the physical Lorentzian Hamiltonian, continuum HDA, physical kernel, or a BQG spacetime null cone.
