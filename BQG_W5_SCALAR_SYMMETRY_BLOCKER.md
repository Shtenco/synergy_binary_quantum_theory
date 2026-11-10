# W5 Lorentzian pilot: scalar-sector symmetry blocker

Date: 2026-11-10

## Status

The first two executions of the W5 node-0 seed pilot both fail before a W5 KKV column is returned. They are not zero-column results, not closure results, and not physical-kernel results.

- Run 42699139272: failed at 2026-11-09 21:16:31 UTC.
- Run 42699764434: failed at 2026-11-09 21:20:15 UTC.
- Both runs use node 0, Jmax2=10, five W5 seed states, max_states=5, max_depth=0.
- Neither run produced the requested JSON artifact.

The adapter/unit-test workflow 42700029814 completed successfully. Those tests validate finite-dimensional adjoint/index algebra and normalization identities; they do not validate scalar SU(2) covariance of the W5 Peter–Weyl Euclidean operator.

## Exact failure

Both pilot logs terminate in:

```text
RuntimeError: scalar H_E violated canonical J,M preservation
at node 0: input J2=0, canonical M2=0, output=(2, 0)
```

The exception is raised by `restore_scalar_charge_M` in
`scripts/peter_weyl_covariant_K_composition_cached_gate.py`, after
`complete_HE_sine` has assembled the sine-ordered Euclidean action.

The input has source-node J=0, M=0. The assembled output contains a source-node
J=1, M=0 component. A genuinely gauge-scalar local Euclidean Hamiltonian must
preserve the source-node total-J sector (and M). Therefore the current W5
sine-ordered action has not yet passed its scalar-sector covariance gate.

## Quantified follow-up: the rejected sector is tiny, but not yet cleared

A later five-seed pilot (run 42702071445, node 0, Jmax2=10) reached the
same scalar-restoration guard and reported the complete weight distribution:

- expected source sector J=0, M=0: 6.641081550043974;
- rejected J=1, M=0 weight: 5.390512613173637e-20;
- rejected component count: 15.

Thus the rejected-to-expected weight ratio is
8.116919770602643e-21, and the corresponding relative amplitude level is
9.00939496892141e-11. This is numerically tiny, but the existing guard is an
exact nonzero test. It is not valid to call this either a physical symmetry
violation or harmless roundoff until the independent full-magnetic diagnostic
has tested its cutoff and aggregation dependence.

A second run after the initial diagnostic additions (42702605894) produced
the same values. It still produced no W5 columns or JSON pilot artifact.

## Diagnostic implementation and preservation of frozen replay

The independent runner
`scripts/peter_weyl_w5_scalar_covariance_diagnostic.py` evaluates the full
magnetic output before highest-weight restoration, separates T and T-dagger,
and compares H-plus and H-sine for the same W5 seed. It reports per-oriented-
spec sector weights and cutoff-to-cutoff sparse-map differences. Intermediate
state accumulation uses zero pruning; support is reported separately at
1e-13, 1e-11, and 1e-9.

Completed cutoffs are now checkpointed to the requested JSON file as they
finish, so a later timeout does not erase earlier numerical evidence. The
workflow is still computing the 10/12/14 cutoff sweep; no cutoff-convergence
result is claimed yet.

The diagnostic expansion of
`peter_weyl_covariant_K_composition_cached_gate.py` temporarily invalidated
the preregistered numerical-continuation backend hash. That file has now been
restored byte-for-byte from its frozen source commit
`54c5d907338679cdafe1b9c612ae02c1ee684186`; the diagnostic information is
preserved in the runner and CI logs rather than by modifying the frozen
backend. The post-restoration replay is pending.

## Interpretation and non-claims

This is a hard structural blocker. It must not be bypassed by:
- disabling the preservation exception;
- dropping the J=1 component by hand;
- changing the amplitude threshold until the component disappears;
- declaring the resulting KKV columns to be physical.

The failure could originate in the W5-spin-2 representation path, the recoupling/
orientation conventions, the assembled oriented-spec sum, the adjoint-sequence
construction, or cutoff/projection interactions. The present log does not
distinguish these causes. It is not yet justified to identify a unique bug.

The failure does not establish that the underlying continuum operator violates
gauge invariance. It establishes that this implementation and finite-cutoff
pipeline have not demonstrated the required invariance on the tested W5 input.

## Required diagnostic sequence before another full KKV pilot

1. Add a scalar-covariance diagnostic that records the complete output weight
   by source-node J and M before the highest-weight restoration check.
2. Run it on one W5 seed (start with K=0) and split the assembled
   `H_E^sine=(T-T†)/(2i)` contribution into:
   - each oriented node triple/spec;
   - T and T† separately;
   - their signed sum.
3. Check the source-node SU(2) selection rule at each aggregation boundary.
   A non-scalar individual primitive term may contain J=1; the complete
   gauge-invariant oriented sum must be tested before concluding that the
   scalar operator itself fails.
4. Compare the sine-ordered and previously accepted Euclidean conventions on
   exactly the same W5 input, cutoff, and charge basis. Do not compare runs
   with different link spins or cutoffs.
5. Only after the scalar covariance residual is below a declared tolerance,
   restore the highest-weight Wigner–Eckart cache and rerun the five seed
   columns. Then inspect support, column norms, direct-adjoint residual, and
   component leakage.
6. Repeat at Jmax2=10,12,14 with identical amplitude tolerance and complete
   closure conditions. Only then expand to all five source nodes.

## Frozen second-stage replay after backend restoration

The numerical-continuation workflow on the restored backend (run 42704881970)
completed successfully. This is an integrity/replay result, not a completed
second-stage physical calculation. Its inventory reports:

- 419 saved records and 50 second-stage charge checks;
- 50 available direct-H_E records against 352 required;
- zero new H_E records generated in this run;
- zero second-stage C_K actions completed;
- direct reverse-H_L column incomplete;
- full error not certified;
- physicalization incomplete.

The weak reverse-H_L diagnostic has components
`[-3.322481021219022e-17, 0.08463553933259263]`, and the reported
conjugacy difference is `1.9429155614935588e-16`. The small conjugacy
difference does not overcome the incomplete direct-reverse column or certify
the full error. This result is kept separate from the W5 Lorentzian pilot.

## Gate ordering remains frozen

No all-five-node closure, physical Hermitian ordering selection, physical
coefficient fit, full Gram, or joint-kernel calculation is authorized by this
pilot. The previously derived conversion
`K_H=K_comm/(i hbar)` and two-K factor `-1/hbar^2` remains an algebraic
normalization identity; it does not repair this scalar-covariance failure.


## Repeated five-seed pilot: blocker reproduced (2026-11-10)

The latest five-seed pilot, workflow run
[42781449559](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42781449559),
completed with **FAIL** at 2026-11-10 08:28:24 UTC. It again failed inside
`restore_scalar_charge_M`, reached through
`he_sine_reduced -> apply_HE_complete_key -> apply_K_complete_custom ->
C_K_component -> triple_forward`. The exact exception is unchanged:

```text
RuntimeError: scalar H_E violated canonical J,M preservation
at node 0: input J2=0, canonical M2=0, output=(2, 0)
```

The requested `W5_ADJOINT_NODE0_J10_SEED_PILOT.json` was not created; the
artifact upload step explicitly reported that no file existed. Thus this
run returned **zero certified W5 seed columns** and no adjoint-closure
diagnostic. This is a reproduction of the earlier blocker, not a new
independent measurement of its amplitude: the exception occurs before the
runner can return the complete output-sector weights.

The run used node 0, doubled-spin wall Jmax2=10, five seed values
K=(0,2,4,6,8), max_states=5 and max_depth=0. It does not test closure beyond
the seeds.

At the same time, the independent scalar-covariance diagnostic run
[42781449616](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42781449616)
and the refreshed jhalf-to-j1 RG run
[42781449706](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42781449706)
were still running when this note was updated. No result from either is
inferred from the five-seed failure.


## Latest frozen-backend replay confirmation (2026-11-10)

After the status-documentation commits, numerical-continuation run
[42782890490](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42782890490)
completed **SUCCESS** on commit
`5b3f35b83348bdd8b854aee133e9987825963323`. The finite-evidence job
repeated the corrected finite K5 control, replayed the saved weak element
without amplitude regeneration, regenerated the three exact local volume
blocks, and verified the resumable direct-stage inventory.

The final inventory still says: 50 of 352 required direct-H_E records are
available, zero new H_E records were generated, zero second-stage C_K actions
were completed, the direct reverse-H_L column is incomplete, the full error
is not certified, and physicalization is incomplete. Therefore this is a
successful frozen-backend integrity/replay check only; it does not advance
the W5 scalar gate or certify the second-stage physical operator.

## Diagnostic scope correction and exact failure-context probe (2026-11-10)

The partial artifact from scalar diagnostic run
[42781449616](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42781449616)
was recovered from artifact 13746947390. It contains the completed
K=0, node-0, Jmax2=10 slice and all 12 oriented specifications. The 12/14
cutoffs did not complete in that run.

A numerical issue was found in the diagnostic reducer: forbidden-sector weight
was computed as total_weight minus expected_sector_weight. When the forbidden
weight is around 1e-30 and the expected sector is O(1), this subtraction can
round to zero. Commit
[4d017688370ba2da33944a1b20f2cadede8bbc95](https://github.com/Shtenco/synergy_binary_quantum_theory/commit/4d017688370ba2da33944a1b20f2cadede8bbc95)
changes the diagnostic to sum forbidden J,M entries directly and records them
separately. This changes only diagnostic arithmetic.

The recovered one-node K=0 slice gives direct forbidden weights:
- T: 8.684036835017698e-30;
- T-dagger: 7.814664597259229e-30;
- H-plus: 4.529569123612589e-30;
- H-sine: 3.719781592525875e-30.

These are not yet a reproduction of the five-seed failure. An important
correction was established by reading the complete traceback of run
[42781449559](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42781449559)
against the actual call graph. The exception occurs in the FIRST, direct branch
of C_K_component:

C_K_component -> direct_K_covariant -> apply_K_complete_custom
-> apply_HE_complete_state_custom -> restore_scalar_charge_M.

The failing call is therefore the one-endpoint direct K branch with
charged_nodes=(source_v,), BEFORE inverse_complete and before the two-endpoint
K branch. The earlier interpretation that the exception first occurred after
inverse_complete was incorrect and is superseded. The one-endpoint standalone
diagnostic is consequently closer in endpoint structure than previously stated,
but it still does not establish equivalence of the exact canonical input key,
the H_E input state, the full sequence path, or the per-sequence coherent sum.
Its tiny residual must not be used to explain or dismiss the five-seed
5.390512613173637e-20 wrong-sector weight.

The probe was added in commits
[1dbcce60b53f2e64bfc4757a787d785c0d45c15f](https://github.com/Shtenco/synergy_binary_quantum_theory/commit/1dbcce60b53f2e64bfc4757a787d785c0d45c15f)
and
[49a539f24077c860894dab9269356f84a0011272](https://github.com/Shtenco/synergy_binary_quantum_theory/commit/49a539f24077c860894dab9269356f84a0011272).
It captures the exact input, canonical key, charged endpoint, K call kind,
H_E position relative to V, full output J,M distribution and largest components.

Commit
[c20132781af557c529efb41b721edec15ccbce7f](https://github.com/Shtenco/synergy_binary_quantum_theory/commit/c20132781af557c529efb41b721edec15ccbce7f)
extends the probe with passive per-sequence tracing inside the exact
complete_HE_sine loop. It records each oriented specification, adjoint choice,
T sequence, branch count and J,M leakage before and after oriented-spec
coherent accumulation. The physical output accumulator retains the original
addition order and no projection, threshold, operator coefficient or guard is
changed.

This probe is diagnostic only. It does not alter the operator, projection,
tolerance, or acceptance criteria. Its CI result is pending; no physical cause
is assigned yet.


## Exact-key sequence trace and pruning hypothesis (2026-11-10, follow-up)

The per-sequence trace artifact from probe run
[42787724085](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42787724085),
artifact 13750701830, reproduces the exact failing canonical input:
`spins=(4,...,4)`, source label `C(J2=0,M2=0,K12=0,K34=0)`, all other
nodes `G(K=0)`, node 0, K seed 0, Jmax2=10. It records 768 primitive
sequence rows and 12 oriented-specification aggregate rows.

An important qualification: the 12 aggregate rows report zero wrong-sector
weight, but the probe's auxiliary accumulators use the production
`add_state(..., tol=1e-11)`. They therefore do **not** establish exact
cancellation below 1e-11. Primitive sequence rows contain J2=2 (J=1) weights
of order 0.037--0.066, and some sequences also contain J2=4 (J=2) weights.
Individual primitive sequences are not scalar; only their coherent sums are
subject to the scalar selection rule.

The independent unpruned diagnostic was rerun on the same K=0, node-0,
Jmax2=10 seed in workflow
[42787724099](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42787724099),
artifact 13750821413. It explicitly calls
`project_branch_complete_charges(..., tol=0.0)` and accumulates with
`add_state(..., tol=0.0)`. Its full H-sine output has:

- wrong-sector weight: 4.016710228941321e-30;
- total H-sine norm squared: 6.641081686796968.

The exact-key production failure probe from run
[42787724085](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42787724085)
has:

- wrong-sector weight: 5.3905126131736365e-20;
- allowed J2=0,M2=0 weight: 6.641081550043974.

The production wrong-sector weight is about 1.34e10 times the unpruned
diagnostic value. The unpruned H-sine norm squared and production allowed
weight differ by 1.36752995e-7 (about 2.06e-8 relative). Thus the discrepancy
is not only a tiny forbidden tail: the pruning path also changes the
allowed-sector result at a small but measurable level.

Source inspection identifies two independent default pruning layers in the
production path: `project_branch_complete_charges` defaults to
`tol=1e-11`, and `add_state` defaults to `tol=1e-11`. The standalone
diagnostic sets both to zero. The projection routine also retains its
separate internal basis-coefficient cutoff at 1e-13 in both paths. The current
evidence therefore points to an interaction between branch-output pruning,
state-accumulator pruning and summation order; it does not yet identify which
of those is sufficient to produce the 5.39e-20 residual.

A dedicated ablation has been committed:
- script: `scripts/w5_scalar_tolerance_ablation.py`;
- workflow: `.github/workflows/w5-scalar-tolerance-ablation.yml`;
- run: [42795399991](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42795399991).

It compares projection tolerance and accumulator tolerance independently,
flat interleaved accumulation against per-oriented-spec accumulation, reversed
specification order, and a `math.fsum` real/imaginary reference. It uses the
same captured input and operator sequences. It does not modify the operator,
the production guard, or the physical acceptance thresholds. No result is
claimed until this run produces its artifact.

**Correction to the immediately preceding trace interpretation:** the
per-spec zero reported by the probe is a thresholded diagnostic result, not
proof of exact per-spec scalar covariance. The unpruned independent diagnostic
is the relevant evidence for the small residual; the new ablation is required
to determine the production-path cause.
