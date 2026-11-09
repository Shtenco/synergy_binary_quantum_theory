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

## Gate ordering remains frozen

No all-five-node closure, physical Hermitian ordering selection, physical
coefficient fit, full Gram, or joint-kernel calculation is authorized by this
pilot. The previously derived conversion
`K_H=K_comm/(i hbar)` and two-K factor `-1/hbar^2` remains an algebraic
normalization identity; it does not repair this scalar-covariance failure.
