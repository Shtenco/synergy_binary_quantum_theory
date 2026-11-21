# Erratum to repository contradiction audit — CORE CLOSED semantics

Status: **correction to `REPOSITORY_CONTRADICTION_AUDIT_2026_09_05.md`.**

The earlier audit incorrectly classified `core_theory_closed_declared=true` as a contradiction with open physicalization/generative-extension questions.

That classification is withdrawn.

## Correct reading

The canonical repository explicitly defines three different layers:

```text
STRUCTURAL CANDIDATE CORE   = CLOSED in declared scope
PHYSICALIZATION             = OPEN
EXPERIMENTAL CONFIRMATION   = NO
```

`theory_gates.json` defines core closure over registered gates whose legal statuses are

```text
proved
tested_finite
conditional
```

while stronger uniqueness/universality results are explicitly assigned `closure_role=extension` and do not reopen the declared structural core.

The verifier `scripts/verify_theory_gates.py` makes the scope unambiguous: `core_theory_closed` is a backward-compatible alias for **structural internal candidate only**, and physical projector/history, interacting TT kernel and scale calibration are tracked separately in `physicalization_gates.json`.

Therefore

```text
core_theory_closed_declared = true
```

is internally consistent with

```text
physical_projector_history_closed   = false
physical_TT_kernel_frozen           = false
IR_six_wilson_vector_frozen         = false
common_physical_scale_calibrated    = false
experimentally_confirmed            = false
```

There is no logical contradiction because these booleans answer different questions.

## What remains open without reopening CORE

The following are stronger extensions/generalizations rather than missing registered core gates:

```text
unique attraction from broad microscopic ensembles
arbitrary-graph / arbitrary-habitat HDA theorem
uniform unbounded refinement theorem
broad universality theorem
unique derivation of every blocking prescription
```

They can strengthen or narrow the universality class, but by the repository's declared semantics they are not prerequisites for structural candidate closure.

## Physicalization remains a real blocker for physical claims

Separately, the following remain required before claims about an interacting physical graviton kernel, dark sector, Maxwell dynamics or experiment:

```text
theory-specific physical projector / history
connected interblock physical history
physical Gamma^(2)
physical TT kernel
frozen IR six-Wilson vector
one common physical scale
physical scalar/background cosmology
blind external comparison
```

Thus the corrected status is

```text
CORE CLOSED                          = TRUE
STRUCTURAL CANDIDATE CLOSED          = TRUE
FULL PHYSICALIZATION CLOSED          = FALSE
EXPERIMENTALLY CONFIRMED             = FALSE
```

This erratum supersedes section C3 of the earlier contradiction audit. No numerical or algebraic result is changed.