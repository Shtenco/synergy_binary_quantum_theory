# Open physicalization frontier, stronger extensions and experimental tests

This filename is retained for compatibility with historical links.  The repository now distinguishes three logically different layers that must not be collapsed into one Boolean:

```text
structural_candidate_closed              = false
candidate_mathematical_framework         = true
physical_projector_history_closed        = false
physical_TT_kernel_frozen                = false
IR_six_wilson_vector_frozen              = false
common_physical_scale_calibrated         = false
experimentally_confirmed                 = false
```

The 2026-09-29 declaration of `structural_candidate_closed=true` below is a historical snapshot. It is superseded by the 2026-11-13 machine validation: the 55-gate structural ledger still contains the core `PURE_ANCILLA_NODE_REFINEMENT` no-go, so the current verifier reports `structural_candidate_closed=false`. The bilinear refinement route passes its finite gate, but the rejected pure-ancilla construction has not been silently reclassified.

The first statement refers only to the declared exact/finite/conditional **structural candidate** recorded in `theory_gates.json`.  The remaining statements describe the stronger physicalization and experimental frontier recorded in `physicalization_gates.json`.

A GREEN structural CI therefore means that the declared mathematical/computational candidate has reproduced its registered internal arrows.  It does **not** mean that the theory-specific physical inner product, interacting graviton propagator or experimental truth has been established.

See:

- `THEORY_STATUS.md` — current human-readable status;
- `CANONICAL_THEORY_PACKAGE.md` — structural evidence index;
- `physicalization_gates.json` — fail-closed physicalization truth ledger;
- `MASTER_CONSTRAINT_PHYSICAL_PROJECTOR.md` — finite master-projector theorem and its boundary;
- `Q2_RELATIONAL_HISTORY_PROJECTOR.md` — finite relational-history positive control;
- `Q2_RELATIONAL_METRIC_SOURCE_GENERATING_FUNCTIONAL.md` — finite projector/source/Gamma2 positive control;
- `PREDICTIONS_AND_EXPERIMENTAL_TESTS.md` — observational/falsification layer;
- `theory_gates.json` — structural closure roles.

---


## 0. Immediate finite master-kernel frontier — 29 September 2026

Before the continuum physicalization programme below can be promoted from a legal interface to a theory-specific construction, the current corrected $K_5$ depth-6 master-kernel calculation must be finished honestly.

Machine source of truth:

- \`depth6_frontier.json\`
- \`scripts/verify_depth6_frontier.py\`

Current status:

$$
\boxed{\dim\mathcal H_{d=6}=3\,111\,637}
$$

with CLOSED sectors

$$
\boxed{[1^5],\ [5],\ [4,1]}.
$$

Structural support is closed for all remaining mixed irreps, but structural support is not the same as numerical master-kernel closure.

The $S_4$-sign route relevant to $[2,1,1,1]$ has 16 exact $H_0$-null directions. They are genuine local-geometric null directions, not numerical artifacts. They are, however, fully lifted by $H_1$:

$$
\boxed{\operatorname{rank}H_1|_{16}=16/16},
\qquad
\boxed{\sigma_{\min}=1.1304521906426823}.
$$

The remaining giant component has a rank-aware Hall flow

$$
\boxed{130007/130007},
$$

and a planned sparse square minor

$$
\boxed{130007\times130007}
$$

with $47\,543\,521$ expected nonzero scalar entries. Its rank-revealing factorization has **not yet been completed**.

Therefore:

$$
\boxed{[2,1,1,1]_{d=6}\ \text{ACTIVE / NOT YET CLOSED}}
$$

and

$$
\boxed{\ker M^{(d=6)}=\operatorname{span}\{|0\rangle\}\ ?}
$$

remains a target rather than a theorem.

After that sector, the remaining mixed sectors are

$$
[3,2],\quad[3,1,1],\quad[2,2,1].
$$

A future finite depth-6 PASS will still not imply continuum closure; the refinement/rigging/history gates below remain independently necessary.

---


## 1. Open physical gates that DO matter

These are not optional universality embellishments.  They are required before the repository may claim an interacting physical graviton pole or a frozen physical six-Wilson prediction.

### 1.1 Theory-specific physical projector / history

At finite regulator, the exact operator theorem

```text
M_G = C_A^dagger G^AB C_B >= 0
ker(M_G) = intersection_A ker(C_A)
P_phys^(epsilon) = 1_{0}(M_G)
```

is available for every positive-definite constraint metric `G` when the zero sector is isolated.

What is **not** yet established is the candidate-theory refinement/rigging limit built from the actual graph-changing gravitational constraint family, together with the physical boundary/history prescription or a genuinely derived relational clock.

The finite theorem is therefore a solved reference theorem, not a solved continuum gravity projector.

### 1.2 Connected interblock physical history

The finite C8 relational control proves the legal algebraic possibility

```text
combined constraint projector
-> gauge-invariant relational observables
-> Z[J]
-> W[J]
-> connected response
```

without destroying nontrivial conditional evolution.

The current control deliberately uses `R=J` and an externally declared finite clock.  The gravity theory still needs connected multi-block/refinement metric cumulants generated by its own physical amplitude, with vacuum-disconnected pieces removed by the generating functional.

### 1.3 Physical TT kernel

The legal order remains

```text
{C_A}
-> M
-> P_phys / rigging or boundary history
-> Z[J_g]
-> W[J_g]
-> Gamma[g]
-> Gamma^(2)_metric(omega,k)
-> TT projection
-> K_TT(omega,k).
```

The following shortcut is forbidden:

```text
(z - H_constraint)^(-1)  !=  G_physical(omega)
```

without an independently derived physical-history/time construction.

The reduced TT positive control is valuable, but it is not the interacting theory-specific 1PI graviton kernel.

### 1.4 Physical six-Wilson vector

The generic parity-even S4 quartic TT quotient and its six-observable extractor are closed algebraically.  What remains open is the microscopic physical value

```text
c_IR = (c1,c2,c3,c4,c5,c6)_IR
```

read from the derived physical TT pole.  No lower-dimensional `eta/zeta` or single-`Q_tet` ansatz may be substituted after seeing the microscopic result; those remain nested hypotheses only.

### 1.5 One common physical scale

The unit/observable translator is already executable.  The physical calibration is not frozen.

After the dimensionless physical six-vector is derived, either:

1. a microscopic principle derives the common scale, or
2. exactly one preregistered physical datum fixes it.

The same value must then be held fixed for all remaining observables.  Per-observable fitting is forbidden.

---

## 2. Solved reference/positive-control layer

The following results strengthen the route but must not be promoted to the stronger physical claims above:

```text
constraint Feshbach / block-Krylov identities       = PROVED reference
finite master-constraint zero-sector theorem         = PROVED reference
three-node graph-changing HDA finite prerequisite    = FINITE PASS
C8 combined relational projector                     = FINITE POSITIVE CONTROL
projected relational metric-source / Gamma2 tangent  = FINITE POSITIVE CONTROL
```

In particular, the relational positive controls establish that the correct mathematical order is executable; they do not assert that `R=J`, the C8 clock, or their finite source ensemble are the gravitational dynamics of nature.

---

## 3. Stronger non-blocking extensions

These strengthen generality, uniqueness or universality of the already declared structural candidate but are conceptually different from the open physicalization gates above:

```text
an arbitrary-graph / arbitrary-habitat HDA theorem
a uniform theorem over unbounded refinement and collective spin
a universality theorem over broad alternative microscopic rule families
a derivation of every coarse-blocking assumption from a unique microscopic measure
independent alternative discretizations with the same continuum structural class
```

Failure of one such stronger extension can narrow the universality class without retroactively changing a finite theorem into a false theorem.  Conversely, success does not by itself solve the physical projector/history problem.

---

## 4. Experimental tests

The external programme begins only after the physical kernel, six-vector and one-scale rule are frozen.

Examples:

```text
blind gravitational-wave dispersion / birefringence comparison
external likelihood constraints on the frozen six-Wilson TT pole response
held-out observables after one common scale calibration
independent numerical reimplementation of expensive certificates
```

Internal near-hits, selector identities or oracle reconstructions are not external observations.  Spatial dimension near 3, `z` near 1, 4D-like history scaling and the unit-S4 Lambda≈3 control remain internal construction evidence unless a genuinely blind external protocol says otherwise.

---

## 5. Claim discipline

Use the following rule throughout the repository:

```text
finite/exact structural theorem
    != continuum universality theorem
    != physical projector/history
    != interacting physical propagator
    != frozen microscopic observable
    != experimental confirmation.
```

A future contributor may change an `open_physical` gate to `frozen` only by adding the corresponding reproducible evidence and updating the fail-closed machine ledger.  GREEN CI must continue to mean **truthful status reproduction**, not “all questions of physics are solved.”
