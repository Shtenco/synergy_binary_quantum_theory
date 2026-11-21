# BQG physical spectral-history certificate

Status: **fail-closed certification architecture. The actual theory-specific `PHYSICAL_PROJECTOR_HISTORY` remains open until real full-master, spectral-closure and refinement evidence are linked and pass.**

Date: 2026-09-06.

## 1. Why this certificate exists

A finite block-Krylov spectral quotient can be mathematically exact and still fail to be a physical BQG history.

The finite statement is

\[
\mathbb M Q=QJ,
\qquad
V^\dagger f(\mathbb M)V
=E_0^\dagger f(J)E_0
\]

on a certified invariant cyclic seed sector.

The physical statement is stronger. It additionally requires that the operator is the complete theory-specific constrained master on the declared graph-changing habitat, that HDA/target-diffeomorphism requirements are satisfied, and that the candidate physical projector/history survives the declared regulator/refinement limit.

Therefore the repository now separates four closure levels.

---

## 2. Closure level S0 — finite spectral quotient

Input:

```text
certified moments or direct master-Krylov data
```

Gate:

```text
scripts/bqg_spectral_history_graph_gate.py
```

Required result:

```text
finite_spectral_history_closed = true
```

This certifies only the seed-visible finite master functional calculus.

It may expose

```text
candidate_zero_spectral_weight
heat_history
weighted positive spectral zeta
constraint spectral dimension
```

but it is forbidden to emit a physical BQG projector.

### Critical hardening introduced on 2026-09-06

Earlier versions accepted four inline booleans from the moment packet:

```text
domain_complete
master_constraint_certified
quantum_hda_or_explicit_dtarget_certified
source_seed_complete_for_claim
```

and could set

```text
physical_history_closed = true
```

when all four were manually asserted.

That was not sufficiently fail-closed because the spectral gate itself could not verify the claimed external evidence.

The gate has now been hardened:

```text
physical_history_closed   = false
physical_projector_emitted = false
```

unconditionally at this layer.

Even when all inline declarations are `true`, the strongest status available from S0 is

```text
FINITE_SPECTRAL_HISTORY_CLOSED_PHYSICAL_CERTIFICATE_REQUIRED
```

The self-test explicitly attacks this bypass.

---

## 3. Closure level S1 — finite complete full-master spectral history

A separate evidence-linked gate is used:

```text
scripts/bqg_physical_spectral_history_certificate_gate.py
```

It binds the S0 result to the output of

```text
scripts/bqg_constraint_master_assembler_gate.py
```

through the same frozen identifiers:

```text
habitat_hash
domain_hash
convention_hash
master_pencil_hash
```

The full-master artifact must independently prove:

```text
assembler passed
domain complete
quantum HDA closed
matching HDA certificate valid for this exact master
physical finite zero projector emitted by the full-master assembler
all identity hashes present
```

The spectral artifact must independently prove:

```text
finite spectral history closed
same master/habitat/domain/convention hashes
spectral gate did not self-promote physical history
spectral gate did not self-emit P_phys
```

Only then may the certificate say

```text
finite_full_master_spectral_history_certified = true
finite_stage_passed = true
```

This is still not continuum physicalization.

---

## 4. Closure level S2 — physical projector/history under refinement

Repository-level

```text
PHYSICAL_PROJECTOR_HISTORY
```

may close only after an independent refinement/rigging-map certificate with schema

```text
BQG_PHYSICAL_REFINEMENT_HISTORY_CERTIFICATE_V1
```

passes.

The certificate must bind to the finite anchor master and a frozen refinement family, and must establish at least:

\[
\lambda_r(\epsilon)\to0,
\qquad
\frac{\lambda_r(\epsilon)}{\lambda_{r+1}(\epsilon)}\to0
\]

for an asymptotic low sector when no exact zero exists at every finite regulator;

projector convergence under frozen embeddings,

\[
\|P_r(\epsilon')-I P_r(\epsilon)I^\dagger\|\to0;
\]

boundary-history convergence,

\[
B_\epsilon^\dagger P_r(\epsilon)B_\epsilon
\to G_{phys};
\]

and HDA/target residual convergence on that same candidate sector.

The candidate rank must be preregistered from independent structure or survive genuinely held-out refinement levels without changing the selection rule.

Required machine fields are:

```text
schema = BQG_PHYSICAL_REFINEMENT_HISTORY_CERTIFICATE_V1
passed = true
refinement_family_hash != empty
anchor_habitat_hash
anchor_domain_hash
anchor_convention_hash
anchor_master_pencil_hash
low_cluster_scale_separation = true
projector_converged_under_embeddings = true
boundary_history_converged = true
hda_residual_converged = true
rank_rule_preregistered_or_heldout = true
```

Only after S0 + S1 + S2 does the evidence-linked certificate set

```text
physical_projector_history_closed = true
passed = true
physicalization_gate_update_allowed = true
```

This is intentionally stronger than a successful finite calculation.

---

## 5. Closure level S3 — connected source-dressed physical history

Physical projector/history is still not the connected metric/scalar/TT generating functional.

For connected history claims an additional source certificate is required:

```text
BQG_SOURCE_DRESSED_HISTORY_CERTIFICATE_V1
```

with the exact same finite master identity plus

```text
source_operator_set_complete = true
source_dressed_history_converged = true
connected_W_not_Z_used = true
```

The last condition forbids reading disconnected raw amplitudes as connected correlators. The legal chain remains

\[
P_{BQG}
\to Z[J]
\to W[J]=\log Z[J]
\to \Gamma[g]
\to \Gamma^{(2)}.
\]

Only then may the certificate report

```text
connected_source_history_closed = true
```

---

## 6. Hash-link theorem implemented by the certificate gate

Let `A_M` be the full-master artifact and `A_S` the spectral-history artifact.

Define

\[
I_M=(h_{hab},h_{dom},h_{conv},h_M).
\]

S1 requires exact equality

\[
I_S=I_M.
\]

A different master hash, habitat, domain or convention is not approximately accepted. The certificate fails closed.

The refinement certificate may contain later master hashes because the regulator changes, but it must identify the exact finite anchor through

```text
anchor_habitat_hash
anchor_domain_hash
anchor_convention_hash
anchor_master_pencil_hash
```

and separately identify the complete frozen sequence by

```text
refinement_family_hash.
```

This prevents a favorable refinement plot from a different microscopic family being attached to an unrelated finite master.

---

## 7. Explicit negative controls

The certificate self-test rejects all of the following:

### Missing refinement

```text
finite spectral closure = true
full master + HDA linked = true
refinement certificate = missing
```

Expected:

```text
finite_stage_passed = true
passed = false
physical_projector_history_closed = false
```

### Wrong master hash

One characteristically different `master_pencil_hash` in the spectral provenance must make the finite linked stage fail.

### Invalid HDA linkage

A full master whose attached HDA certificate is not valid for that exact master may not authorize S1.

### Forged spectral self-promotion

If a spectral result itself claims

```text
physical_history_closed = true
physical_projector_emitted = true
```

it is rejected by the evidence-linked gate. A lower layer is not allowed to promote itself.

---

## 8. Consequence for the current actual BQG calculation

The current Euclidean production run is only an S0 precursor.

The actual frozen calculation is

\[
Y_i=\mathbb M_E b_i,
\qquad i=0,\ldots,31,
\]

then

\[
\mu_2=Y^\dagger Y,
\qquad
R_1=\mu_2-\mu_1^\dagger\mu_1=B_1^\dagger B_1.
\]

A nonzero `R1` generates the first actual master-Lanczos block `Q1`; it is not a failure.

Subsequent steps are

\[
A_1=Q_1^\dagger\mathbb M_EQ_1,
\]

\[
R_2=\mathbb M_EQ_1-Q_0B_1^\dagger-Q_1A_1.
\]

Continue until a direct residual closes or the declared finite habitat is exhausted.

Even an exact finite Euclidean closure does not reach S1 because the physical full-master layer still requires the common-habitat Lorentzian/HDA construction.

---

## 9. What the spectral graph can genuinely close

The precise answer is:

\[
\boxed{
\text{yes, a derived spectral graph can exactly close the entire finite seed-visible master history}
}
\]

when the next block residual is certified zero.

And, with additional independent evidence,

\[
\boxed{
\text{the same spectral representation can be the computational core of a physical BQG history certificate.}
}
\]

But the logical implication is

\[
\text{spectral closure}
\not\Rightarrow
\text{physical history closure}
\]

without full constraints, HDA/target-diffeomorphism control and refinement/rigging convergence.

This distinction is now enforced by code rather than prose alone.

---

## 10. Current scientific claim boundary

The strongest present claim is:

> BQG now has an executable, operator-derived and fail-closed architecture capable of representing an entire finite seed-visible master-constraint history through a spectral Krylov quotient. The architecture also defines exactly what additional independently hashed evidence would be necessary to promote such a finite result to a physical refinement-compatible history. The actual production BQG physical-history certificate has not yet passed.

This is the status to preserve until the real S0/S1/S2 artifacts exist and agree.
