# Physical sine-ordering two-node Euclidean × route HDA result

Status: **preregistered finite PASS**.

The protocol and thresholds were frozen in `PETER_WEYL_TWO_NODE_SINE_HDA_PREREGISTRATION.md` before the calculation. No threshold, channel normalization, subtraction, route coupling or phase was changed after observing the result.

## Provenance

GitHub Actions:

```text
workflow: two-node-sine-hda
run:      31855735615
head:     c02c1892cd9d18b9a62413d6c06abf1211e6b7f6
artifact: 9239198320
digest:   sha256:21e2da508fd583d9007a5bd400d074e8cee39990656e6c75e5968d2601323526
```

Both the independent `H_E^sine/K_sine` prerequisite and the preregistered two-node gate completed successfully.

## Physical Euclidean ordering

The calculation uses

\[
\boxed{H_E^{sine}=(T-T^\dagger)/(2i)}
\]

throughout, with

```text
Jmax=5/2
L=48
carrier=8
all ten links j=1/2
all five K=0
zero-aware volume convention
state prune tolerance=1e-8.
```

The one-node actions are

```text
support(H0)=37
support(H1)=37
||H0||=2.1712581763270546
||H1||=2.171258176327055.
```

The genuine two-node sine Euclidean commutator is nonzero:

```text
support([H0,H1])=514
||[H0,H1]||=2.8794538147049544.
```

## Regulator sequence

| epsilon | route-only | cross / D | pure GG / D | joint / D |
|---:|---:|---:|---:|---:|
| 1/4  | 1.3221741194e-5 | 0.3257885267 | 0.03600350128 | 0.3277718969 |
| 1/8  | 6.6115292991e-6 | 0.1614789446 | 0.008898535351 | 0.1617239423 |
| 1/16 | 3.3058470645e-6 | 0.08038526165 | 0.002211784339 | 0.08041568437 |
| 1/32 | 1.6529339109e-6 | 0.04010403764 | 0.0005513364039 | 0.04010782728 |
| 1/64 | 8.2646874425e-7 | 0.02002986590 | 0.0001376326747 | 0.02003033878 |

Frozen power fits give

\[
\boxed{p_{cross}=1.0056948923496356},
\]

\[
\boxed{p_{GG}=2.007490390559045},
\]

\[
\boxed{p_{joint}=1.0076444430189475}.
\]

Thus the physical sine ordering preserves the expected hierarchy

\[
\boxed{C_{cross}/D=O(\epsilon)},
\qquad
\boxed{C_{GG}/D=O(\epsilon^2)},
\qquad
\boxed{\Delta_{joint}=O(\epsilon)}.
\]

The final frozen endpoint is

\[
\boxed{\Delta_{joint}(1/64)=0.020030338775070305<0.05}.
\]

All preregistered conditions pass.

## Why this result matters

The historical two-node gate used

\[
H_+=(T+T^\dagger)/2
\]

and obtained a smaller finite endpoint, but the subsequent Lorentzian stack uses `K_sine=[V,H_E^sine]`. The present calculation removes that ordering mismatch.

The finite coefficients are genuinely different:

```text
||H_sine|| / ||H_plus||             ~= 1.3870
||[H0,H1]_sine|| / ||[H0,H1]_plus|| ~= 1.7124
joint_endpoint_sine / plus           ~= 1.3619.
```

Nevertheless the asymptotic exponents are essentially unchanged. Therefore the `O(epsilon)` / `O(epsilon^2)` HDA hierarchy is robust to the physically required Euclidean ordering, while the finite calibration is ordering sensitive.

## Remaining scope

This closes the **physical sine-order Euclidean two-node finite calibration** on the frozen expectation-metric route habitat.

It still does not include:

1. the phase- and real-normalization-completed `H_L` amplitudes in the same two-node commutator;
2. the full operator-first flux-metric square root on the complete geometry-changing habitat;
3. multiple independent WKB/habitat probes;
4. collective-spin/refinement scaling;
5. the uniform simultaneous-cutoff theorem.

The next gravity calculation is therefore the ordering-consistent

```text
H_E^sine + (1+beta^2) H_L + R_operator-first
```

two-node HDA, with all relative normalizations fixed upstream rather than fitted to the commutator.
