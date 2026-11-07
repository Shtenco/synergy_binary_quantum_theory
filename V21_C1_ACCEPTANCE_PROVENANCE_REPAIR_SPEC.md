# V21 C1 acceptance and provenance repair specification

Status: audit specification only. No worker or guard code has been modified by this document.

Target branch under review: `research/bqg-v21-c1-checkpoint-repair-20261107`  
Target code commit: `2c28123720db59306890daa0dfba1ad21c2e7eb1`

## 1. Current observed defects

1. `validate_physical_stage()` uses a generic set of physical leakage keys for every stage and requires only one matching key. It does not bind required diagnostics to the actual operation signature.
2. The guard permits an empty state and returns `stage_physical_basis_pass=true`. The test suite explicitly codifies this.
3. `read_stage()` trusts `meta.guard.stage_physical_basis_pass`; it does not recompute the guard from `meta.diagnostics`.
4. Stage 2/3 copy `previous["input_state_sha256"]` rather than using `previous["npz_sha256"]` as the digest of the immediate input checkpoint. `predecessor_ok()` then compares inherited input hashes for equality. This is not a direct stage-to-stage chain attestation.
5. Exact-zero V4/V2 states can compare with relative mismatch zero. The collector reports pair norm but does not distinguish a verified exact zero from a nonzero signal.
6. Stage metadata is not externally bound by a trusted manifest that pins its JSON digest, NPZ digest, source commit, backend digest, and operator signature.

These are workflow acceptance/provenance defects. They do not, by themselves, prove that any previously computed amplitude is mathematically wrong.

## 2. Required stage-specific diagnostic schema

The required schema must be selected from the actual operator signature, not inferred from a user-supplied stage number alone.

| Operator | Required physical-basis diagnostics | Separate diagnostics |
|---|---|---|
| `CV` | `CV_complete_basis_leakage` | finite amplitude check |
| `CK` | `outer_complete_basis_leakage`, `internal_volume_sector_leakage` | `complete_charge_basis_leakage` must be present and finite, but is not silently treated as a physical-basis leakage bound |

For `CK`, require both outer and internal-volume leakage keys. A single one of them is insufficient. For `CV`, require its own key. Reject unexpected operation signatures and reject a mismatch between the stage metadata and the operation derived from the frozen path specification.

Every supplied numeric diagnostic must be finite and nonnegative. The threshold comparison must use a documented strict convention; retain the current production threshold only if the physics protocol explicitly justifies it.

Charge leakage remains a separate unresolved covariance diagnostic. Until a charge-covariance theorem/gate is implemented, the result must say `charge_covariance_verified=false`; it must not be silently promoted to a passed physical condition.

## 3. Separate computation, acceptance, and signal statuses

Do not collapse these logically different outcomes into one boolean.

- `computation_status`: `COMPUTED` or `COMPUTATION_FAILED`.
- `physical_acceptance_status`: `PASS`, `FAIL`, or `NOT_EVALUABLE`.
- `amplitude_status`: `NONZERO` or `EXACT_ZERO`.
- `pair_status`: `VERIFIED_NONZERO_PAIR`, `VERIFIED_EXACT_ZERO_PAIR`, `MISMATCH`, or `INCOMPLETE`.

A zero intermediate state must not be used as a passed precursor for later stages: it should terminate that path as `EXACT_ZERO`. A zero final pair may be a legitimate mathematical cancellation, but it is not evidence of a nonzero physical action. It must be recorded distinctly and must not count toward a nonzero C1 signal.

Do not require every individual ordered contribution to be nonzero. Exact zeros can occur by symmetry. The nonzero condition belongs to the aggregate whose physical role is being claimed.

## 4. Correct direct checkpoint provenance

For each stage metadata object, store:

- `seed_state_sha256`: immutable digest of the original seed state.
- `input_state_sha256`: digest of the exact NPZ file consumed by this stage.
- `output_state_sha256`: digest of the exact NPZ written by this stage (same value as the saved artifact digest).
- `source_commit`: full frozen source commit.
- `volume_backend_sha256`: digest of the exact backend implementation.
- `operator_signature`, stage index, path identity and cutoff.
- `diagnostics`, state support and norm.
- A schema/version identifier for the metadata.

For stage 1, `input_state_sha256` is the canonical seed serialization digest. For stage N>1, enforce:

```text
current.input_state_sha256 == previous.output_state_sha256
```

Do not compare `current.input_state_sha256` to `previous.input_state_sha256`. Those values are expected to differ after an operation.

The collector must independently load and hash each NPZ and verify the complete predecessor chain. It must also verify that all stages share the same frozen source commit, graph, source node, omitted slot, cycle, cutoff, and backend hash, while allowing stage-specific input/output hashes to differ.

## 5. Trusted manifest and metadata integrity

Use a separate manifest, not a self-hash field inside the JSON being hashed. The manifest must pin, for every stage:

- path key (mode/order/aux/stage);
- JSON metadata SHA-256;
- NPZ SHA-256;
- source commit;
- backend SHA-256;
- direct input and output checkpoint SHA-256;
- operator signature.

The collector must compare actual file bytes against the manifest before parsing/using the metadata. The manifest itself must be retained as a workflow artifact and committed or otherwise anchored to a trusted run/commit. A checksum detects accidental or uncoordinated mutation; it does not authenticate a manifest if an attacker can replace both the manifest and all artifacts. If adversarial provenance is in scope, add a detached signature or an immutable trusted-run attestation.

## 6. Collector acceptance order

1. Validate manifest schema/version and trusted manifest identity.
2. Validate exact expected path census: 32 unique path keys.
3. Validate each metadata and NPZ digest against the manifest.
4. Parse metadata only after digest validation.
5. Validate source commit, graph, backend, cutoff, stage signature, and direct input/output hash chain.
6. Recompute stage diagnostics with the stage-specific guard. Ignore any stored `stage_physical_basis_pass` boolean except as a consistency field.
7. Reject failed stages and terminate exact-zero intermediate paths.
8. Reconstruct the four aggregate V4 states and project Gauss/scalar channels under the frozen rules.
9. Compare each against a distinct, independently produced V2 primitive; pin V2 metadata and NPZ hashes and source/orientation/cutoff.
10. Compute the Hermitian pair and report its norm and exact-zero/nonzero status.
11. Emit a machine-readable certificate with all 32 path statuses, all 4 V2 comparisons, residuals, support sizes, norms, and provenance digests.

## 7. Mandatory adversarial regression tests

- CV stage missing `CV_complete_basis_leakage` -> reject.
- CK stage with only outer leakage -> reject.
- CK stage with only internal leakage -> reject.
- CK stage missing charge diagnostic -> reject or explicitly mark charge check unevaluable; never claim full acceptance.
- CV stage supplied only CK diagnostics -> reject.
- Stage 1 metadata claims CK while path requires CV -> reject.
- Empty intermediate state -> `EXACT_ZERO`, not accepted precursor.
- Corrupt NPZ -> reject before loading.
- Modified diagnostics JSON with unchanged manifest -> reject.
- Modified metadata and modified NPZ but unchanged trusted manifest -> reject.
- Wrong source commit/backend/cutoff/operator signature -> reject.
- Stage 2 input hash not equal stage 1 output hash -> reject.
- Stage 3 input hash not equal stage 2 output hash -> reject.
- Two V2 references point to same file -> reject.
- V4 and V2 both exactly zero -> `VERIFIED_EXACT_ZERO_PAIR`, not nonzero signal.
- One side zero and other nonzero -> mismatch.
- Nonzero complex phase mismatch -> reject even if magnitudes match.
- Repeated path key or missing one of 32 paths -> reject.
- All 32 paths complete but aggregate pair zero -> report exact-zero pair and do not close nonzero C1.

## 8. Closure criteria

C1 is not closed by unit tests alone. Require all of:

- all 32 original physical paths completed on the frozen production source;
- all stage metadata and checkpoint chains pass independent collector verification;
- all four V2 references are distinct, independently generated, provenance-pinned calculations;
- V4/V2 complex-state comparison passes at the preregistered tolerance;
- the final pair is classified explicitly as nonzero or exact zero;
- full logs and artifacts are retained with reproducible hashes;
- a fresh replay from a clean checkout reproduces the result;
- no synthetic or mocked amplitudes appear in the physics evidence (mocks remain allowed in unit tests clearly labelled as such).

Even a successful C1 certificate is only one physical Lorentzian pair/precursor. It does not prove the full HDA, physical projector, TT graviton pole, or GR limit.
