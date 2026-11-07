# V21 C1 external stage attestation format

Status: implementation contract for the fail-closed V21 collector. This document does not claim that C1 has been computed or that the physical HDA is proved.

## Why an external attestation is required

A stage JSON file cannot authenticate itself. The collector therefore requires a separate manifest entry that pins the exact bytes of both the stage metadata JSON and its NPZ checkpoint. The pin is checked before the metadata is trusted. The metadata's own `npz_sha256` and `output_state_sha256` must agree with the actual checkpoint bytes.

A plain hash manifest protects against accidental substitution only when the manifest is stored in a separately controlled location. If the threat model includes an actor able to replace both stage files and the manifest, the manifest must additionally be signed or placed in immutable, access-controlled storage. Hashes alone do not authenticate their own authority.

## Per-path manifest schema

Every entry in the collector's `paths` array must include:

```json
{
  "stage1": "artifacts/path-A/stage1",
  "stage2": "artifacts/path-A/stage2",
  "stage3": "artifacts/path-A/stage3",
  "stage_attestations": {
    "stage1": {
      "metadata_sha256": "<sha256 of exact stage1 JSON bytes>",
      "npz_sha256": "<sha256 of exact stage1 NPZ bytes>"
    },
    "stage2": {
      "metadata_sha256": "<sha256 of exact stage2 JSON bytes>",
      "npz_sha256": "<sha256 of exact stage2 NPZ bytes>"
    },
    "stage3": {
      "metadata_sha256": "<sha256 of exact stage3 JSON bytes>",
      "npz_sha256": "<sha256 of exact stage3 NPZ bytes>"
    }
  }
}
```

The top-level manifest must also contain `seed_state_sha256`, `source_commit`, `source_node`, `omit`, `cycle`, `Jmax2`, and `volume_backend_sha256`, `volume_backend`, and `volume_backend_sha256`, alongside the existing `independent_v2` references.

## Required chain invariants

1. Stage 1 `input_state_sha256 == seed_state_sha256`.
2. For stage N>1, `input_state_sha256 == previous.output_state_sha256`.
3. Each stage `output_state_sha256 == npz_sha256 == SHA256(actual NPZ bytes)`.
4. Each attestation pins the SHA256 of the exact metadata JSON and exact NPZ bytes.
5. The collector recomputes the operation-specific guard from diagnostics and the loaded complex amplitudes. A stored `guard` flag alone is never accepted.
6. CV requires `CV_complete_basis_leakage`. CK requires `outer_complete_basis_leakage`, `internal_volume_sector_leakage`, and `complete_charge_basis_leakage`. Charge leakage is reported separately and is not silently treated as proof of charge covariance.
7. An exact-zero intermediate stage is recorded as a stopped path, not a passed precursor. An exact cancellation in the final assembled pair remains mathematically possible and is not rejected merely for being zero.

## Scope

Passing these software checks would establish provenance and diagnostic consistency for the accepted finite computation. It would not establish that the computation is physically complete, prove the full Hamiltonian constraint algebra, establish a physical TT sector, or derive Einstein dynamics. The 32 ordered paths and four independent V2 comparisons must still be run on the real backend and compared at the full complex-amplitude level.
