# Depth-6 numerical completion status

Branch: `bqg/depth6-numerical-completion`

Base: `3d48d8ba6c190efe23ee655528365f5045216af7`

## [3,2] execution pipeline

The isolated workflow `.github/workflows/bqg-depth6-32-rematerialize.yml` now performs:

1. Extract the immutable Stage-A bundle from commit `887e76fba5210524650002da0b54653d2c5b37fd` and verify bundle blob `cbbaa80f6b6b61d353458db1fa1f03f824730532`.
2. Deterministically regenerate the 2,755-block / 130,903-column assignment identity.
3. Cross-check all 15 persisted identities and all 102 recoverable historical shard records.
4. Run 112 KEYMETA-first shards. Each shard writes KEYMETA first, then persists the corresponding raw master matrices.
5. Verify all shard manifests, all raw file SHA-256 digests, exact orbit coverage, matrix dimensions, and q-key consistency.
6. Run the fail-closed unique-q peeling rank gate.

The complete assignment identity and raw master matrices are not themselves a rank certificate. The peeling stage certifies only blocks whose unique-q submatrix has full column rank. If a residual coupled system remains, the current gate deliberately refuses to assemble a global SVD until the frozen engine's cross-block q-row coordinate identity is independently established. Equal row counts alone are not accepted as proof of row alignment.

## Current claim boundary

No numerical result is recorded as passed until a workflow artifact contains the replayed certificate. The current finite depth-6 theorem remains `NOT_YET_PROVED`.

## Independent sectors still pending

- `[3,1,1]`: 2,719 blocks, 153,455 columns.
- `[2,2,1]`: 2,749 blocks, 130,503 columns.

These require separate assignment identities, independent KEYMETA/raw replay, and their own rank witnesses. A [3,2] result cannot be reused as evidence for either sector.

## Integrity note

The previous Stage-A workflow artifacts are no longer available through the artifact API. The 15 persisted raw records and the 102-row historical salvage ledger remain available as cross-check inputs. The full assignment identity is regenerated from the pinned engine, not reconstructed from the incomplete persisted ledger.
