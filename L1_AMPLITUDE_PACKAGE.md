# L1_AMPLITUDE_PACKAGE v1

## Purpose and status

This schema binds sparse complex amplitudes to their operator implementation,
Hilbert carrier, background, shard provenance, and aggregation procedure.
Version 1 is deliberately an **E-only regression precursor**. It retains all
72 source shards and all 72 second-action shard states, as well as the three
reconstructed source and response states. It is not a complete physical L1
package and cannot be consumed as evidence for the joint E+S+R carrier, metric
calibration, Lorentzian Schur reduction, or HDA closure.

## Files

- `L1_AMPLITUDE_PACKAGE.json`: manifest, provenance, diagnostics, and status.
- `L1_AMPLITUDE_PACKAGE.npz`: lossless sparse-state arrays for the three
  representative source states (u_e) and reconstructed responses (v_e).

The JSON manifest records the SHA256 of the NPZ archive and producer source.
The package digest is computed over canonical JSON before adding the
`package_sha256` field. The source commit is caller-supplied and must be the
exact `git rev-parse HEAD` at production time; the producer does not independently
query Git. The legacy shard metadata also does not self-attest its generating
source blob, so worker identity is pinned by the schema's known source blob but
is not cryptographically attested by each shard.

## Sparse state encoding

Each state is represented by three arrays:

- `spins`: int16 matrix, one full global spin key per row;
- `Ks`: int16 matrix, one full intertwiner key per row;
- `amp`: complex128 vector, one raw complex amplitude per row.

Rows are unique and sorted by the byte encoding of the complete state key.
No phase alignment, real projection, orbit averaging, or intermediate
amplitude pruning is permitted in the raw archive.

## Provenance invariants

A consumer must reject the package if:

1. schema or version is unknown;
2. the NPZ SHA256 does not match the manifest;
3. any source shard is absent, duplicated, malformed, non-finite, or has a
   mismatched edge/node identity;
4. shard backend, parent block, spin cutoff, or carrier dimensions disagree;
5. source states are not identical within the frozen tolerance;
6. aggregation regression against the historical collector exceeds tolerance;
7. any required operator is marked missing when a full physical package is
   requested.

The producer recomputes the Git blob SHA-1 of its worker and collector source
files and rejects any mismatch with the frozen source identities. The manifest
retains source Git commit, worker/collector/producer blob identifiers, and
SHA256 hashes of all 216 consumed shard JSON/NPZ files. The verifier checks the
216 source files against their SHA256, the producer Git blob and SHA256
identities, the complete archive member set, and both package and archive
hashes.

## Aggregation

All contributions for a state key are accumulated before thresholding.
Real and imaginary components use deterministic compensated summation.
A single global threshold is applied after complete aggregation; both the
discarded amplitude norm and its ratio to the full pre-threshold aggregate norm
are recorded. The package fails closed if this relative loss exceeds the frozen
regression tolerance (3\\times10^{-7}). A consumer must not silently
re-threshold the archive.

The legacy collector is retained as a regression oracle only. Its intermediate
threshold behavior is not the canonical aggregation rule.

## Execution

Run the deterministic aggregation self-test:

```bash
python scripts/collective_l1_amplitude_producer.py --self-test
```

Build the E-only package from a complete directory containing all 72 shard
triplets:

```bash
python scripts/collective_l1_amplitude_producer.py \
  --root <shard-directory> \
  --source-commit <full-40-character-git-sha> \
  --output-dir verification_results/l1_amplitude_package
```

Verify the emitted manifest and archive:

```bash
python scripts/collective_l1_amplitude_producer.py \
  --verify-manifest verification_results/l1_amplitude_package/L1_AMPLITUDE_PACKAGE.json
```

The producer compares each stable aggregate with the historical incremental
collector and also compares the full (K,A,B) moment matrices. A successful
self-test is only a unit-level numerical control; the E-only regression status
requires the real 72-shard input set.

## Reproduction commands

After all 72 shard triplets have been produced into `SHARD_DIR`, run from the
repository checkout:

```bash
SOURCE_COMMIT="$(git rev-parse HEAD)"
python scripts/collective_l1_amplitude_producer.py \
  --root SHARD_DIR \
  --output-dir verification_results/l1_amplitude_package \
  --source-commit "$SOURCE_COMMIT"

python scripts/validate_l1_amplitude_package.py \
  --manifest verification_results/l1_amplitude_package/L1_AMPLITUDE_PACKAGE.json \
  --archive verification_results/l1_amplitude_package/L1_AMPLITUDE_PACKAGE.npz \
  --producer-source scripts/collective_l1_amplitude_producer.py \
  --shard-root SHARD_DIR
```

The producer's regression comparison is against the historical collector and
uses the frozen relative tolerance (3\\times10^{-7}). The stable assembler
uses a single post-sum absolute amplitude threshold (10^{-10}), matching the
legacy cutoff value but not its intermediate-pruning semantics. The discarded
global L2 norm is reported. Neither tolerance is adjusted after observing the
result.

## Status semantics

- `PASS_E_ONLY_REGRESSION_PRECURSOR`: the E-only ingestion and regression
  checks passed. This is not a full physical PASS.
- `FAIL_E_ONLY_REGRESSION`: input, provenance, aggregation, or regression
  failed.
- A future complete package must use a new status and include E, S, and R on
  one shared carrier and one shared (W_l), plus independently measured metric
  calibration, full Schur blocks, and direct HDA quantities.

## Required future extension

A complete v2 package must include:

- raw (E,S,R) and all ordered depth-two words on the same Hilbert space;
- the common target-independent (W_l), all singular values, rank and
  condition number;
- (eta_E,eta_S,eta_R) and every ordered-pair leakage;
- independent (M_{hq}^{(l)}) or (B_l) with full six-value spectrum;
- full (P/Q) Lorentzian block, (D_l) spectrum, gap, conditioning, Schur and
  direct-solve residuals;
- microscopic and compressed HDA residuals, generator leakage, and structure
  function discrepancy;
- one immutable provenance manifest binding every component to the same
  background, carrier, operator versions, and source commit.

No field may be populated by importing a result from a different carrier or
refinement level.
