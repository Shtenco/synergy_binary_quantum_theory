# L1_AMPLITUDE_PACKAGE v1

## Purpose and status

This schema binds sparse complex amplitudes to their operator implementation,
Hilbert carrier, background, shard provenance, and aggregation procedure.
Version 1 is deliberately an **E-only regression precursor**. It is not a
complete physical L1 package and cannot be consumed as evidence for the joint
E+S+R carrier, metric calibration, Lorentzian Schur reduction, or HDA closure.

## Files

- `L1_AMPLITUDE_PACKAGE.json`: manifest, provenance, diagnostics, and status.
- `L1_AMPLITUDE_PACKAGE.npz`: lossless sparse-state arrays for the three
  representative source states (u_e) and reconstructed responses (v_e).

The JSON manifest records the SHA256 of the NPZ archive. The package digest is
computed over canonical JSON before adding the `package_sha256` field.

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

The manifest must retain source Git commit, worker and collector source blob
identifiers, and hashes of all consumed shard JSON/NPZ files.

## Aggregation

All contributions for a state key are accumulated before thresholding.
Real and imaginary components use deterministic compensated summation.
A single global threshold is applied after complete aggregation; discarded
amplitude norm is recorded. A consumer must not silently re-threshold the
archive.

The legacy collector is retained as a regression oracle only. Its intermediate
threshold behavior is not the canonical aggregation rule.

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
