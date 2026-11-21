# Depth-6 numerical completion workspace

This branch isolates all new finite depth-6 numerical work from `main`.

Base commit: `3d48d8ba6c190efe23ee655528365f5045216af7`.

## Frozen problem

The finite depth-6 shell and S5 decomposition are already fixed:

- Gauss-admissible assignments: 264,962
- Hilbert dimension: 3,111,637
- S5 spin orbits: 2,757
- Remaining numerical sectors:
  - `[3,2]`: multiplicity 130,903; 2,755 blocks
  - `[3,1,1]`: multiplicity 153,455; 2,719 blocks
  - `[2,2,1]`: multiplicity 130,503; 2,749 blocks

Structural support is closed for all seven S5 irreps. Do not recompute shell enumeration, multiplicities, Jucys selectors, branch identities, or the structural matching/capacity proofs.

## Required numerical certificates

Each sector must have an independent, reproducible certificate generated from pinned source inputs. At minimum, preserve:

1. Exact input/operator definition and source commit/blob hashes.
2. Complete block/column coverage and an independently checked assignment identity.
3. Actual-q coverage evidence; identity rematerialization alone is insufficient.
4. Rank evidence from a documented sparse QR/SVD or equivalent certified procedure.
5. Explicit smallest singular-value / rank-threshold policy, residual norms, and backward-error bounds.
6. A kernel witness check: either an empty kernel with a certified lower bound, or explicit basis vectors with residual checks.
7. Deterministic replay instructions and machine-readable JSON output.

The `[3,2]` identity-rematerialization gate currently has status `PASS_IDENTITY_REMATERIALIZATION_ONLY`. It does **not** certify actual-q coverage, numerical rank, a master kernel, or closure. Do not promote it to a numerical certificate.

## Completion order

1. Recover surviving Stage-A numerical evidence for `[3,2]`; compute only genuinely missing master-map pieces.
2. Produce independent numerical witnesses for `[3,1,1]` and `[2,2,1]`.
3. Run cross-sector and full-frontier consistency checks.
4. Update the frontier only after the certificates independently replay.

## Claim boundary

Until all three numerical witnesses pass, keep the full finite depth-6 theorem at `NOT_YET_PROVED`. Finite depth-6 closure, if later established, does not establish continuum physicalization or experimental confirmation.

## Execution environment limitation in this session

The GitHub repository was reachable through the GitHub integration, but direct local `git clone` was unavailable because this execution environment could not resolve `github.com`. No local numerical computation is claimed by this branch initialization. The dedicated branch is ready for the isolated computation workflow; numerical status remains unchanged until reproducible certificates are generated and replayed.
