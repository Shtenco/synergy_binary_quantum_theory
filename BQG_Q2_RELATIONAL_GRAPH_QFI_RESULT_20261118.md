# Q2 general finite graph relational-QFI control (2026-11-18)

## Scope and construction

A generic finite graph is represented by four-valent spin-1/2 SU(2)
intertwiners at vertices and normalized epsilon-singlet contractions on
selected pairs of half-links. Each internal edge receives one independent
three-axis relative twist, inserted at its first endpoint. All remaining
half-links are retained as the boundary state. The gate contracts the full
tensor network and computes its complete pure-state multiparameter QFI.

The first fixtures are a three-node path and a four-cycle. This tests graph
contraction beyond a single shared face or a single star, while retaining
explicit edge identity and all cross-edge QFI blocks.

## Acceptance checks

For each graph and random normalized invariant node states:
- boundary dimension and contraction are consistent;
- the full QFI is symmetric and positive semidefinite;
- each diagonal 3-by-3 edge block agrees with the isolated local identity
  control within numerical tolerance;
- off-diagonal cross-edge blocks are not discarded and are measured by
  their Frobenius norm.

The output deliberately does not turn the matrix into scalar edge
conductances. That requires a separately justified gauge quotient and
coarse-graining map.

## Limitations

This is a finite spin-1/2 tensor-network control on the explicit path and
four-cycle fixtures. It is not yet a proof for arbitrary graphs. It does not
include mixed-state boundary tracing, a complete multi-node gauge reduction,
higher-spin sectors, the BQG history-amplitude measure, the causal-PL
pushforward, or a continuum spectral-dimension result.
