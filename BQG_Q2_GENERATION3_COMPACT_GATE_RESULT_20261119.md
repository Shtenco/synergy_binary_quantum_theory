# Q2 generation-3 compact causal-to-PL gate — result (2026-11-19)

## Execution status

The compact generation-3 implementation was added to the mandatory
`.github/workflows/q2-pl-graph-export.yml` workflow and executed successfully.

- CI run: https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44434874944
- Commit: `018d73bf3a8edc3b8bc099ed6e53e2f9ce6c5765`
- Result: **PASS**
- Regression suite: 121 passed across 18 pytest invocations; zero failures.

The workflow now runs `tests/test_q2_pl_generation3_stream.py`, exports
`Q2_PL_GENERATION3_COMPACT_GATE.json`, and uploads it with the other Q2
verification artifacts.

## Exact finite combinatorial result

| Quantity | Generation 3 |
|---|---:|
| Parent tetrahedra | 9,216 |
| Child tetrahedra | 221,184 |
| PL dual edges | 442,368 |
| Unique PL dual edges | 442,368 |
| Causal fibres | 512 |
| Within-parent dual edges | 331,776 |
| Across-parent dual edges | 110,592 |
| Connected components | 1 |
| Pushforward counting-measure total | 16 |

The four route classes each contain exactly 110,592 dual edges. Every child
has degree four and exactly one neighbour in each route class. No duplicate
dual edges were emitted.

The direct construction was compared with the incidence-based reference:
- Generation 1: 768/768 edge-and-route signatures agree; zero missing and
  zero extra signatures.
- Generation 2: 18,432/18,432 edge-and-route signatures agree; zero missing
  and zero extra signatures.

All 14 generation-3 gate checks passed, including parent face pairing,
four-regularity, route incidence, connectedness, causal-fibre surjectivity,
and parent pushforward mass conservation.

## Nonuniform causal-fibre distribution

The 512 fibre sizes have 31 distinct values. The exact histogram is:

| Fibre size | Number of causal edges |
|---:|---:|
| 64 | 8 |
| 128 | 8 |
| 144 | 8 |
| 192 | 32 |
| 208 | 8 |
| 240 | 24 |
| 256 | 16 |
| 288 | 16 |
| 304 | 8 |
| 320 | 24 |
| 336 | 16 |
| 352 | 40 |
| 384 | 40 |
| 416 | 40 |
| 432 | 8 |
| 448 | 8 |
| 480 | 48 |
| 512 | 16 |
| 528 | 24 |
| 544 | 16 |
| 576 | 16 |
| 608 | 8 |
| 624 | 16 |
| 640 | 8 |
| 672 | 8 |
| 704 | 8 |
| 752 | 8 |
| 832 | 8 |
| 848 | 8 |
| 928 | 8 |
| 992 | 8 |

Thus the minimum is 64, maximum 992, and mean 432. The max/min ratio is
15.5. This heterogeneity is explicitly reported; no inverse-fibre
reweighting or smoothing is introduced.

With the finite counting pushforward convention,
`mu_3(e) = |Phi_3^{-1}(e)| / 24^3`, the total mass is exactly 16 and is
conserved from every generation-2 parent fibre to its eight causal children.
This is a finite combinatorial measure only, not a physically selected
continuum measure.

## Implementation and regression changes

- `scripts/q2_pl_generation3_stream.py`: compact fixed-width NumPy arrays
  for child tetrahedra, dual endpoints, route colours, and causal-fibre IDs;
  exports the full fibre histogram and measure total.
- `tests/test_q2_pl_generation3_stream.py`: checks the generation-3 graph,
  exact small-level incidence agreement, fibre histogram normalization,
  mean fibre size, and pushforward mass.
- `.github/workflows/q2-pl-graph-export.yml`: runs the new regression and
  includes the compact generation-3 JSON in the artifact.

## Strict scope boundary and next gate

This result establishes a finite generation-3 combinatorial dual graph and
causal-fibre correspondence for the frozen q=2 and barycentric refinement
rules. It does **not** establish a physical QFI conductance, physical
probability measure, continuum operator, spectral-dimension limit, causal
spacetime, or Einstein dynamics.

The next technical gate is a compact matrix-free generation-3 Feshbach
calculation. It must independently validate the fibre projection
(`U^T U=I`, `P^2=P=P^T`), compare matrix-free products against the existing
generation-1/2 dense controls, establish solver residual and gap bounds, and
report the correction and self-energy derivative without allocating a dense
221,184-by-512 embedding. No generation-3 Feshbach result is claimed here.


## Matrix-free Feshbach primitive controls (2026-11-19)

A generic matrix-free operator layer has now been added:

- `scripts/q2_matrix_free_feshbach.py`
- `tests/test_q2_matrix_free_feshbach.py`

It applies the normalized fibre isometry (U), its adjoint (U^T), the
orthogonal projectors (P=UU^T) and (Q=I-P), the leakage maps
(B=QLU), (B^T=U^TLQ), and the Galerkin block (A_{PP}=U^TLU) without
allocating the fine-by-coarse dense embedding. It solves with (QLQ+P)
using sparse-graph matrix-vector products and conjugate gradients, then
exposes matrix-free applications of (B^T(QLQ|_Q)^{-1}B) and
(B^T(QLQ|_Q)^{-2}B).

The six new tests compare the matrix-free implementation against explicit
dense controls on a four-vertex cycle. They cover isometry/projector
identities, adjointness of (B) and (B^T), the Galerkin block, the
Feshbach correction, the second-resolvent self-energy slope, and invalid
graph/fibre inputs.

- CI run: https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44435219895
- Result: **PASS**, 127 tests across 19 pytest invocations, zero failures.

**Strict boundary:** these are exact small-graph controls for the generic
matrix-free primitives. They are not yet a generation-3 Feshbach computation.
The 221,184-vertex generation-3 graph has not yet been passed through the
matrix-free solver. Before that step, the implementation still needs
generation-1/2 cross-checks against the existing Feshbach diagnostics,
generation-3 residual/gap controls, solver iteration reporting, and a
scalable coarse-operator norm/eigensolver strategy.
