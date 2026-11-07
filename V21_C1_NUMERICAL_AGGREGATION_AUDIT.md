# V21 C1 numerical aggregation audit

Status: code-level comparison of V21 staged aggregation with the historical V2 ordered-triple implementation. No physical computation was run by this audit.

## 1. Historical pruning points

The historical ordered-triple path uses `TOL = 1e-11` and `add(dst, src, scale, tol)` when combining the eight auxiliary `(i,j,k)` contributions. The add operation removes a key whenever the current partial sum has magnitude less than or equal to the threshold. A later auxiliary contribution can then re-create the key from zero.

The independent scalar Gauss projection also uses `LGA.TOL = 1e-11`. The V21 collector must not use the triple worker's looser `W.TOL = 1e-10` for this comparison. It now explicitly uses `LGA.TOL`.

This ordered thresholding is not algebraically equivalent to summing all eight contributions first and pruning once. Example: two contributions of magnitude `0.75e-11` to the same key are each dropped by the historical incremental accumulator, while an unpruned sum retains `1.5e-11`.

## 2. Two comparison lanes

The V21 collector now calculates both:

- **Unpruned lane:** stable complex summation using separate real/imaginary `math.fsum`; no support threshold during auxiliary aggregation; then scalar Gauss projection at `LGA.TOL`.
- **Legacy-aligned lane:** ordered incremental accumulation reproducing the historical `add(..., tol=LP.TOL)` behavior; then scalar Gauss projection at `LGA.TOL`.

The four independent V2 references are gated against the legacy-aligned lane at relative delta below `1e-10`. Unpruned-versus-V2 deltas are separately reported. Both lanes must satisfy the nonscalar rejection gate.

The final V21 physical pair is assembled from the unpruned lane, with stable summation and no tolerance pruning. Its `VERIFIED_EXACT_ZERO_PAIR` label means zero in the resulting floating-point complex sum only; it is not a symbolic proof of exact cancellation.

## 3. Why this matters

If only the legacy-aligned lane is examined, a sub-threshold nonzero residue can disappear by construction. If only the unpruned lane is compared to a legacy V2 artifact, threshold-induced differences can be misdiagnosed as an operator mismatch. Reporting both lanes separates these two questions.

## 4. Remaining caveats

- The two lanes still inherit floating-point arithmetic inside the underlying SU(2) operator kernels.
- The unpruned lane is not arbitrary precision and does not establish a rigorous interval bound.
- The independent V2 artifact metadata must remain externally pinned and semantically matched.
- The 32-path computation, four V2 comparisons, and all structural tests still require actual execution on a working runner.
- A successful software comparison is not a proof of the full Hamiltonian constraint algebra, physical TT closure, or Einstein dynamics.

## 5. Underflow and overflow handling

State and pair L2 norms use iterative `math.hypot`, not `sqrt(sum(abs(z)**2))`. This prevents a nonzero amplitude such as `1e-200` from being mislabeled as zero because its square underflows. Pair zero classification is based on whether the assembled amplitude map is empty after exact floating-point zero removal, not whether a squared norm happens to evaluate to zero.

Relative V2 comparison no longer uses an arbitrary `1e-15` denominator floor, which could conceal a complete mismatch between tiny states. If any comparison norm or delta overflows to a non-finite value, the relative delta is set to infinity so the acceptance comparison fails closed. Projection and final pair norms are also required to be finite.
