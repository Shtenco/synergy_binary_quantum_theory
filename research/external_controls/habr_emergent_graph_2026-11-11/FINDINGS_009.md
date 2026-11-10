# Findings 009 — independent finite-size control of the proposed small-world scaling

## 1. Why this is worth testing

The article proposes a 3D local cubic graph and a distance-dependent nonlocal edge probability p_r=p/r², with p=1/(K N^(1/3)). Unlike the problematic spectral-zeta algebra, this scaling has a plausible graph-theoretic derivation once the ensemble is made explicit.

For a 3D lattice, the number of candidate vertices in a shell of radius r and width dr scales as C r² dr. If each candidate pair is linked with probability p/r², the expected number of long links contributed by each shell is approximately C p dr. Summing shells to linear size L=N^(1/3) gives:

    E[d_nonlocal] ~ C p L.

To keep the expected nonlocal degree finite and nonzero as N grows, require pL=O(1), hence:

    p ~ L^(−1) = N^(−1/3).

Thus p∝N^(−1/3) is not arbitrary under these assumptions: it is the unique power-law scaling that keeps the mean shortcut degree O(1) for a 3D shell measure r² and a per-pair kernel r^(−2). This does not derive d=3 or the r^(−2) kernel; those remain model inputs.

## 2. Explicit operationalization simulated

Because the article does not fully specify boundary conditions or the exact per-pair sampling law, the audit uses one natural implementation:

- periodic L×L×L cubic torus;
- six nearest-neighbour edges per vertex;
- for every nonlocal unordered pair at toroidal Euclidean distance r, add an independent edge with probability p/r²;
- p=1/(6L) for the claimed exponent;
- compute sampled mean shortest-path length and the combinatorial Laplacian's first positive eigenvalue.

Reproducible implementation: simulate_smallworld_scaling.py. Sizes L=8,12,16 (N=512,1728,4096), two independent seeds per size and per scaling exponent.

## 3. Finite-size result for the article's exponent alpha=1

For p=1/(6L), the six simulations give:

| L | N | Mean degree (range over seeds) | Mean path length (range) | Spectral gap (range) | gap·(ln N)² (range) |
|---:|---:|---:|---:|---:|---:|
| 8 | 512 | 6.95–7.00 | 3.99–4.09 | 0.9706–0.9756 | 37.77–37.97 |
| 12 | 1728 | 7.04–7.05 | 4.87–4.90 | 0.6872–0.7072 | 38.19–39.30 |
| 16 | 4096 | 7.10–7.16 | 5.37–5.40 | 0.6103–0.6274 | 42.22–43.41 |

The quantity gap·(lnN)² is roughly stable (mean ≈39.8, standard deviation ≈2.2 across these six finite graphs). This is consistent with a spectral-gap scaling of the form:

    Δ₁(N) ~ C / (ln N)²

for this particular graph ensemble and the tested size range. Average shortest paths also grow much more slowly than the cubic-lattice baseline 3L/4.

## 4. Scaling-exponent control experiment

Repeat the same construction with p=1/(6L^alpha), for alpha=0.5, 1.0, 1.5, 2.0. Across L=8,12,16 and two seeds:

| alpha | Mean degree (approx.) | Mean path (approx.) | Mean gap·(lnN)² | Trend with L |
|---:|---:|---:|---:|---|
| 0.5 | 8.6–10.5 | 3.3–4.0 | 111.0 (large growth) | shortcuts proliferate; gap decays more slowly than ln⁻² |
| 1.0 | 6.95–7.16 | 4.0–5.4 | 39.8 (stable) | consistent with ln⁻² gap |
| 1.5 | 6.28–6.30 | 4.7–7.0 | 21.9 (declining) | shortcuts vanish; gap decays faster |
| 2.0 | 6.07–6.08 | 5.2–8.6 | 17.6 (declining) | approaches local-lattice behaviour |

Within this finite-size operationalization, alpha=1 is the scaling that keeps mean degree O(1) and produces the clearest ln⁻² spectral-gap collapse. This is a positive, nontrivial result worth extending, not proof of the article's full theory.

## 5. Important caveats

- Only L≤16 and two seeds per alpha/size were tested; finite-size corrections and realization variance remain.
- The article does not explicitly specify toroidal boundary conditions, Euclidean distance, independence of edge trials, or the exact cutoff/normalization. Alternative choices must be tested.
- The observed coefficient C≈40 is not derived by the article's formulas and is not claimed to be universal.
- The simulation supports one graph-theoretic scaling claim only. It does not validate the spectral-zeta formula for lnN, particle masses, constants, TSCO, or emergent Lorentzian geometry.

## 6. Variation-script variable mismatch

The separate variation_principle.py sets k=xN^(2/3) and K_eff=K+lambda*k, producing k≈2.9×10^79 at the minimum for N≈4.2×10^121. Under the explicit graph ensemble above, the mean nonlocal degree is O(1), not O(N^(2/3)). If k denotes total shortcut count, adding it directly to a vertex degree K is not valid; if it denotes per-vertex degree, its N^(2/3) scaling contradicts the sparse p_r law. The article must define k's graph-theoretic meaning and reconcile these two constructions.

## 7. Verdict

- p∝N^(−1/3): DERIVABLE CONDITIONALLY from a 3D shell measure and r^(−2) per-pair link kernel by finite expected shortcut degree.
- Δ₁∼1/(lnN)²: PRELIMINARY FINITE-SIZE SUPPORT in one explicit implementation.
- Unique universality class / asymptotic theorem: NOT ESTABLISHED.
- Compatibility with variation_principle.py's k∼N^(2/3): UNRESOLVED / APPARENTLY INCONSISTENT.
- This is the first substantive positive control result found so far; it should be retained and tested, not discarded with the failed zeta derivation.

## Correction added 2026-11-11 — superseded spectral-gap interpretation

The initial claim in sections 3–4 that the tested gap·(lnN)² was approximately stable was based on only L=8,12,16 and two seeds. It does not survive the larger control. See FINDINGS_011.md and sparse_scaling_control.py: for L=16,24,32,48 (three seeds each), the sampled λ₂ decreases toward a plateau while λ₂·(lnN)² rises; the exact expected-Laplacian gap through L=512 tends to approximately 0.7284. Therefore the preliminary λ₂∝(lnN)⁻² interpretation is withdrawn. The conditional expected-shortcut-degree argument for p∝N⁻¹/³ remains valid under the specified 3D r⁻² ensemble. Do not cite the original small-size gap table as support for a logarithmic gap law.
