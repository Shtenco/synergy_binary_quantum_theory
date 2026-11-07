# BQG v17 — HDA evidence replay and CI blockage (2026-11-07)

Scope: read-only replay of prior K5/W5 finite numerical data; **NO NEW SU(2) LORENTZIAN AMPLITUDES** and **NO full off-shell HDA closure**.

## Confirmed numerical checks

- Refit published K5 five-point `||X(eps)||^2=A+B eps+C eps^2` using rounded X norms yields coefficients `A=2718.6183589770153`, `B=385.3687143979623`, `C=23.667376245836316`. Rounded-input agreement with v0.15 full-precision archived coefficients `[2718.618358974829,385.3687144585335,23.66737605353228]` is within 2e-7 absolute.
- From archived full-precision coefficients: `||X0||=52.14037168044383`; `||X1||=4.864912748809816`; cosine of Hilbert-space overlap `0.7596215557523616`.
- Under the EXPLICIT v0.15 fixed-cutoff parity and bounded-Lorentzian-channel hypotheses: `liminf ||A_full|| >= 34.76024778696255`. From archived `lim eps||D||=40.718968061580426`, conditional `liminf Delta_full/eps >= 0.8536622965099131`. This says nothing unconditional about general refinements or continuum.
- Independent W5 Euclidean spectrum: five singular values `[5.693291183213442,5.728740480166856,5.798524699282443,5.8412914834963,6.074945744003254]`; condition number `1.067035840695293`. For `|lambda|=1`, any `||B||<=5` entails `sigma_min(A+B)>=0.693291183213442`. **B has NOT been measured in the same common habitat**.
- Charged E-only 2-state Gram eigenvalues `[1.463343482394494,4.86283333633768]`, singular values `[1.2096873490263895,2.2051832886038474]`, condition `1.8229365549525482`. It is a distinct input space from W5; never concatenate those Grams.
- Checked S5 irrep multiplicity-weighted dimension sum: `3,111,637`. On `main/depth6_frontier.json` FOUR numerical sectors are closed at stated tolerances `[5],[4,1],[2,1,1,1],[1^5]`; their weighted dimension is `883,877`, and the remaining `[3,2],[3,1,1],[2,2,1]` weigh `2,227,760`. These ratios are NOT fractions of physical proof completion.

## V16 CI triage

The existing `shell: bash` workflow contained THREE PowerShell `if ($LASTEXITCODE -ne 0) { exit 2 }` statements. Fixed at commit `f68c4d07d38e843e08b9fa64814512aec84e0a05` on `research/bqg-v16-first-leg-20261107`. It triggered run [42367951178](https://github.com/Shtenco/binary_quantum_theory/actions/runs/42367951178), which nevertheless ended `failure` with `jobs[].steps=null` and ZERO artifacts. The reason remains unconfirmed. **No physical first-stage C(V) result was obtained.** Do not repeatedly restart expensive runs before establishing why jobs cannot expose steps.

## Mandatory remaining gates

1. Obtain 32/32 original true V4 paths for one physical pair on a single frozen SHA; independently compare all sparse complex V2 coefficients; no synthetic/mock outputs as science.
2. Build full Lorentzian B on the exact W5 used for A and compute `M(lambda)=(A+lambda B)^dagger(A+lambda B)` at a predeclared lambda, including output sectors and leakage.
3. Compute E/S/R compositions on a common graph-changing habitat; compare their full HDA commutators with an **independent** D, absolute AND relative residuals, and controlled cutoff/geometry refinement.
4. Only then physical projector, relational clock and interacting TT tests.

Sources: `main/depth6_frontier.json`, `main/BQG_DEPTH6_2111_MASTER_CLOSED_2026-10-01.json`, `research/node1_fullstate_gate_2026_10_31/BQG_DEEP_COHERENT_HE_FULL_AUDIT_20261106.md` in `research/bqg-node1-lorentzian-columns-20261031`, plus v0.14/v0.15 preserved reports.
