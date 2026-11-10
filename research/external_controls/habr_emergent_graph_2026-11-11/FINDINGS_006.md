# Findings 006 — independent recomputation of the 82-entry constant dictionary

Source: physics_constants/stress_test.py, blob SHA b617ac8f42bec5e7916dfefce576665bb68e81fc.

## 1. Recalculation scope

The formula dictionary was independently transcribed and evaluated for all four N conventions. The source contains 83 formula outputs and 82 target constants; B_s_meson is an unmatched formula output. The target values below are copied from the author's own constants dictionary, not yet independently refreshed from CODATA/PDG or scheme-specific metrology.

## 2. Full-sample fit of N

Reimplementing the stated sum of squared log-ratios over the 82 target entries gives a fitted N of approximately 4.1967710718×10^121 (lnN≈280.04711169), with objective ≈0.00118113. This is very close to the JSON's hard-coded 4.1977×10^121 and materially below the geometric-resonance root 4.4759473527×10^121.

This is direct evidence that the code's fit calibrates N against the same target constants. It is a legitimate calibration operation, but it is not an independent derivation of N from graph geometry.

## 3. Aggregate relative-error metrics at fixed N

Independent formula transcription gives the following preliminary metrics against the author's 82-entry target dictionary:

| N convention | Mean relative error | Median relative error | <0.1% | <0.5% | <1% |
|---|---:|---:|---:|---:|---:|
| Exact geometric R=1 root | 1.5051% | 1.8488% | 7/82 | 22/82 | 24/82 |
| stress_test.py N_phys | 0.2811% | 0.2210% | 26/82 | 68/82 | 80/82 |
| mass_relation.json N | 0.2841% | 0.2182% | 19/82 | 69/82 | 80/82 |
| variation-script N | 0.2856% | 0.2142% | 19/82 | 69/82 | 80/82 |

The exact R=1 root performs substantially worse against the author's own constants than the corrected/fitted N values. This is an important distinction: the advertised high accuracy is associated with the adjusted N, not the direct geometric resonance root.

The article reports mean error 0.26%, median 0.17%, 24 entries below 0.1%, 62 below 0.5%, and 79 below 1%. The independent transcription does not reproduce those counts. A final reproduction must execute the frozen script in a clean environment and reconcile every row, but the discrepancy is already material.

## 4. Recomputed leave-one-out structure

Re-fitting N separately after excluding each of the 82 target entries gives preliminary results:

- fitted N range: 4.1943889790×10^121 to 4.1990906348×10^121;
- mean held-out relative error: 0.2877%;
- median held-out relative error: 0.2185%;
- 80/82 held-out values below 1%;
- 69/82 below 0.5%;
- 20/82 below 0.1%.

This procedure is a genuine leave-one-observable-out calibration test, but it is not a prediction from a fixed theoretical N: every held-out target is predicted after N has been re-estimated from the other 81 targets. It also leaves highly correlated observables in the training set.

## 5. Key implication

The direct geometric formula gives lnN≈280.1115. The code applies a correction to lnN and the full-data optimizer returns lnN≈280.0471. The formula accuracy rises sharply after this adjustment. Unless the correction and fitted value are derived independently from the graph dynamics, the high-accuracy table cannot be described as a parameter-free prediction.

## 6. Current verdict

- Formula dictionary arithmetic: independently transcribed and evaluated for 82 targets / 83 outputs.
- Full-data fitted N: reproduces a value near the JSON N.
- Exact geometric-root prediction quality: materially worse than adjusted/fitted N.
- Article aggregate accuracy statistics: not exactly reproduced.
- LOO: calibration-based, not blind prediction; correlated observables remain.
- Physical validity of target values, dimensions, and formula derivations: still open.