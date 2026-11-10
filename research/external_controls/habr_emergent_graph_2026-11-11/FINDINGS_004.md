# Findings 004 — neutron lifetime as a candidate falsification test

## 1. Formula and independent numerical values

The article proposes tau_n = π ln N. This is dimensionally incomplete unless an implicit time unit is specified, but its numerical prediction can be evaluated under the article's convention.

| N convention | ln N | π ln N (s, article's numerical convention) |
|---|---:|---:|
| Geometric R=1 root | 280.1115142805 | 879.99627545 |
| stress_test.py N_phys | 280.0442211596 | 879.78486788 |
| mass_relation.json | 280.0473330085 | 879.79464404 |
| variation script | 280.0478807776 | 879.79636490 |

The article's quoted 879.78 s is consistent with the stress-test N_phys / JSON N conventions, not with the exact R=1 geometric root (about 880.00 s).

## 2. Update the experimental reference

The article and stress-test dictionary use 877.8 s with ±0.4 s stated in the article's outlook. The 2026 Particle Data Group listing instead reports the current neutron mean-life average as 878.3 ± 0.4 s. The PDG notes that this average is based on eight of the best nine ultracold-neutron measurements and that including the in-beam measurement gives 878.6 ± 0.6 s. Source: PDG 2026, neutron mean life, pp. 3–4.

Against 878.3 ± 0.4 s, the article's 879.785 s prediction differs by about 1.485 s, or 3.71 quoted experimental standard deviations (before any theory uncertainty). The direct geometric-root value 879.996 s differs by about 1.696 s, or 4.24 experimental standard deviations.

Against the older 877.8 s reference, 879.785 s differs by about 1.985 s, approximately 4.96 times 0.4 s. This explains the article's “5 sigma” wording, but it is not the current PDG average.

## 3. Interpretation

This is potentially a useful falsification target because tau_n is an observable not reducible to a definition like h=2πℏ. But it is not yet an independent prediction in the published computational workflow: neutron lifetime is included in the constants dictionary used by find_optimal_N() and in each leave-one-out refit except when neutron lifetime itself is excluded. A proper test must freeze N from a predeclared non-lifetime derivation, freeze the lifetime formula, and compare to an independent future measurement with its measurement method and systematic uncertainties specified.

The claimed significance must not be called a proof. A single discrepancy can falsify a frozen model, but a model that has been calibrated against the observable cannot claim that same discrepancy as an independent prediction.

## 4. Current verdict

- Numerical evaluation of tau_n=πlnN: PASS.
- Dimensionally complete physical equation: NOT ESTABLISHED.
- Article's use of 877.8±0.4 as current world average: OUTDATED / METHOD-SPECIFIC; PDG 2026 average is 878.3±0.4 s.
- Independent blind prediction status: NOT ESTABLISHED because N fitting includes neutron lifetime in the full fit and most LOO folds.
- Neutron lifetime: retain as a preregistered high-value held-out falsification observable.