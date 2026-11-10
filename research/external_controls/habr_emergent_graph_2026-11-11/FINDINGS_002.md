# Findings 002 — provenance and constant-table reproduction

## 1. Three incompatible values of the global node count are in circulation

The target's public artifacts do not use one frozen N:

| Artifact | N | ln N | Source |
|---|---:|---:|---|
| Geometric resonance equation at K=6 | 4.4759473527×10^121 | 280.1115142805 | Direct solution of displayed R=1 |
| Stress-test script's N_phys | 4.1846576952×10^121 | 280.0442211596 | stress_test.py main(): N_math=4.475947×10^121, then subtract πK/ln(N_math) from ln(N_math) |
| Article table's N_phys | 4.185×10^121 | actual ln N=280.0443029562 | Article prints lnN≈280.047, which does not equal ln(4.185e121) |
| mass_relation.json | 4.1977×10^121 | 280.0473330085 | JSON constants block |
| Variation script | 4.2×10^121 | 280.0478807776 | variation_principle.py |

These are not rounding-equivalent when formulas contain high powers of lnN and N^(1/3). Every result must specify which N it uses. The article table also prints N^(1/3)=1.61×10^40, whereas the cube root of 4.2×10^121 is 3.4760×10^40. The printed value appears to omit the factor 10^(1/3).

## 2. Recalculation of the Python formula dictionary

The current physics_constants/stress_test.py contains 82 target entries in constants, 83 entries returned by formulas, and one returned quantity B_s_meson without a matching target entry.

At the mass_relation.json value N=4.1977×10^121, independent transcription gives:
- l_P = 1.6173305226×10^-35;
- t_P = 5.3989008875×10^-44;
- m_P = 2.1739375451×10^-8;
- c = 2.9956662593×10^8;
- G = 6.6763419975×10^-11;
- alpha = 0.00729807124271;
- m_e = 9.0793802302×10^-31;
- m_s = 1.6777595530×10^-28.

The stress_test.py main() does not use the JSON N; it computes N0≈4.1846576952×10^121. Re-evaluating the source formulas at this script N (using each formula's explicit powers of N and lnN) yields a preliminary mean relative error of about 0.2811% across 82 target entries, median about 0.2210%, with 80/82 below 1%, 68/82 below 0.5%, and 26/82 below 0.1%. These are preliminary independent calculations, not yet a fully validated reproduction of the script output. They differ from the article's reported 0.26%, 79/82, 62/82, and 24/82. A line-for-line executable reproduction remains required.

## 3. Strange-quark table contradiction

The article's displayed particle table gives m_s=K(lnN)^5/(sqrt(pi) N^(1/3)) and claims m_s=9.64×10^-30 kg. At K=6, N=4.1977×10^121, the displayed formula evaluates to 1.6777595530×10^-28 in the article's numerical convention. The article's printed value is smaller by a factor of about 17.40. The current code's target dictionary instead uses 1.66×10^-28, broadly consistent with the formula. Thus the article's table formula and printed result are internally inconsistent, while the code uses a different reference value.

## 4. Reference-value hygiene

The code's ep_0 target is 8.8725366415×10^-12, while the CODATA/SI vacuum permittivity corresponding to the listed mu_0 and exact c is about 8.8541878×10^-12. The difference is about 0.207%. A target constant must not be accepted merely because it appears in the author's dictionary; provenance and uncertainty must be checked independently.

## 5. Statistical and fitting protocol issues confirmed from source

- find_optimal_N() minimizes the sum of squared log-ratios against the full constants dictionary. It estimates N from the same observables subsequently described as predictions.
- leave_one_out_test() excludes one target at a time but refits N using all remaining targets. This is a leave-one-observable-out calibration exercise, not a blind prediction from an independently derived N.
- scipy minimize_scalar(method='brent', bracket=(a,b)) uses a bracketing interval, not a hard bounded interval. The comment ±10% does not constrain the optimizer to that range.
- The final status is based on a threshold count of errors below 1% and labels the model as derived from first principles; that conclusion does not follow from the fitting procedure.
- Physical dependencies such as h=2πℏ and Planck units defined from c, G, ℏ are counted as separate entries unless explicitly grouped.

## 6. Current verdict

- One fixed global N across article, scripts, and JSON: FAIL.
- Strange-quark displayed formula/value consistency: FAIL.
- Formula dictionary count vs target dictionary count: 83 vs 82; one unmatched output.
- Reported aggregate accuracy: NOT REPRODUCED YET; preliminary independent metrics differ.
- Blind prediction / no-fit status: FAIL AS CLAIMED; the published code explicitly refits N.
- Remaining constant formulas: AUDIT IN PROGRESS, not yet individually passed.