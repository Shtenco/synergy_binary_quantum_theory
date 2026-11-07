# BQG v19 — independent W5 Gershgorin bound and PL16 ten-channel HDA rank obstruction

Date: 2026-11-07. Scope: independently **computed postprocessing of saved real SU(2) numerical evidence**, no new physical Lorentzian amplitudes. Neither full C1, C2 nor C3 is closed.

## 1. Genuine K5 W5 complex Euclidean Gram, independent positivity certificate

Read from the genuine W5 matrix at `research/bqg-node1-lorentzian-columns-20261031/research/node1_fullstate_gate_2026_10_31/V04_CROSS_BRANCH_REAL_EUCLIDEAN_W5_LOG_EVIDENCE.json`, Git blob SHA `0231c54d2e9a7fc9401995c999d3989f3f1f7648`, upstream GitHub Actions run `42121958625`. The upstream calculation completed and artifact upload failed quota; the log-based matrix is preserved in Git.

Gershgorin discs for the 5x5 Hermitian complex Gram G=A†A:

| row | diagonal | sum abs(offdiagonal) | disc lower |
| --- | ---: | ---: | ---: |
| 0 | 34.040571789990 | 3.054557084101 | 30.986014705889 |
| 1 | 33.260658846437 | 1.241219698185 | 32.019439148252 |
| 2 | 33.665423938263 | 2.659392028544 | 31.006031909719 |
| 3 | 33.499459067520 | 0.898496026483 | 32.600963041037 |
| 4 | 35.414459019687 | 2.279583736927 | 33.134875282760 |

Therefore for this printed Hermitian matrix `lambda_min(G) >= 30.986014705888927 > 0`, `sigma_min(A) >= 5.566508304663609`. Independent diagonalization reproduces archived `lambda_min=32.41356449685593`, `sigma_min(A)=5.693291183213444`; spectral condition number `1.0670358406952927`.

For `H_v = a E_v + Q_v`, `a=-2/3`, `Q_v=-32/9 S_v+R_v`, with `A:u -> (E_0 u,...,E_4 u)` on **the same K5 W5 and common output basis**, any nonzero vector in `ker(a A + Q)` needs `||Q||_op >= |a| sigma_min(A)=3.795527455475629`. The Gershgorin-only necessary bound is `3.7110055364424057`.

Stronger nodewise necessary conditions (EVERY one must hold to allow a common master-kernel zero) from genuine archived node spectra:

`||Q_0||>=1.4744574358967717;`
`||Q_1||>=1.658063169176699;`
`||Q_2||>=1.7243767365832534;`
`||Q_3||>=1.7681182322608326;`
`||Q_4||>=1.768118232260833.`

These are **necessary only**, not sufficiency, and no actual Q_v amplitude has been produced. Gershgorin is rigorous for the floating-point entries treated as exact numbers; rigorous SU(2) error enclosures still needed.

## 2. Separate PL16 finite E-only ten-commutator HDA rank theorem

From saved `v13_real_results.json` and `BQG_V12_DEEP_REPORT_RU.md`: ten linearly independent complex full-space `C_r=[E_0,E_r]Omega`, active `r=[1,2,3,4,5,6,8,9,10,12]`; measured `rank(C)=10`, `sigma_min(C)=1.324189319007516`, `||C||_F=7.608464409020511`. This PL16 habitat is **NOT THE K5 W5 HABITAT** above. No matrices are concatenated across the two graphs.

Let `a=-2/3`. For each lapse pair `(delta_0,delta_r)`, let `K` be the aggregate of ALL non-EE commutator channels and the independently defined D target, with exactly the same output space and conventions. The full HDA matrix restricted to these ten probes is `R=a² C+K`.

- Exact closure `R=0` forces `K=-a² C`, hence `rank(K)=10`, `||K||_F=a²||C||_F=3.3815397373424494` and `sigma_min(K)=a² sigma_min(C)=0.5885285862255626`.
- If `rank(K)<=9`, then Eckart–Young yields `||R||_op>=a² sigma_min(C)=0.5885285862255626`.
- If `||K||_op<0.5885285862255626`, then reverse triangle similarly forbids exact closure.
- One-hop (graph distance 2) six measured commutators have `a²||C_{distance2}||_F=2.07535394654993`: required aggregate compensation norm for exact zero, not a measured anomaly.

**Actual K unknown**: in particular no full Lorentzian S, R or independent quantum D on this PL16 habitat has been computed. Thus this is a sharp *conditional finite-family necessary gate*, not a no-go for BQG gravity.

## 3. GitHub Actions infrastructure admission blocker

Audited latest 100 workflow runs: 2026-11-07 has 21/21 `failure` across many independent workflows; 2026-11-06 in the same page had 31 `success`, 27 `failure`, 21 `cancelled`. Minimal V18 workflow run `42371971150` failed after ~4 s with `steps=[]`, no artifacts, and no job logs. This reproduces failure without checkout, dependencies or BQG. **Exact cause not proven**. Upstream prior W5 experiment reported artifact-upload quota exhaustion; inspect account GitHub Actions usage/budget/storage policy and workflow-run check annotations; if blocked pre-runner, contact GitHub Support with IDs `42371971150`, `128401781450`.

Do NOT attribute zero-step failures to null physical Lorentzian tensors; do NOT retry 32 costly V4 paths before runner admission is demonstrably fixed.

## 4. Reproducible local bundle

`BQG_V19_HDA_GATES.zip` includes `SOURCE_SNAPSHOT.json` (full numerical 5x5 W5 complex Gram, original blob SHA), `v19_science_audit.py`, `RESULTS.json`, `test_v19.py` (12/12 PASS) and legacy V18 validators. No new SU(2) amplitudes. The packaging step verifies real saved numbers, valid eigenvalues and Hermitian structure. Synthetic rank matrices appear ONLY in unit tests of the analytic rank lemma.

Follow-up: revive admitted Actions runner; complete 32/32 **genuine** final V4 paths and 4 independent original V2; then construct consistent Q and D and evaluate the two numerical gates above on their own matching Hilbert spaces.
