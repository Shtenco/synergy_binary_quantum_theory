# BQG node1: точный G20 и NO-GO глобальной оценки остатка — 2026-11-05

**Фактический источник:** успешно завершённые независимые GitHub Actions runs 42007105732 (четыре настоящих K5 KA outputs), 42073256951 (SHA-верификация всех G20), 42073665150 (точная норма остатка), 42073876985 (невзвешенные столбцы оператора и нижняя граница его нормы). Нет синтетического HDA/TT.

## 1. Фактические вычисленные результаты

- Исходный frozen KB: source node1, ordered triple (0,2,3), auxiliary index channel 000, column 0, Jmax2=7; input SHA256 f486e2250d4a5678a34438467b9e338d730ebf6e1d1ec956d2de7f8c2521b293.
- SHA-verified 20 of 120 K5 covariant input states. Output support union: 1052, G20 matrix rank20, norm of coherent output 0.06524151188803094; quadrature norm0.06933112886898857; signed interference squared norm -0.0005503505568122168.
- G20 spectrum positive: min4.1695033622212527e-7, max0.0014839366479810132. All 20 mapped weighted source directions independent. This is *not* 20 physical propagating modes.
- Input square-norm fraction processed 0.5105542994964184, remaining 0.4894457005035816; full KB norm0.12845640035659833; exact remaining input norm 0.08986860634523096 (100 basis inputs).

## 2. Rigorous conditional residual upper-bound formula

For A=C_a(K_sine), P onto the 20 processed orthonormal K5 KB basis sources, Q=I-P and frozen input ψ,
```
δ_out=A Q ψ
||δ_out|| ≤ ||A Q|| ||Q ψ||
        =0.08986860634523096 ||A Q||.
```
To certify δ_out ≤10% of measured ||APψ||=0.06524151188803094, need ||AQ||≤0.07259655461597475; for 5% need ≤0.036298277307987376.

A loose factorized estimate ||A|| ≤ 4||V|| ||HE_sine|| would be valid ONLY for precisely bounded shared-domain V, H_E, and contractive holonomy/charge/intertwiner/projections; current CK primitive charge-leak max=1.0 has not been fully resolved by global post-sum check. No unconditional upper norm claimed.

## 3. New numerically proved NO-GO for global-norm certificate

Divide the *real measured* 20 outgoing columns by the **original complex coefficients of their 20 KB basis vectors**, then form Hermitian unweighted M20†M20. The largest eigenvalue's square root gives an exact finite-subspace norm:
```
||A|| ≥ ||A P20|| = 1.0773367231979318
max_j ||A e_j|| = 0.9342608558335371.
```
Since 1.0773367231979318 > 0.07259655461597475, it is mathematically impossible to certify 10%-precision from ANY global upper bound ||A|| based on the generic inequality, at this stage. Such an upper bound must be ≥1.0773, producing right-hand side ≥0.09681874987833598, i.e. ≥148.40% of known partial norm. This is a NO-GO for that **particular proof method**, NOT a lower bound on actual unknown residual ||AQψ|| and NOT a no-go for BQG.

## 4. Automatically selected next computation

Replace global bound by **Q-restricted** estimate or directly compute coherent A Q ψ. Reuse existing cached SU2 primitives but make all subsequent 100 source columns real and SHA-verifiable; group sources by actual cost, save complete outgoing vectors per key. Independently measure charge conservation after full HE summation, not only max primitive-branch leak. If a certifiable Q-restricted norm ≤0.07259655 cannot be established, complete all 100 actual outputs; avoid advertising fractional input coverage as an error guarantee. Beyond this one channel, still seven auxiliary channels, all oriented terms/adjoints, two-node covariance, ten independent off-shell HDA pairs and physical TT pole.

Source code: scripts/bqg_node1_ka_gram20_from_logs.py, scripts/bqg_node1_residual_operator_bound_gate.py, scripts/bqg20_ck_norm_lower_gate.py. Main branch unchanged; physicalization gates remain OPEN.
