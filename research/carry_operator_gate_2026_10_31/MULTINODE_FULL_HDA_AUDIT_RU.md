# Аудит многоузловой BQG / Full HDA — 31 октября 2026

**Итог:** настоящие многоузловые результаты и исторический `Full HDA` есть, но полного микроскопического лоренцева off-shell HDA для актуального эрмитова оператора нет. `main` HEAD = `90a5acde961c09cde565f942c9e21a7d7db1376b`, 613 blob entries; все 9 обязательных физических ворот `open_physical`.

## Уже выполнено и не должно пересчитываться

- [3-node graph-changing Euclidean HDA](../../THREE_NODE_GRAPH_HDA_RESULT.md): 3 пары, `Delta_joint(1/64)=0.0252238079`, exponent 1.0064429; не полная лоренцева HDA.
- [Oriented K5 ten pairs](../../K5_ORIENTED_QUANTUM_HDA_RESULT.md): все 10 пар, норма 3906.62207, ориентированная S5-ковариантность, но `Jmax=1/2` ниже safe wall `5/2`; status finite diagnostic, FAIL as physical HDA.
- [16-cell staggered orientation assembly](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/LORENTZIAN_16CELL_GLOBAL_ASSEMBLY.md): условная сборка 16 одноузловых Y, не полные 16-cell outgoing graph-changing columns.
- [Multi-neighbor raw Lorentzian Walsh](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/verification_results/PETER_WEYL_LORENTZIAN_ENVTRACE_WALSH_NODE012.json): настоящие 24 oriented source0 terms, four diagonal environments K1,K2 with K3=K4=0, leakage < 7e-16; *not* offdiagonal environment flips.
- [Serial Euclidean multiplicity](../../BQG_SERIAL_MULTIPLICITY_MASTER_RESULT.md): five tested j=2,5/2,3,7/2,4 positive master channel gaps, not global full Lorentzian physical projector.

## Почему «Full HDA PASS» не равно доказательству

[Исторический сертификат](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/FULL_OPERATOR_FIRST_HDA_CERTIFICATE.md) имеет `passed:true` на **условной WKB/lapse/finite-cutoff теореме**: N=Nbar+eps*n; Omega~eps^-1; Dnorm~eps^-1; смешанный относительный дефект O(eps), геометрический O(eps²). Полные микроскопические EL/LE/LL workers прерваны шестичасовым лимитом без итоговых амплитуд. [Октябрьский независимый аудит](../../research/audit_2026_10_10/FIVE_ROUTES_AUDIT.md) на fixed physical momentum наблюдал cross/D~0.266767 без убывания на проверенной последовательности.

Ещё позже найдено, что `L_raw` не всегда антиэрмитовен. Обязательный кандидат [Hermitian completion](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/LORENTZIAN_HERMITIAN_COMPLETION.md) фиксирует `S=-i/2(Lraw-Lraw†)`, `G=-2E/3-32S/9`. Старый certificate script использовал `+32i Lraw/9` и поэтому **не является прямым проверенным результатом нового физического оператора**. Актуальная [пререгистрация](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/PETER_WEYL_HERMITIAN_FULL_HDA_PREREGISTRATION.md) требует `[G0,G1]=(4/9)EE+(64/27)(ES+SE)+(1024/81)SS`, общие выходы, независимый D, пять epsilon и cutoff leakage. [Коллективный HDA intertwining](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/COLLECTIVE_HDA_INTERTWINING_THEOREM.md) даёт точную условную bound по leakage, НЕ измеренную production collective HDA.

Большая [16-cell V3 orbital reduction](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/PL_16CELL_HERMITIAN_LORENTZIAN_PREREGISTRATION_V3_ORBIT.md) разрешает 12 held-out direct workers + 36 точных симметрийных реконструкций вместо 48 тяжёлых terms; сам документ датирован до первого полного исправленного science column.

## Новый actual-data 3-node replay (внешний воспроизводимый пакет, без новой micro-amplitude fit)

Для `Hcorr,node0 = -32/9 ImPauli(Lraw)`, K3=K4=0:
```
-1.1943161007640883 Y0
+0.02499063901326094 Y0Z1
-0.008313353267907534 Y0Z2
-0.016626706535814353 Y0Z1Z2
```
Neighbor amplitudes Z1,Z2:
```
(+,+) -1.1942655215545492
(+,-) -1.1443854019471054
(-,+) -1.2109933865094424
(-,-) -1.2276200930452563
```
The 8x8 compressed Hermitian operator has rank 8, each |amplitude| appears twice as singular value, sigma_min 1.1443854019471054, condition 1.0727330940752406. Setting Z2=+1 exactly recovers earlier 4x4 two-node result, an independent replay control. This is **only a compressed source0 lower bound**; full source1 and off-code Q output missing.

## Новый расчёт, выбранный аудитом

**Начать с одного held-out source1 Lorentzian term** at same K5 face orientation and Peter-Weyl recoupling; test node-swap covariance before reconstructing 24 terms. If the transport test fails, compute source1 24-term directly (no post-hoc sign). Then obtain node0 and node1 OFFDIAGONAL neighbor-flip elements and full outgoing 4+4 columns (forward + adjoint), and only then run latest ES/SE/SS 10-pair off-shell HDA with independent D, physical rigging/history and interacting TT. Do not rerun 16-cell 48 direct terms until this small gate is clear. No physicalization ledger promotion.

Scope: repo branches/critical code, stored JSON and GitHub Actions run statuses audited; the full historical high-cost amplitude workers were NOT re-executed.