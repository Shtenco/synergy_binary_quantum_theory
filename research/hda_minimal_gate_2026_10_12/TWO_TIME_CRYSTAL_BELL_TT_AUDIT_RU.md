# BQG: два времени + кристалл времени + Bell + TT — единый проход (2026-10-12)

**Статус: математические положительные/отрицательные контроли выполнены; физическая БКГ НЕ закрыта.**

Источники: [главный 2T frontier](../../README_2T_FRONTIER.md); [2T scanner](../../BQG_2T_CLOSURE_SCAN.md); [Gauss-синглеты](../../BQG_MINIMAL_ENTANGLEMENT_VOLUME_RESULT.md); [C8 реляционная история](../../Q2_RELATIONAL_HISTORY_PROJECTOR.md); [Bell-reference](../../BELL_CHSH_Q2_REFERENCE_CONTROL.md); в разошедшихся ветках — [ориентированная фаза](https://github.com/Shtenco/binary_quantum_theory/blob/research/q2-history-fourier-complex-structure/Q2_HISTORY_FOURIER_REAL_COMPLEX_STRUCTURE.md) и [carry history C4→C8](https://github.com/Shtenco/binary_quantum_theory/blob/research/q2-history-phase-refinement-winding/Q2_HISTORY_PHASE_REFINEMENT_AND_WINDING.md). Ни одна не объявляет доказанными физические два времени или устойчивый Floquet DTC.

На двух Gauss-логических узлах с J²=-I задаём R1=J⊗I, R2=I⊗J, два C_N часовых регистра и *externally chosen* commuting constraints G1=S1⊗R1, G2=S2⊗R2. Для N=8,16,32 строим точную изометрию V_N=(1/N)Σ_(s,t)|s,t>⊗R1^sR2^t; G1V_N=G2V_N=V_N; fixed sector rank4. Два набора dressed наблюдаемых коммутируют с G1/G2 и между собой, в одной Bell seed |Φ+> дают CHSH=2√2. Но это **не** общий полный микроскопический H_E+λH_L constraint, а engineered positive control.

## Новая точная развилка двух часов
(R1R2)|Φ+>=|Φ+>, R1|Φ+>=−R2|Φ+>. Conditional pure orbit имеет ранг2, зависит от разности s−t. Для непрерывных углов тождественно:
F_Q=π²[[1,−1],[−1,1]], spec(F_Q)={0,2π²}.
Синхронный сдвиг часов информационно нулевой; независимого второго физического времени в данной подготовленной модели не получено.

## Floquet π и ложный TT-полюс
Те же Gauss shape-источники Z_A,X_A имеют Heisenberg C_n=(-1)^n I2: C0=I, C1=−I, C2=I, Hankel determinant/Schur defect =0. **C1 отрицательно**, так что это NOT 0≤T_E≤I Euclidean TT correlation and no massless graviton. Периодическое изменение только унитарное: ε=0.08 в π-пульсе даёт C8=0.802096, C16=0.286715, C32=−0.835589, т.е. не наблюдается устойчивой subharmonic rigidity. Требуется interacting many-body, stability, order-parameter correlator и physical time. Тест НЕ является настоящим DTC.

## Ретрансляция старой q2-фазовой истории
Independent Z4×Z2 update has order 4, не Z8. Coupled carry (k,0)→(k,1),(k,1)→(k+1,0) has order8. Исходная ветка доказывает математическое первое уточнение; физическая бесконечная башня carry остаётся conditional, circle phase/ winding ≠ literal timelike direction/ Floquet time crystal.

## Bars 4+2 reference / no-go finite Hilbert
В *отдельном*, вручную определённом классическом parent 4+2, Q11=X²,Q12=X·P,Q22=P² дают {Q11,Q12}=2Q11, {Q12,Q22}=2Q22, {Q11,Q22}=4Q12. SymPy PASS. Это НЕ извлечённые BQG ограничения.

Для эрмитовых конечномерных A,B,C с [A,B]=iC, [B,C]=−iA,[C,A]=iB, цикличность следа даёт Tr(A²)=−Tr(C²). Все три матрицы нулевые. Следовательно, exact nontrivial *unitary positive-Hilbert* sp(2,R) невозможно найти как простую 2×2/finite matrix triple. Нет запрета на соответствующие неограниченные бесконечномерные представления или полуклассическую Poisson closure. Отличать родительскую 4+2 сигнатуру для физического 3+1 shadow от неоднозначной записи 3+2 в старом README.

## Следующий расчёт автоматически выбран
**Не добавлять новые независимые часы и не выдавать π/Floquet за физическое ω.** Заморозить complete Hermitian Lorentzian BQG, вычислить images all nodes на 2 соседних коллективных q4 метрических источниках; одновременно independent D (10-парный HDA), projector/rigging с полным Q leakage, одетые коммутирующие block observables. Лишь после положительного physical-time construction вычислять C0,C1,C2, спектральную положительность/положительный massless TT pole. Для физической clock-ветви проверять microscopic graph-change carry selection (Z4×Z2 не даёт Z8) и отдельно устойчивость коллективного Floquet периода под возмущениями. Sp(2,R) искать по actual constraints в controlled refinement, не копировать Bars prototype.

**Запуск на клоне**: python research/hda_minimal_gate_2026_10_12/two_clock_time_crystal_reference_gate.py. Полный matrix-free N=8,16,32 рабочий расчёт, JSON, исследовательский отчёт включены в выдаваемый ZIP-пакет диалога. Physicalization registry не изменялся.
