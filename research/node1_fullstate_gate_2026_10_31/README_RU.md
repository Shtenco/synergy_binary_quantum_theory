# BQG: соседний узел 1 — подготовленный реальный fullstate gate (31.10.2026)

## Проверено

Архивный source0 двухсоседний diagonal environment и асимптотический **условный** Full HDA разобраны ранее в `MULTINODE_FULL_HDA_AUDIT_RU.md`; восстановлен 8×8 compression, но **полного физического полного HDA нет**.

Устаревший `scripts/peter_weyl_lorentzian_fullstate_triple_worker.py` в ветке `research/bcqg-core-candidate-v1` жестко выбирал `PW.basis_full_jhalf()[0]`, поэтому повторный запуск узлов 0/1 сам по себе не вычислял четырех исходных двухузловых колонок. Исправленный script в этой исследовательской ветке принимает `--column-index` с диапазоном [0,31], записывает `input_Ks`, `column_index` и `source_pair_K01` в JSON для последующего аудита. Реальный contraction `LP.ordered_triple_state`, basis, `Jmax2=7`, covariance/Gauss closure и thresholds **не изменены**.

Доказанная аналитически и перепроверенная в локальном Python таблица исходных 4 K-пар:

| column_index | K0 | K1 | K2,K3,K4 |
|---:|---:|---:|---|
| 0 | 0 | 0 | 0,0,0 |
| 8 | 0 | 2 | 0,0,0 |
| 16 | 2 | 0 | 0,0,0 |
| 24 | 2 | 2 | 0,0,0 |

Перестановка source0↔source1 при том же `column_index` **не** гарантирует равенство комплексных amplitudes. В локальных intertwiner conventions знаки исходящих/входящих линков и SU(2) recoupling меняются; до прямой held-out node1 проверки нельзя назначать относительный физический знак Y0±Y1. Простая перестановка чисел без unitary intertwiner map была бы ошибкой.

`scripts/bqg_lorentzian_2node_column_preflight.py` предрегистрирует exact `2 nodes ×4 columns×24 ordered epsilon triples=192` **forward raw** science packets, checks 24 distinct orientation tags per source (+12/-12), column ID, K-list, source ID, epsilon, `passed`, both JSON and NPZ, returns EXIT2 BLOCKED if data missing. Это не готовые результаты и не обещание их получения в этой сессии. Для действительного исправленного эрмитова `S=-i/2(Lraw-Lraw†)` нужны также independent adjoint-channel witnesses либо полные outgoing matrix elements. **192 forward packets не составляют автоматически полное Hermitian (S_v) на общей области.**

## Следующий самый дешёвый микроскопический расчёт

1. Один direct node1 held-out ordered triple (например, первые K=(0,0,0,0,0), source1, triple 0,2,3; epsilon из независимого preflight), плюс соответствующий source0 transported reference. Ввести явный SU(2) oriented node-swap intertwiner `U_swap` на полных Gauss states и проверить residual `||L1 U_swap - eta U_swap L0||` на измеренном столбце, отдельно outer leakage.
2. Если residual проходит frozen tolerance, раскрыть covariance по full orbit; если FAIL — отказаться от этой экономии и сделать все 24 term/column непосредственно.
3. Получить недостающие full source0/source1 outgoing raw/adjoint packets с frozen support; только затем `ES/SE/SS`, independently defined `D` 10 node pairs on common habitat, projector/rigging, TT pole.
4. Отдельно контролировать representative-state Jmax=7/2 и second hit Jmax≥13/2 согласно archived Lorentzian HH support bound.

## Репродукция команд

```bash
python scripts/peter_weyl_lorentzian_fullstate_triple_worker.py \
  --a 0 --b 2 --c 3 --coefficient -1 --source 1 \
  --column-index 0 --json-output /tmp/bqg_node1_T023.json \
  --state-output /tmp/bqg_node1_T023.npz

python scripts/bqg_lorentzian_2node_column_preflight.py --root artifacts/lorentzian_columns
```

Эта запись показывает **план исполнения**; первый ресурсоёмкий реальный расчёт в данном сеансе не выполнен. Локально проверены только исходные индексы [0,8,16,24], полнота 24 ε-меток и баланс знаков 12+12 по каждому source. Результат физической верификации: **BLOCKED**, а не PASS.

Ничего не изменено в основной ветке, физические gates остаются открыты.
