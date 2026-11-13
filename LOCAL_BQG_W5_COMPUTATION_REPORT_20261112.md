# Локальный пересчёт BQG/W5: branch cancellation, exact fixtures и физические пределы

Дата: 2026-11-12. Все результаты ниже получены локальным исполнением файлов из `/mnt/data`, без использования GitHub Actions как вычислительного backend.

## 1. Выполненные локальные исполнения

Синтаксически скомпилированы через Python `compile(source, filename, "exec")` 10 модулей:
- `w5_all_node_carrier.py`
- `w5_all_node_operator_adapter.py`
- `w5_exact_closure_audit.py`
- `w5_operator_audit.py`
- `w5_source_gap_audit.py`
- `w5_branch_matrix_exporter.py`
- `w5_v10_dual_threshold_runner.py`
- `w5_v10_json_to_matrix.py`
- `bqg_physical_history_decoherence_engine.py`
- `bqg_interference_causality_bridge.py`

Локально выполнены:
- five-node carrier kinematic self-test: `ALL_NODE_CARRIER_KINEMATIC_SELFTEST_PASS`;
- all-node operator adapter self-test: `PASS`;
- physical history/decoherence reference engine: `PASS_DECOHERENCE_FUNCTIONAL_REFERENCE`;
- interference/causality bridge: `PASS_REFERENCE_IDENTITIES`;
- exact SymPy full-covariant gate regression: все заявленные regression controls прошли;
- exact 36x36 fixture audit и повторный 8x8 raw-operator audit.

Системное сообщение spreadsheet-runtime warmup в stderr возникало при старте Python, но не меняло exit code у перечисленных успешных расчётов. Это не ошибка самих расчётных скриптов.

## 2. 36x36 all-J/M fixture: независимый пересчёт

Взяты байты из `w5_audit_toolkit_v2.zip`:
- `w5_exact_allJ_multinode_selftest.npz`;
- `w5_exact_allJ_multinode_basis.json`;
- `w5_exact_selftest_ledger.json`.

Размерности: `M,D = 36x36`, `branches = 5x36x36`.

Численно:
- `||M-D†||_F = 0`;
- `||M-M†||_F = 0`;
- `M=D=I_36`;
- все 36 собственных значений M равны 1; ранг 36;
- weighted branch sum совпадает с M с Frobenius residual 0;
- все пять отдельно сохранённых adjoint branches совпадают с независимо экспортированными adjoint matrices с residual 0.

На каждом из двух локальных узлов supplied basis содержит J=0, 1/2, 1; это 36 совместных двухузловых состояний. Это не пять узлов полного W5.

### 2.1 Branch-resolved SU(2): когерентное восстановление

Локально построены генераторы Jz, J+, J− из supplied basis labels (J,M,multiplicity). Проверены коммутаторы на обоих узлах.

Для полной матрицы M=I_36 все шесть остатков (3 генератора × 2 узла) строго равны нулю.

По ветвям:
- branch 0: все шесть коммутаторов равны нулю;
- branches 1–4: Jz-коммутаторы равны нулю, но J+ и J− дают Frobenius residual `5.477225575051661` на каждом узле;
- `branch_1 + branch_2 = 0` точно;
- `branch_3 + branch_4 = 0` точно;
- следовательно, `branch_1 + ... + branch_4 = 0`, а полный weighted sum равен branch 0 и M.

Это демонстрация реального механизма coherent cancellation в предоставленной матричной fixture. Однако fixture содержит пять ветвей и два локальных узла; она не заменяет экспорт 24 epsilon-ветвей × 5 узлов из production W5.

### 2.2 Closure: исправление интерпретации

Повторный `w5_exact_closure_audit.py` с единственным seed 0 даёт reached=1 из 36, max depth=0. Причина: M=D=identity, поэтому отдельный seed не переходит в другие состояния.

Запуск со всеми 36 seeds даёт 36/36, но это тривиально, поскольку все basis IDs уже были начальными seeds. Поэтому прежнее значение 36/36 не является доказательством нетривиального динамического замыкания. Оно только подтверждает отсутствие дополнительных переходов в этой identity fixture.

## 3. 8x8 raw non-Hermitian fixture: повторный audit

Независимо пересчитан `w5_audit_control.npz`:
- `epsilon_adjoint = 0`;
- `epsilon_H = 0.39816652969059774`;
- `rank(M)=8`;
- `||M||_F = 3.622154055254967`;
- `||H||_F = 3.5496478698597698`, где H=(M+M†)/2;
- `||A||_F = 0.7211102550927979`, где A=(M−M†)/2;
- Pythagorean decomposition residual `1.2212453270876722e-15`.

Собственные значения Hermitian part H в этой конкретной fixture положительны: приблизительно
`[0.7, 0.7, 0.7, 0.96849271, 0.96849271, 1.0, 2.03150729, 2.03150729]`.

Это не делает M эрмитовым и не доказывает положительность физического master constraint BQG. Adjoint consistency и raw Hermiticity остаются разными требованиями.

### 3.1 Локальный ledger

- stage-local intermediate-K discarded norm: `1.1180339887498949e-12`;
- относительно `||M||_F`: `3.086654989529967e-13`;
- spin-wall coherent L2: `1.1180339887498949e-9`;
- относительно `||M||_F`: `3.0866549895299675e-10`.

Это stage-local диагностические величины, не глобальная ошибка оператора. Нужна downstream propagation отброшенных вкладов с coherent aggregation по destination sectors.

## 4. Exact SymPy full-covariant regression

Локально распакован `W5_FULL_COVARIANT_EXACT_GATE_20261112.zip` и запущены:
- `w5_full_covariant_gate_regression.py`;
- `w5_full_covariant_gate.py` на supplied exact input.

Проверки:
- exact base case PASS;
- branchwise non-scalar / coherent scalar distinction PASS;
- perturbation `sqrt(2)*1e-30` обнаружена exact arithmetic;
- branch-local adjoint defect отвергнут;
- adjoint consistency отделена от raw Hermiticity;
- внешний state за пределами supplied carrier найден, и global closure withheld.

Это exact algebraic regression на экспортированной fixture, не выведение матриц из production W5 source.

## 5. Физические reference checks

Из локальных запусков:
- decoherence functional Hermiticity/strong positivity reference tests прошли;
- exhaustive class-operator normalization прошла;
- environment records подавляют off-diagonal cross terms в toy factorization;
- Sorkin I3 = `1.1102230246251565e-16`;
- interference-aware history depth: support depth 1, coherent depth 2 после cancellation;
- synthetic threshold fronts: photon=matter=tensor=2, spread 0.

Все эти значения относятся к supplied finite reference inputs. Они не доказывают физическую BQG factorization, Born rule, физический clock или continuum causal cone.

## 6. Production W5: жёсткая локальная проверка наличия исходников

Выполнен прямой вызов:
```bash
python /mnt/data/w5_all_node_operator_adapter.py --scripts-dir /mnt/data --real-su2
```
Он завершился `ModuleNotFoundError: k5_peter_weyl_safe_hda_column`.

Проверены локальные пути шести обязательных production-модулей:
- `k5_peter_weyl_safe_hda_column.py`;
- `peter_weyl_covariant_K_leg_gate.py`;
- `peter_weyl_covariant_volume_leg_gate.py`;
- `peter_weyl_covariant_composition_gate.py`;
- `peter_weyl_covariant_K_composition_gate.py`;
- `peter_weyl_lorentzian_ordered_triple_gate.py`.

Все шесть отсутствуют в `/mnt/data`. Следовательно, production `C(V) -> C(K) -> C(K)` amplitudes, полный all-node/all-J,M production commutators и production repeated closure остаются NOT RUN. Локальный archive — это audit toolkit и fixtures, а не полный production source checkout.

## 7. Точный следующий вычислительный шаг

Необходимо локально получить полный production `scripts/` source tree, после чего:
1. экспортировать каждую из 24 epsilon-ветвей и direct-adjoint независимо для всех пяти узлов;
2. построить полный all-node/all-J,M carrier с multiplicity labels;
3. прогнать реальные `C(V)->C(K)->C(K)` пути, сохраняя intermediate charge sectors;
4. измерить coherent discarded norm после каждого реального threshold/projection/spin-wall stage;
5. проверить branchwise и coherent SU(2) commutators на всех колонках, разрешая output support вне input carrier;
6. повторить closure с одним seed per irreducible multiplet и отдельно со всеми seeds; любой cap означает INCOMPLETE.

Никакой production W5 или full-BQG physical PASS не объявляется до этих шагов.
