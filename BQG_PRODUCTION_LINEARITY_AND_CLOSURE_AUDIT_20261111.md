# BQG production linearity, closure loss and evidence ledger

**Дата аудита:** 11 ноября 2026 года  
**Репозиторий:** `Shtenco/synergy_binary_quantum_theory`  
**Ветка:** `research/lemaitre-bqg-causal-cone-20261110`  
**Проверенный HEAD до фиксации этого документа:** `4d328b9c934f97d055cd679e5c0348261e9f1523`  
**Baseline:** `7c4bd56c7d78bdba333246a7910f6161cac7d28a`  
**Назначение:** консолидировать подтверждённые изменения V09/V10/V11, новый аудит численной линейности и поддержки, ограничения интерпретации и воспроизводимый план следующих gates.

> **Статус документа:** технический аудит исходного кода и спецификация следующих проверок. Он не объявляет новые тесты запущенными, CI пройденным, физический оператор нормированным или причинный конус доказанным.

## 1. Executive evidence ledger

| Область | Подтверждённый результат | Что результат НЕ доказывает |
|---|---|---|
| V09 raw W5 column | Реализованы Euclidean A и raw Lorentzian B через K-K-V, epsilon-orbit из 24 triples, scalar projection и leakage diagnostics | Полный физический (H_L), абсолютную нормировку или физическое ядро |
| V10 adjoint closure | Прямой reversed-order adjoint строится отдельно; closure исследует forward и adjoint; matrix adjoint residual сравнивается с (B^\dagger) | Линейность всех production-функций на произвольных суперпозициях |
| V10 acceptance | Исправлен conjunctive gate: adjoint/closure PASS больше не может обходить LSINE gate, spin wall, leakage и scalar projection | Что LSINE/physical normalization уже завершены |
| V11 common carrier | Node-local V09 contributions когерентно суммируются; единый threshold применяется после coherent sum; проверяется (B=B_{sym}+B_{anti}) и V09/V10 cross-check | Что внутренние пороги KKV или closure discovery без потерь |
| Joint-Gram hardening | Проверка basis identity, явных (c_E,c_L), непустых carriers и корректное дополнение структурного правого ядра при rectangular SVD | Физический смысл коэффициентов, если они не выведены независимо |
| Causal-cone audit | Сформулированы независимые operational, retarded, spectral и universality gates | Наличие физического времени, ретардированного оператора или доказанного светового конуса |

### Коммиты, уже присутствующие в истории ветки

- `a4fb1cfdfc5a5c99415a77ccb0eba933f2ed20fc` — V11: threshold после когерентной суммы.
- `44feba50f09dbe18e03237d93e9dcaf730c8db8a` — V10: conjunctive LSINE/spin-wall/adjoint-closure gate.
- `4d328b9c934f97d055cd679e5c0348261e9f1523` — документирование acceptance contract (docstring-only follow-up).

`main` не изменялась; слияние не выполнялось.

## 2. V11: корректное место coherent threshold

Для node-local вкладов (b_v) неверно в общем случае сначала пороговать каждый узел, а затем складывать:
[
\sum_v \mathbf 1_{|b_v|>\varepsilon}b_v
\ne
\mathbf 1_{|\sum_v b_v|>\varepsilon}\sum_v b_v.
]
Контрольный случай пяти вкладов (0.6\varepsilon) даёт суммарно (3\varepsilon), хотя каждый вклад меньше порога.

Исправленный V11:
1. сохраняет сериализованные node-local contributions;
2. когерентно складывает их по общему state key;
3. применяет один общий `amp_tol` к сумме A и B;
4. сравнивает суммарный V09 B с V10 forward map;
5. проверяет разложение (B=B_{sym}+B_{anti});
6. сохраняет basis/provenance hashes.

У serializer V09 есть собственный floor (10^{-13}) на отдельный элемент. Поэтому условие `amp_tol > N_nodes * 1e-13` ограничивает worst-case сумму именно serializer-omitted node-local элементов. Оно **не** ограничивает потери от всех внутренних порогов KKV и не превращает общий V11 carrier в exact/no-pruning carrier.

Тест `test_threshold_is_applied_after_coherent_node_sum` использует FakeV09/FakeV10 и проверяет coherent node accumulation. Это полезный регрессионный тест контракта V11, но он не проверяет внутренний KKV production path на линейность.

## 3. V10: closure and adjoint acceptance

V10 вычисляет forward (B) и независимо упорядоченный direct adjoint. Для принятой конвенции K-sine anti-Hermitian, V Hermitian; два K-adjoint знака минус сокращаются. В direct adjoint auxiliary indices транспонируются явно:
[
K_a[j,i]\to K_b[k,j]\to V_c[i,k].
]

Исправленный `closure_acceptance` требует одновременно:
- forward norm (>10^{-10});
- adjoint norm (>10^{-10});
- relative adjoint residual (<10^{-8});
- ноль scalar-projection rejected terms;
- max component leakage (<10^{-9});
- max spin не выше (3.5+10^{-12});
- словарь результата LSINE и `passed is True`.

Это устраняет логическую дыру, в которой малый adjoint residual мог дать PASS при проваленном LSINE или spin-wall gate.

**Ключевое ограничение:** матрица, собранная из действий на отдельных базисных векторах, может быть корректно проверена на (B_{direct}^\dagger=B^\dagger), но это ещё не означает, что исходная production-функция на произвольной суперпозиции равна действию этой матрицы. Для такого вывода сначала нужна production-linearity gate.

### 3.1. Closure threshold и недоучтённые потери

В `run_closure` для каждой входной базисной колонки:
- суммируются node maps;
- затем (f_{total}) и (ad_{total}) фильтруются по `amp_tol`;
- только оставшаяся поддержка добавляется в closure.

`omitted_new_support_keys` учитывает только прошедшие порог новые ключи, заблокированные `max_states` или `max_depth`. Он не считает переходы, удалённые по `amp_tol`.

При этом `covariant_to_gauss_map` отдельно возвращает rejected count и L2 norm для scalar projection; V10 их накапливает. Это отдельный канал потерь и не заменяет амплитудный discard ledger.

Следовательно, текущие поля `omitted_new_support_keys` и `scalar_projection_rejected_norm` не описывают все причины потери поддержки.

**Обязательное расширение diagnostics до следующего научного PASS:**
- `below_threshold_forward_count`;
- `below_threshold_adjoint_count`;
- `below_threshold_forward_l2`;
- `below_threshold_adjoint_l2`;
- `max_discarded_forward_amplitude`;
- `max_discarded_adjoint_amplitude`;
- отдельно — внутренние pruning counters каждого production stage.

Нормы ниже порога следует считать по когерентно суммированным output map перед closure pruning. Если внутренний primitive уже удалил вклад, его необходимо учитывать отдельным счётчиком на том уровне, где удаление произошло. Нельзя восстановить такие потери только по финальному словарю.

## 4. Классификация порогов: fixed matrix sparsification vs state-dependent pruning

Не всякий порог автоматически нарушает линейность. Необходимо различать два класса.

### 4.1. Фиксированный порог матричного коэффициента

Пример: отбрасывание заранее вычисленного элемента (V_{row,index}), если его модуль меньше фиксированного порога, при неизменной маске для заданного блока.

Это задаёт приближённую фиксированную матрицу (A_\tau). Действие (x\mapsto A_\tau x) остаётся линейным. Ошибка относительно исходной матрицы должна быть измерена, но сама фиксированная маска не ломает additivity/homogeneity.

### 4.2. State-dependent amplitude pruning

Пример:
[
z\leftarrow dst[k]+scale\cdot a,qquad
\text{сохранить }z\text{ только если }|z|>\tau.
]
Здесь решение зависит от входной амплитуды, порядка накопления и промежуточного состояния словаря. Это уже не просто фиксированная sparse matrix.

Для hard threshold:
[
P_\tau(x)_i=
\begin{cases}
x_i,&|x_i|>\tau,\\
0,&|x_i|\leq\tau,
\end{cases}
]
в общем случае:
[
P_\tau(x+y)\ne P_\tau(x)+P_\tau(y),\qquad
P_\tau(\alpha x)\ne\alpha P_\tau(x).
]

Контрпример аддитивности: две отдельные компоненты (0.6\tau) удаляются, но в одной и той же координате их сумма (1.2\tau) должна сохраниться. Контрпример зависимости от порядка: вклады (1.1\tau,-0.9\tau,0.9\tau) могут дать разные результаты при разных последовательностях промежуточного удаления.

### 4.3. Подтверждённые production hotspots

### 4.3.1. Подтверждённый активный V10 call path (с привязкой к строкам)

В V10 `node_kkv_maps` передаёт `RAW.COMP.C_volume_component` и `RAW.KCOMP.C_K_component`. Для этого пути подтверждены следующие активные точки:

| Файл и строки | Операция | Классификация |
|---|---|---|
| `peter_weyl_covariant_composition_gate.py:48-54` | `add(dst,src,scale,tol)` порогует после каждого добавленного вклада и удаляет ключ, если промежуточная сумма опустилась ниже `TOL=1e-11` | State-dependent и order-dependent |
| `peter_weyl_covariant_composition_gate.py:99-110` | `inverse_complete` вызывает `CK.project_branch_complete_charges`, затем складывает результат через `add` | Ветвевые state-amplitude losses плюс incremental pruning |
| `peter_weyl_covariant_composition_gate.py:113-120` | `close_complete` вызывает `CV.project_covariant_branches([hb], source_v)`, затем `add` | Branch-amplitude pruning плюс incremental pruning |
| `peter_weyl_covariant_composition_gate.py:142-157` | `C_volume_component` складывает прямую и (hVh^{-1}) части через тот же `add` | State-dependent, order-dependent |
| `peter_weyl_covariant_K_leg_gate.py:53-60` | `add_state` с `tol=1e-11` | State-dependent, order-dependent |
| `peter_weyl_covariant_K_leg_gate.py:89-134` | `project_branch_complete_charges`: (val=amp\cdot\prod c), порог до и после сборки output | State-dependent branch pruning |
| `peter_weyl_covariant_K_leg_gate.py:168-185` | `apply_volume_complete_state`: fixed (V_{row,index}) mask (10^{-13}), затем финальный state-amplitude filter (10^{-11}) | Первый порог — fixed matrix sparsification; финальный — nonlinear state threshold |
| `peter_weyl_covariant_K_leg_gate.py:207-221` | `apply_HE_complete_key` суммирует completed branches через `add_state` | Активный путь построения basis column с промежуточным pruning |
| `peter_weyl_covariant_K_composition_gate.py:43-51` | `apply_HE_complete_state_custom` масштабирует cached basis-column contribution коэффициентом `amp` и сразу передаёт в `COMP.add` | Прямое доказательство state-dependent pruning при накоплении произвольной суперпозиции |
| `peter_weyl_covariant_K_composition_gate.py:54-64,75-99` | Коммутатор (K=[V,H_E]), (C_e(K)) и их суммы используют те же thresholded helpers | Активный K-leg production path |
| `peter_weyl_covariant_volume_leg_gate.py:98-145` | `project_covariant_branches` вычисляет (amp\cdot\prod c) и фильтрует ветви; вызывается активным `COMP.close_complete` | State-dependent branch pruning |

Особенно критична связка строк `K_composition_gate.py:43-50`: для произвольного входа `state` каждая basis column умножается на фактический коэффициент `amp`, после чего накопитель с порогом решает, сохранить или удалить сумму. Это не эквивалентно умножению входного вектора на заранее зафиксированную матрицу.

### 4.3.2. Смежные/legacy функции, не подтверждённые как часть V10 KKV route

В `k5_peter_weyl_safe_hda_column.py` найдены `add_dict` (191-195), `apply_primitive_oriented` (231+), `apply_T_cached_state` (277-282) и `prune_state` (292-293). Однако прямой поиск вызовов в активных V10 modules не обнаружил обращений к `PW.add_dict`, `PW.apply_primitive_oriented`, `PW.apply_T_cached_state` или `PW.prune_state`. Модуль `k5_peter_weyl_safe_hda_column.py` действительно предоставляет базовые геометрические функции, используемые текущим путём (например, `apply_path_branch`), но перечисленные legacy helpers нельзя записывать в список доказанных V10 pruning hotspots без отдельного caller path. Они остаются отдельной задачей аудита других HDA execution routes.

Аналогично, `apply_HE_complete_state` в `peter_weyl_covariant_K_leg_gate.py:230-236` имеет финальный amplitude filter, но текущий V10 K-composition path использует `apply_HE_complete_key` через `apply_HE_complete_state_custom`; поэтому этот wrapper не следует автоматически считать активным V10 вызовом.

Важная оговорка для column builders: порог, применённый только при построении фиксированной basis column, сам по себе не доказывает нелинейность уже зафиксированного столбца. Но active `apply_HE_complete_state_custom` порогует уже сумму basis-column contributions, умноженных на произвольные input coefficients; следовательно, риск для полной production map реален и должен быть измерен L2/L3 тестами.

Таким образом, аудит теперь отделяет **подтверждённые активные V10 точки**, **fixed matrix masks**, **смежные legacy функции** и **необходимые эмпирические проверки**. Он не подменяет L2/L3 test каждого production primitive.

## 5. Новые acceptance gates: L1 / L2 / L3

### L1 — deterministic threshold algebra

Изолированно проверить primitive accumulators:
1. (0.6\tau+0.6\tau) в одной координате;
2. все перестановки набора (1.1\tau,-0.9\tau,0.9\tau);
3. cancellation и последующее re-accumulation;
4. complex phases;
5. homogeneity для нескольких (α\in\mathbb C);
6. сравнение обычного, pairwise и compensated summation.

Проверять `COMP.add`, `CK.add_state`, `PW.add_dict` и локальные accumulator helpers отдельно. Ожидаемое поведение должно быть задано заранее: либо strict linear mode без state-amplitude pruning, либо явно маркированная нелинейная approximation policy.

### L2 — production primitive linearity

Для каждого реально вызываемого production primitive:
- (F(x+y)) против (F(x)+F(y));
- (F(\alpha x)) против (\alpha F(x));
- прямое действие на суперпозиции против суммы basis-column actions;
- permutation stability при изменении порядка входных sparse keys.

Минимальные targets: `apply_HE_complete_state`, `apply_K_complete_custom`, `C_K_component`, `C_volume_component`. Тестировать тот же monkey-patch/zero-aware environment, что используется V10/LSINE, и фиксировать cache state.

### L3 — complete KKV production map

На замороженном конечном секторе построить basis matrix (M) из действий на (e_j), затем независимо вычислить production action:
[
F\left(\sum_j c_j e_j\right)
\stackrel{?}{=}M c
=\sum_j c_jF(e_j).
]

Проверить случайные deterministic complex coefficient vectors и специально подобранные near-threshold vectors. Report:
[
\delta_{add}=
\frac{\|F(x+y)-F(x)-F(y)\|}
{\max(\|F(x+y)\|,\|F(x)\|+\|F(y)\|,\epsilon)},
]
[
\delta_{hom}=
\frac{\|F(\alpha x)-\alpha F(x)\|}
{\max(\|F(\alpha x)\|,|\alpha|\|F(x)\|,\epsilon)}.
]

Порог PASS должен быть preregistered и согласован с precision/conditioning. Не выбирать его после просмотра желаемого результата.

## 6. Рекомендуемая архитектура строгого linear mode

1. Заморозить конечные input/output bases и их content hashes.
2. Построить каждую basis column без state-dependent pruning после умножения на произвольный input coefficient.
3. Использовать deterministic pairwise или compensated complex accumulation.
4. Собрать sparse matrix (M).
5. Если нужен pruning, применить одну объявленную fixed matrix-element mask к (M), а не к произвольному состоянию на каждом промежуточном шаге.
6. Применять (M) к произвольным state vectors без amplitude pruning.
7. Получить (M^\dagger) из уже замороженной матрицы и независимо проверить direct reversed-order adjoint.
8. Отдельно измерять support leakage и matrix sparsification error.

Если точное полное closure слишком дорого, допустим finite sector с явной spin wall и измеренным leakage. Нельзя называть это exact closure без уточнения domain.

## 7. Closure convergence protocol

Минимальная последовательность:
[
\varepsilon_{amp}\in\{10^{-8},10^{-10},10^{-12},10^{-14}\}
]
плюс strict/no-state-pruning reference, если вычислительно достижим.

Для каждого запуска сохранять:
- closure size, depth histogram, число вычисленных колонок;
- state-cap/depth-cap hits;
- forward и direct-adjoint support/norm;
- below-threshold counts, L2 losses, max discarded amplitude;
- scalar projection rejected count/norm;
- component leakage;
- adjoint residual;
- rank, nullity, singular/eigen spectra;
- exact basis hash и source hashes.

Сравнивать спектры разных запусков только после embedding в один common basis или на явно определённом общем подпространстве. Разные basis dimensions и разные basis identities не дают прямого сравнения eigenvalue arrays.

Важно: уменьшение `amp_tol` не обязано давать монотонную последовательность, пока state-dependent internal filters остаются включёнными. Поэтому convergence study должен фиксировать все пороги и отделять их эффекты.

Для строго линейной цепочки операторов можно использовать оценку:
[
\|\Delta_{final}\|\leq
\sum_r\left(\prod_{s=r+1}^{n}\|A_s\|\right)\|e_r\|.
]
Но она не применяется автоматически к текущей state-dependent thresholded production map: сначала требуется линейность либо независимая линейная reference realization.

## 8. Hermitian ordering and normalization ledger

При принятой конвенции:
[
H_E^{sine}=\frac{T-T^\dagger}{2i},\qquad
K_{comm}=[V,H_E^{sine}].
]
При Hermitian (V,H_E^{sine}), (K_{comm}) anti-Hermitian.

Преобразование двух K legs:
[
K_H=\frac{1}{i\hbar}K_{comm},\qquad
B_H=-\hbar^{-2}B_{comm}.
]
Поэтому symmetric candidate:
[
B_{H,sym}=-\frac{1}{2\hbar^2}
(B_{comm}+B_{comm}^\dagger).
]
Знак минус обязателен в этой конвенции. Это algebraic conversion, не завершённая физическая нормировка.

Классическая derivative-free identity:
[
H_E^{kin}=-\beta^2Q_{DW},\quad
H_L^{corr}=(1+\beta^2)Q_{DW},\quad
H_E^{kin}+H_L^{corr}=Q_{DW}.
]
Она фиксирует относительный множитель в указанном секторе, но не абсолютный W5 coefficient. Нормировка
[
\mathcal N^{reg}_{W5}(\beta,\kappa,\hbar,\ell_P,\ell,\text{trace},\text{orientation},N)
]
остаётся открытой и не должна подгоняться под желаемое ядро.

## 9. Physical kernel and causal cone: evidence boundary

Даже устойчивое конечномерное W5 nullspace не является автоматически полным физическим constraint kernel. Необходимы:
- выведенные (c_E,c_L) и физическая нормировка;
- полный физический (H_L), ordering и domain;
- symmetric/self-adjoint extension analysis;
- physical inner product / rigging map;
- relational observables и physical clock;
- regulator/cutoff removal;
- сравнение common kernel и sum kernel.

Физический causal cone также не доказан. Нужны независимые реконструкции:
1. operational causal influence из микроскопических binary events, а не корреляция сама по себе;
2. causal partial order: acyclicity, transitivity, locality;
3. conformal metric из causal order при явных предпосылках;
4. volume density для фиксации conformal scale;
5. retarded Green function и support propagation;
6. независимое сравнение со spectral metric;
7. universality для scalar/vector/tensor/matter sectors.

Критический cross-check:
[
\mathcal C_{binary}
\stackrel{?}{=}
\mathcal C_{retarded}
\stackrel{?}{=}
\mathcal C_{spectral}.
]

Точные дисперсионные выражения reduced model
[
E_\pm(k)=\pm\sqrt{4m^2+9\cos^2k},\quad
\Delta_{sp}=2|m|,\quad
\xi^{-1}=\operatorname{arsinh}(2|m|/3),\quad
v_{front}=\sqrt{4m^2+9}-2|m|
]
относятся только к этому reduced model и не могут называться полной BQG light cone.

## 10. CI provenance snapshot

При подготовке этого аудита GitHub branch API подтвердил HEAD `4d328b9c934f97d055cd679e5c0348261e9f1523`. Последний опрос на 11 ноября 2026, около 00:10 UTC:

| Workflow run | Job(s) | Состояние | Artifacts |
|---|---|---|---|
| [42908462792](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42908462792) | `w5-common-carrier` | queued | none |
| [42908462954](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42908462954) | `adapter-tests` | queued | none |
| [42908462865](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42908462865) | Jmax2=10,12,14 | queued | none |

Эти runs привязаны к предшествующему коду commit `44feba50...`; текущий HEAD `4d328b9...` отличается docstring-only commit. До получения jobs conclusions и artifacts запрещено писать PASS/FAIL или заявлять numerical convergence.

Workflow V11 имеет `concurrency.cancel-in-progress: true`; новый push на ту же ветку способен отменить активный run. Поэтому последующие изменения следует отправлять после завершения/разблокировки этих проверок либо сознательно отменять их только при критическом исправлении.

## 11. Next actions — без пропуска gates

1. Повторно опросить все три CI runs; сохранить job conclusion, steps, logs и artifacts.
2. Если CI завершился, проверить артефакты V11: пять nodes, все пять W5 K seeds, complete closure/no caps, common basis hash, source hashes, cross-check и decomposition residuals, leakage, finite matrices.
3. Проверить V10 cutoff jobs Jmax2=10/12/14: closure size, cap hits, residuals, losses, basis comparability, rank/spectrum.
4. Реализовать L1 accumulator tests.
5. Реализовать L2 production primitive tests с фактическими вызовами, не synthetic fakes.
6. Реализовать L3 KKV matrix-reconstruction test.
7. Добавить before-threshold loss counters и только затем повторять cutoff convergence.
8. Лишь после production-linearity и closure gates обсуждать adjoint как операторное утверждение.
9. Продолжить физическую нормировку, kernel и causal-cone gates независимо друг от друга.

## 12. Final status

Доказуемый прогресс: исправлены конкретные acceptance/threshold defects; V09/V10/V11 связаны в воспроизводимый structural pipeline; сформулированы недостающие numerical gates; установлено, почему adjoint matrix residual недостаточен без production-linearity test; выявлен неполный closure loss ledger.

Открытые вопросы: полная линейность production KKV, точность support при снижении порогов, cutoff stability, абсолютная физическая нормировка, физический kernel, физическое время и причинный конус.

**Следующий научный приоритет:** не увеличивать сложность физической интерпретации, пока L1–L3 не установят, что вычисляемый production map действительно представляет линейный оператор на объявленном конечном секторе.
