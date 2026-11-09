> **ИСТОРИЧЕСКОЕ НАУЧНО-ПОПУЛЯРНОЕ ПОВЕСТВОВАНИЕ.** Сохранено без удаления авторских глав. Не является действующей ведомостью физической доказанности. Некоторые заголовки и интерпретации, особенно spin-2, HDA, mirror-force и antigravity, требуют сверки с актуальными `theory_gates.json` и `physicalization_gates.json`. Канонический научный обзор — в корневом `README.md` и `THEORY_STATUS.md`.

---

# От бита к пространству-времени

## Увлекательный урок «на пальцах»: как из двух вариантов могут вырасти пространство, время, гравитация, квантовая пена и зеркальная хиральность

Представьте, что в самом глубоком слое природы **ещё нет привычного пространства**. Нет метров, координатной сетки, готового трёхмерного фона и даже заранее заданной геометрии.

Есть только элементарные различия: один переход может пойти одним способом или другим.

Самая маленькая единица такого различия — **бит**:

```text
0   или   1
```

Главная идея проекта очень проста по формулировке и очень трудна по содержанию:

> если огромное число бинарных различий соединять по одному и тому же локальному правилу, может ли их коллективное поведение само стать пространством, временем и гравитацией?

Мы не начинаем с уравнений Эйнштейна и не вставляем готовое четырёхмерное пространство в модель. Мы пытаемся пройти путь в обратную сторону:

```text
бит
 -> локальная комбинаторика
 -> сфера вокруг точки
 -> глобальное 3D-пространство
 -> релятивистское время
 -> гладкий 4D-like предел
 -> квантовая SU(2)-геометрия
 -> spin-2 коллективный сектор
 -> гамильтонова гравитация
 -> HDA
 -> вакуумные флуктуации и квантовая пена
 -> зеркальная ориентация / chirality
 -> 16-cell mirror order
 -> здоровый mirror-force candidate
```

Ниже — этот маршрут как последовательный урок: сначала школьная интуиция, затем строгая математика.

---

# Часть I. Как из битов появляется пространство

## Урок 1. Один бит — это ещё не пространство

Бит умеет различать только два состояния.

Сам по себе он не знает, что такое «лево», «право», «далеко», «рядом» или «три измерения».

Поэтому первый настоящий вопрос звучит так:

> сколько независимых бинарных различий нужно локальному переходу, чтобы отношения между ними начали вести себя как геометрия?

Обозначим это число через `q`.

Для `q` бинарных выборов существует

```text
2^q
```

различных маршрутов.

---

## Урок 2. Почему бинарное правило само выделяет q = 2

У локальной route-оболочки есть два специальных полюса и `2^q` бинарных route-состояний.

Каждый route отличается одним битом от `q` соседних routes и связан ещё с двумя полюсами. Поэтому степень route-вершины равна

```text
q + 2
```

А каждый полюс связан со всеми `2^q` маршрутами:

```text
2^q
```

Если потребовать локальную однородность, получаем простое уравнение:

```text
q + 2 = 2^q
```

Для целых `q >= 1` единственное решение:

```text
q = 2
```

Проверка буквально школьная:

```text
q=1:  3 != 2
q=2:  4  = 4
q=3:  5 != 8
q=4:  6 != 16
...
```

После `q=2` экспонента `2^q` растёт быстрее линейной части `q+2`, поэтому нового пересечения уже не появляется.

Это первый неожиданный момент проекта:

> **мы не назначаем два бинарных направления вручную — минимальная однородная route-оболочка сама выбирает q=2.**

---

## Урок 3. Из четырёх маршрутов возникает сфера вокруг точки

![Walsh-тетраэдр](figures/02-walsh-tetrahedron.svg)


При `q=2` routes имеют вид

```text
00, 01, 10, 11
```

Соединим маршруты, отличающиеся ровно одним битом. Получается квадрат `C4`.

Теперь добавим два полюса, каждый из которых связан со всеми четырьмя route-состояниями.

Получается граф октаэдра:

```text
Sigma Q2 = octahedral graph
```

У его поверхности:

```text
V = 6
E = 12
F = 8
chi = 2
```

То есть это обычная сфера `S2`.

Почему это так важно?

В хорошем трёхмерном многообразии маленькая сфера вокруг обычной внутренней точки должна быть именно двумерной сферой `S2`.

Значит локальный бинарный объект неожиданно имеет **правильную топологическую ссылку для точки трёхмерного пространства**.

Есть ещё одна красивая связь: октаэдральный граф равен line graph тетраэдра `K4`. Его шесть вершин можно читать как шесть рёбер тетраэдра. Поэтому route-комбинаторика и минимальная тетраэдральная квантовая геометрия встречаются в одной шестисостоянийной структуре.

---

## Урок 4. Как локальные сферы складываются в глобальное 3D-пространство

![От q=2 к 16-cell](figures/03-q2-to-16cell.svg)


Локальной `S2` мало. Нужен глобальный объект без края, дыр и сингулярных швов.

Каноническое PL-завершение проекта — граница четырёхмерного cross-polytope, то есть **16-cell**.

Для seed-комплекса:

```text
(V,E,F,T) = (8,24,32,16)
Betti     = (1,0,0,1)
```

Проверяются:

- vertex links;
- edge links;
- face links;
- двухсторонность треугольников;
- ориентируемость;
- `boundary^2 = 0`;
- гомологии;
- сохранение этих свойств после recursive PL refinement.

На проверенных уровнях всё остаётся согласованным с `S3`.

### Почему minimal+flag gluing почти неизбежно даёт именно 16-cell

Если дополнительно использовать фактическую simplicial semantics проекта — **минимальные 8 вершин + flag/clique closure** — появляется короткое доказательство уникальности этой глобализации с точностью до переименования вершин.

Каждая вершина должна иметь октаэдральную ссылку, значит среди остальных семи вершин она соединена ровно с шестью.

Следовательно у каждой вершины есть ровно **одна antipodal вершина**, с которой ребра нет.

Граф отсутствующих рёбер состоит из четырёх независимых пар:

```text
4 antipodal pairs
```

То есть 1-skeleton глобального комплекса — это

```text
K8 minus 4 antipodal edges
```

А flag closure разрешает тетраэдр тогда и только тогда, когда мы выбираем по одной вершине из каждой antipodal пары.

Число таких тетраэдров:

```text
2^4 = 16
```

И снова получаем

```text
(V,E,F,T) = (8,24,32,16)
```

— границу 16-cell.

Точный scope результата важен: **minimal 8-vertex flag globalization уникальна**, но произвольная bare-causal nonflag completion без этих правил всё ещё может быть неединственной.

---

## Урок 5. Но действительно ли пространство трёхмерное?

![Размерность и время](figures/04-emergent-dimensions.svg)


Топология говорит, как клетки склеены. Размерность требует отдельной проверки роста объёма и диффузии.

Если пространство трёхмерное, число состояний внутри радиуса `r` должно масштабироваться как

```text
V(r) ~ r^3
```

Frozen held-out результат:

```text
d_H = 2.999229782
```

Независимая spectral-проверка пространственного среза:

```text
d_s(slice) = 3.004393867
```

Один измеритель смотрит на рост объёма, другой — на распространение диффузии, а оба дают практически три измерения.

---

## Урок 6. Откуда появляется время

Теперь системе нужно не только «быть пространством», но и меняться.

Пусть пространственный масштаб увеличивается в `lambda_l` раз, а характерный временной масштаб — в `lambda_t` раз. Вводится динамический показатель `z`:

```text
lambda_t ~ lambda_l^z
```

Для релятивистского scaling нужен

```text
z -> 1
```

Frozen тест даёт

```text
z = 0.998281156
```

То есть пространство и время масштабируются почти одинаково.

Добавляя такое causal/history направление к трёхмерному spatial slice, получаем spectral dimension истории

```text
d_s(history) ~ 4.004393867
```

И первая большая цепочка замыкается:

```text
binary routes
 -> q=2
 -> local S2
 -> global PL S3
 -> d_slice ~ 3
 -> z ~ 1
 -> d_history ~ 4
```

---

# Часть II. Почему дискретный мир может выглядеть гладким

## Урок 7. Пиксели пространства исчезают при coarse-graining

Представьте цифровую фотографию. Вплотную видны отдельные пиксели. Издалека — гладкая поверхность.

Так же и здесь: макроскопический наблюдатель не различает отдельные микросостояния, а видит блоки размера `b`.

В вычислительной модели дефекты гладкой геометрии убывают почти идеальными степенями:

```text
delta g       ~ b^-2.001707
grad(delta g) ~ b^-3.001458
delta R       ~ b^-4.000524
simplicity    ~ b^-1.994838
Urbantke g    ~ b^-2.019746
```

То есть крупномасштабная гладкость не вводится командой «сделать пространство гладким». Она возникает потому, что множество микроскопических различий взаимно усредняются.

Но здесь появляется новый вопрос:

> если средняя геометрия становится гладкой, обязаны ли микроскопические флуктуации исчезнуть полностью?

Нет.

Можно иметь одновременно

```text
<delta g> = 0
<delta g^2> > 0
```

И это уже дверь к квантовой пене.

---

## Урок 8. Квантовая пена: не «хаос пространства», а флуктуации вокруг гладкого среднего

![Coarse-graining и квантовые флуктуации](figures/13-coarse-graining-foam.svg)


Правильная картинка выглядит так:

```text
smooth mean geometry
+ zero-mean quantum fluctuations
= microscopic quantum foam
```

Важное различие:

- **бит / qubit** задаёт дискретное внутреннее различие или поляризацию;
- **bosonic oscillator / Fock sector** задаёт амплитуду и число квантов;
- вакуум может иметь нулевое среднее поле, но ненулевую дисперсию.

Поэтому «пространство кипит» не означает, что макроскопическая метрика хаотична. Микроскопические fluctuations могут быть сильными, но коррелированными так, что на больших масштабах почти исчезают.

### Условное новое предсказание из уже измеренного smoothing exponent

Если — и это дополнительная физическая гипотеза — измеренный закон

```text
delta g_RMS ~ R^-2.001707
```

интерпретировать именно как RMS стационарной трёхмерной квантовой metric fluctuation, а её low-k spectrum записать как

```text
P_delta_g(k) ~ k^n
```

то в трёх пространственных измерениях

```text
RMS ~ R^-(3+n)/2
```

и поэтому

```text
n = 2*2.001707 - 3
  = 1.003414
```

Получается кандидатное предсказание:

```text
P_foam(k) ~ k^1.003414      as k -> 0
```

Это **не белый шум**. Для white noise было бы `n=0` и RMS exponent `3/2`.

Если эта интерпретация подтвердится прямым вычислением quantum two-point function, вакуум проекта окажется hyperuniform-like: микроскопические флуктуации есть, но длинноволновый шум сильно подавлен.

Именно так квантовая пена может быть бурной на микроуровне и почти идеально гладкой на больших масштабах.

Этот вывод в проекте имеет статус **conditional prediction**, а не доказанного свойства физического вакуума.

---

# Часть III. Как из qubits появляется квантовая геометрия и spin-2

## Урок 9. Где появляется SU(2)-геометрия

![SU(2) и Peter–Weyl](figures/05-su2-peter-weyl.svg)


Гравитации недостаточно знать только граф связей. Нужны площади, объёмы, ориентации и генераторы деформаций.

Здесь появляется `SU(2)` и Peter-Weyl язык.

На link/vertex-секторах строятся:

```text
flux operators
intertwiners
volume V
Euclidean Hamiltonian H_E
K = [V, H_E]
C(V) = h [h^-1, V]
C(K) = h [h^-1, K]
```

Finite gates отдельно проверяют:

- gauge covariance;
- Gauss sectors;
- настоящий volume operator;
- spin-support;
- charged recoupling;
- отсутствие запрещённой leakage;
- exact cutoff staircase.

Это уже не просто комбинаторика графа, а квантовая геометрия.

---

## Урок 10. Четыре микроскопических qubits содержат spin-2

Четыре spin-1/2 qubits разлагаются по полному угловому моменту так:

```text
(1/2)^4 = 2 x j=0 + 3 x j=1 + 1 x j=2
```

По размерностям:

```text
2 + 9 + 5 = 16
```

То есть в четырёх-qubit Hilbert space существует **ровно один j=2 irrep**.

Для него

```text
J^2 = j(j+1) = 6
```

а крайние состояния имеют

```text
m = +2
m = -2
```

Finite gate проверяет, что

```text
|up up up up>
```

и

```text
|down down down down>
```

лежат в `j=2` секторе с нулевой projector leakage.

Почему четыре — минимальное число?

Для `N` spin-1/2 объектов максимальный полный spin равен

```text
j_max = N/2
```

Чтобы получить `j=2`, нужно минимум

```text
N = 4
```

---

## Урок 11. Почему физическая поляризация гравитона — это один logical qubit

У массивного spin-2 объекта было бы пять spin projections. Но у безмассового гравитона gauge redundancy убирает нефизические компоненты и оставляет две распространяющиеся helicity:

```text
h = +2
h = -2
```

Поэтому физическая polarization space одного massless spin-2 кванта двумерна:

```text
|R> = |h=+2>
|L> = |h=-2>
```

То есть **поляризация гравитона — logical qubit**.

Это не означает, что весь гравитон — обычный двухуровневый объект. Полный mode всё равно бозонный:

```text
Fock occupation  x  C^2_helicity
```

Именно Fock sector даёт vacuum occupation statistics, а helicity-qubit говорит, в какой spin-2 polarization живёт возбуждение.

Точная цепочка:

```text
4 microscopic spin-1/2 qubits
 -> collective j=2 sector
 -> massless TT reduction
 -> helicity +2/-2
 -> one logical graviton-polarization qubit
```

---

## Урок 12. Что такое «инфотон» в этой кандидатной теории

Слово `infoton` здесь используется только как **внутреннее название проекта** для возможного бозонного collective mode информационного/route-сектора.

Сам факт бинарной polarization ещё не делает его spin-2.

Есть два случая:

```text
scalar/vector information mode
    -> не может линейно смешиваться с graviton h=+/-2

TT rank-2 information mode
    -> может иметь j=2 и helicity +2/-2
```

Отсюда появляется сильный falsifier:

> если route/information sector не содержит ненулевого TT spin-2 компонента, прямой graviton-infoton helicity mixing запрещён representation mismatch.

То есть идея резонанса не является «всё может взаимодействовать со всем». Она требует конкретного spin-2 канала.

---

# Часть IV. Как из квантовой геометрии получается каноническая гравитация

## Урок 13. Зачем нужен Lorentzian член

Евклидова часть `H_E` ещё не является полной лоренцевой гравитацией для вещественного параметра Барберо–Иммирци.

Полный локальный geometry operator имеет вид

```text
G_v = H_E,v + (1 + beta^2) H_L,v
```

где `H_L` строится из ковариантных `K`- и `V`-legs.

Отдельный classical control фиксирует правильную beta-комбинацию — её нельзя подгонять под HDA.

Support-анализ показывает, что для all-`j=1/2` input полный Lorentzian HH-проход безопасен при

```text
Jmax = 13/2
```

То есть достижимое пространство при таком cutoff конечно и контролируемо.

---

## Урок 14. Самый жёсткий тест: HDA

![HDA gate](figures/07-hda-gate.svg)


Общая теория относительности — это не только одна формула действия. Её гамильтоновы ограничения должны правильно описывать композицию деформаций пространственного среза.

Схематически:

```text
[ H[N], H[M] ]
      ->
 i*hbar * D[ sharp(N dM - M dN) ]
```

Интуиция простая.

Сначала слегка «подтолкните» spatial slice по нормали с профилем `N`, затем с профилем `M`.

Теперь поменяйте порядок.

Разница двух последовательностей не должна создавать новую степень свободы. Она должна быть просто **сдвигом вдоль самого пространства**.

Именно этот tangential shift генерирует `D`.

В проекте объект

```text
sharp(N dM - M dN)
```

строится независимо через discrete cochain/Hodge/flux map, а path-register реализует соответствующий diffeomorphism generator.

Отдельный factorization no-go показывает: если geometry Hamiltonian действует только на geometry sector и совсем не чувствует route register, правильный off-shell HDA получить нельзя.

Поэтому route-normal sector — не украшение, а необходимый канал.

---

# Часть V. Кульминация: почему лишние Lorentzian каналы исчезают

## Урок 15. Последняя проверка fixed-cutoff HDA в трёх строках

На frozen habitat берём

```text
N = Nbar + epsilon*n
M = Mbar + epsilon*m
Omega_Q = epsilon^-1 * OmegaTilde_Q
```

### Шаг A. Geometry x route

На первый взгляд cross-коммутатор содержит опасный член порядка `epsilon^-1`.

Но ведущий коэффициент одинаков в двух противоположных lapse-orderings и сокращается **state by state**.

Поэтому

```text
C_cross = O(1)
```

### Шаг B. Geometry x geometry

Для двух узлов antisymmetric lapse-smear имеет вид

```text
N0*M1 - N1*M0 = O(epsilon)
```

Следовательно

```text
C_GG = O(epsilon)
```

### Шаг C. Правильный diffeomorphism target растёт как inverse regulator

На frozen WKB carrier физическая производная даёт `epsilon^-1`, поэтому

```text
D[ sharp(N dM - M dN) ] = O(epsilon^-1)
```

Делим паразитные geometry-каналы на правильный target:

```text
C_cross / D = O(epsilon)
C_GG    / D = O(epsilon^2)
```

И получаем certificate:

```text
Delta_full(epsilon)
 <=
Delta_route(epsilon)
 + C_cross*epsilon
 + C_GG*epsilon^2

            -> 0
```

То есть при fixed regulator-safe Peter-Weyl cutoff полная локальная геометрия

```text
H_E + (1 + beta^2) H_L
```

не создаёт surviving off-shell anomaly поверх уже независимо проверенного route-normal HDA limit.

Поэтому прямой `11.3M-state [H_L,H_L]` brute-force остаётся полезным regression cross-check, но не является логически обязательным для fixed-cutoff composition theorem.

---

## Урок 16. Можно ли одновременно убирать regulator и поднимать spin cutoff?

Fixed cutoff — сильный результат, но естественный следующий вопрос:

> можно ли сделать `epsilon -> 0` и одновременно `Jmax -> infinity`?

Получается консервативная диагональная оценка.

На fixed-valence SU(2) узле:

```text
||E||    = O(J)
||V||    = O(J^(3/2))
||H_E||  = O(J^(3/2))
||K||    = O(J^3)
||C(K)|| = O(J^3)
||H_L||  = O(J^(15/2))
```

Для nondegenerate fixed-shape family flux metric масштабируется как

```text
Q = O(J^2)
```

поэтому route target имеет scale

```text
D = O(J^2 / epsilon)
```

Из уже доказанного cancellation structure следует грубая верхняя оценка

```text
C_cross / D = O(epsilon * J^(13/2))
C_GG    / D = O(epsilon^2 * J^13)
```

Положим

```text
Jmax(epsilon) = epsilon^-alpha
```

Оба лишних канала исчезают, если

```text
0 < alpha < 2/13
```

То есть существует допустимый simultaneous diagonal path, например

```text
Jmax = epsilon^-1/8
```

для которого

```text
C_cross / D = O(epsilon^(3/16))
C_GG    / D = O(epsilon^(3/8))
```

Это снимает fixed-cutoff ограничение **вдоль целого класса диагональных траекторий**:

```text
Jmax = o(epsilon^-2/13)
```

Но это ещё не uniform theorem для абсолютно любого способа одновременно отправлять `Jmax -> infinity` и `epsilon -> 0`.

---

# Часть VI. Почему гравитационные волны могут «раскачивать» квантовую пену

## Урок 17. Гравитационная волна как pump

![Гравитационная волна как pump](figures/14-gravitational-wave-pump.svg)


Пусть существует бозонный information/route mode с частотой `omega_I` и ненулевым TT coupling к метрике.

Проходящая gravitational wave может периодически менять его effective frequency:

```text
q_ddot
 + omega_I^2 * [1 + xi*h*cos(Omega_GW*t)] * q
 = 0
```

Это Mathieu equation.

Первая parametric instability band находится около

```text
Omega_GW ~= 2*omega_I
```

При слабой модуляции:

```text
mu ~= |xi*h| * omega_I / 4
```

где `mu` — Floquet/squeezing growth exponent, а `xi` — микроскопический coupling, который теория пока **не вывела**.

Finite gate использует тестовую модуляцию

```text
xi*h = 0.02
omega_I = 1
```

и получает на резонансе

```text
mu_numeric = 0.004999941407...
mu_leading = 0.005
relative error ~ 1.17e-5
```

А при `Omega=1.8*omega_I` и `2.2*omega_I` Floquet growth исчезает на gate tolerance.

То есть сам математический механизм резонанса воспроизводится очень точно.

Но это **не доказательство существования физического инфотона** и не измерение `xi`.

---

## Урок 18. Почему это не «бесплатная энергия из вакуума»

GW здесь является pump — источником энергии.

На резонансе эффективный квантовый Hamiltonian имеет squeezing-форму

```text
H_sq ~ i*hbar*kappa/2 * (b_dagger^2 - b^2)
```

Он эрмитов при правильной фазовой записи, поэтому его эволюция унитарна.

Вакуум не превращается в энергию «из ничего». Он превращается в squeezed state, а энергия приходит из внешней gravitational-wave pump.

Идея поэтому звучит так:

```text
existing vacuum fluctuations
 + coherent gravitational-wave pumping
 -> squeezed / amplified fluctuation quadrature
```

---

## Урок 19. Два разных резонанса — важный экспериментальный fingerprint

Если information mode действительно имеет тот же spin-2 TT quantum number, возможны два разных механизма.

### 1. Direct spin-2 conversion

Bilinear mixing вида

```text
a_g^dagger a_I + a_I^dagger a_g
```

резонирует около

```text
omega_g ~= omega_I
```

То есть это **1:1 resonance**.

### 2. Vacuum parametric squeezing

Периодическая modulation частоты даёт

```text
Omega_GW ~= 2*omega_I
```

То есть **2:1 resonance**.

Если когда-нибудь будет построен физический infoton sector и вычислен coupling, эти два разных положения резонанса дадут конкретный спектральный falsifier.

---

# Часть VII. Зеркальная ветвь: что действительно переворачивается

## Урок 20. В логическом geometry-qubit есть точный ориентационный бит

У четырёх spin-1/2 face qubits gauge-invariant singlet sector двумерен. В естественном basis два shape-оператора выражаются через `X_L` и `Z_L`, а ориентированный triple product имеет вид

```text
Q = (sqrt(3)/4) * Y_L
```

Теперь возьмём mirror operation как complex conjugation:

```text
M = K
 i -> -i
```

Тогда точно:

```text
X_L -> +X_L
Z_L -> +Z_L
Y_L -> -Y_L
Q   -> -Q
```

То есть mirror не стирает геометрию. Он меняет **только ориентационный знак**.

Естественный binary label:

```text
chi = sign(Q) = +1 or -1
```

Точные eigenvalues:

```text
Q_+ = +sqrt(3)/4 = +0.433012701892...
Q_- = -sqrt(3)/4 = -0.433012701892...
```

Новый finite gate проверяет, что complex conjugation меняет эти состояния местами с overlap

```text
0.9999999999999998
```

при этом shape observables и absolute-volume information совпадают до machine precision.

Так возникает строгая версия идеи «биты разворачиваются от нуля в другую сторону»:

```text
bit state
 + orientation chi=+1/-1
```

Но `chi=-1` пока означает **mirror orientation**, а не отрицательную энергию.

---

## Урок 21. Что значит «число pi идёт в зеркальном порядке»

Само число

```text
pi = 3.14159265...
```

не становится другим и его десятичные цифры не переворачиваются.

Зато у фазы есть направление обхода:

```text
theta -> -theta
```

или

```text
exp(i*theta) -> exp(-i*theta)
```

Поэтому можно ввести **ориентированный угол**

```text
pi_chi = chi*pi
```

где

```text
pi_+ = +pi
pi_- = -pi
```

Два пути вокруг unit circle противоположны по orientation, но их endpoint одинаков:

```text
exp(+i*pi) = exp(-i*pi) = -1
```

Новый gate проверяет `theta -> -theta` через complex conjugation точно; остаток на `+pi/-pi` равен лишь floating-point roundoff около

```text
2.45e-16
```

То есть физически осмысленная версия «зеркального pi» — **обратная ориентация фазы**, не обратный порядок цифр.

---

## Урок 22. Зеркальный тетраэдр имеет другую ориентацию, но ту же метрику

Пусть три edge-вектора образуют matrix

```text
A = (a,b,c)
```

а intrinsic metric

```text
G = A^T A
```

Возьмём reflection `R`:

```text
R^T R = I
det(R) = -1
A' = R A
```

Тогда

```text
G' = A'^T A'
   = A^T R^T R A
   = A^T A
   = G
```

но

```text
det(A') = -det(A)
```

То есть:

```text
orientation flips
metric stays the same
absolute volume stays the same
```

Новый finite gate прогоняет 256 deterministic random nondegenerate tetrahedra и получает:

```text
max metric-Gram error       = 0
max absolute-volume error   = 0
max orientation-flip error  = 0
max face-flux-Gram error    = 0
```

Это очень важный результат.

**Текущая микрогеометрия умеет различать зеркальную ориентацию, но физическая metric geometry не меняет знак.**

---

## Урок 23. Почему из зеркальности пока НЕ следует антигравитация

![Зеркальная ориентация и no-go](figures/15-mirror-chirality-no-go.svg)


В ADM metric variables

```text
g_00 = -N^2 + q_ab N^a N^b
```

Здесь нет отдельного множителя `chi`.

И Hamiltonian, который выбран HDA/DeWitt gate, строится из

```text
q_ab
pi_ab*pi^ab
pi^2
R[q]
```

то есть из parity-even metric contractions.

Поэтому для mirror-related configurations с одной и той же физической metric data текущая архитектура даёт

```text
g_00(+chi) = g_00(-chi)
```

и, следовательно,

```text
delta g_00^mirror = 0
```

Это **negative result**, но очень полезный:

> ориентационный бит существует, однако сам по себе он не переворачивает гравитационный знак.

В существующей ветке

```text
mirror orientation != negative gravitational mass
```

---

## Урок 24. Какой именно новый эффект потребовался бы для антигравитации

Пусть будущая теория добавит orientation-odd gravitational response:

```text
a_chi = a_even + chi*a_odd
```

Калибруем обычную `chi=+1` ветвь на измеренное Newton gravity:

```text
a_+ = g_N
```

Введём долю

```text
f = a_odd / g_N
```

Тогда противоположная mirror branch имеет

```text
a_- / g_N = 1 - 2*f
```

И сразу получаем точные пороги:

```text
f = 0      -> обе ориентации гравитируют одинаково
f = 1/2    -> полное экранирование mirror branch
f > 1/2    -> mirror branch становится repulsive
f = 1      -> равная по модулю противоположная gravity
```

Текущий mirror-even geometry gate соответствует

```text
f_current = 0
```

Поэтому antigravity не «спрятана» в уже построенной метрике. Для неё нужен **новый оператор**.

Самый простой mirror-odd local geometry operator — `Y_L`:

```text
H_odd = lambda_chi * Y_L
```

Но здесь появляется главный killer-test.

Полный constraint станет

```text
H = H_0 + lambda_chi*H_chi
```

а его algebra:

```text
[H,H]
 = [H_0,H_0]
 + lambda_chi * ([H_0,H_chi] + [H_chi,H_0])
 + lambda_chi^2 * [H_chi,H_chi]
```

Каждый новый канал обязан снова замкнуться в правильный diffeomorphism generator либо образовать новую согласованную first-class constraint.

Поэтому настоящий antigravity criterion такой:

```text
nonzero mirror gravity response
AND
stable Hamiltonian
AND
first-class modified HDA
```

Без последнего условия можно получить красивый знак, но не согласованную гравитационную теорию.

---

# Часть VIII. Зеркальность, материя и аномалии

## Урок 25. Chirality может возникнуть без отрицательной энергии

Mirror orientation особенно естественно связывается не с negative energy, а с **left/right chirality**.

У fermion axial density

```text
J5^0 = psi_bar gamma^0 gamma^5 psi
```

parity меняет знак.

Наш oriented geometry coordinate `Y_L` тоже mirror-odd.

Поэтому product

```text
H_chi-psi = lambda * Y_L * J5^0
```

parity-even.

Это даёт очень естественную кандидатную цепочку:

```text
orientation bit chi
 -> local handedness
 -> left/right fermion preference
```

при этом metric и energy sign могут оставаться обычными.

То есть зеркальная ветвь может оказаться гораздо полезнее для **происхождения хиральности материи**, чем для антигравитации.

---

## Урок 26. Mirror-conjugate pair автоматически меняет знак perturbative gauge anomaly

Для left-handed Weyl fermion в representation `R` cubic gauge anomaly пропорциональна

```text
d_R^(abc) = Tr[ T^a {T^b,T^c} ]
```

Для conjugate representation

```text
T_Rbar^a = -(T_R^a)^T
```

поэтому алгебраически

```text
d_Rbar^(abc) = -d_R^(abc)
```

и точная mirror-conjugate pair даёт

```text
d_R + d_Rbar = 0
```

Новый finite stress-test на deterministic noncommuting Hermitian generators получает residual

```text
~2.3e-15
```

То есть сам **знак anomaly действительно зеркалится**.

Это не вывод Standard Model и не доказательство всех global anomalies. Но это реальный механизм:

```text
chiral sector R
 + mirror-conjugate sector Rbar
 -> perturbative cubic anomaly cancellation
```

Главная открытая задача теперь не «может ли зеркало менять anomaly sign?» — может.

Главный вопрос:

> почему видимый low-energy мир остаётся chiral, а mirror sector не делает всю теорию просто vector-like и где этот sector физически находится?

---

## Урок 27. Антиматерия — не то же самое, что mirror orientation

Здесь особенно важно не смешать разные симметрии:

```text
C   charge conjugation: particle <-> antiparticle
P   parity: spatial mirror / left <-> right
chi project orientation: sign of oriented geometry Q
```

Физическая antiparticle имеет положительную excitation energy. Антиматерия в современной QFT не означает макроскопическую отрицательную энергию.

Есть и прямой experimental reality check: ALPHA Collaboration измерила motion antihydrogen в поле Земли. Антиводород ускорялся **к Земле**, согласованно с обычным притяжением; repulsive gravity magnitude `1g` для antihydrogen была исключена.

Reference:

```text
E. K. Anderson et al. (ALPHA Collaboration)
Observation of the effect of gravity on the motion of antimatter
Nature 621, 716-722 (2023)
DOI: 10.1038/s41586-023-06527-1
```

Поэтому если в CIMFIG когда-нибудь появится mirror-antigravity sector, его нельзя автоматически назвать обычной антиматерией. Он должен быть **другим квантовым сектором**, отличающимся от standard charge-conjugated matter.

---

# Часть IX. Как всё-таки построить здоровую mirror-repulsion

## Урок 28. Почему просто заменить H на -H не работает

Самый простой соблазн — объявить, что для зеркальной ветви

```text
H_chi[N] = s_chi * H_GR[N].
```

Но HDA тут же даёт строгий тест.

Если

```text
{H_GR[N],H_GR[M]} = D[beta],
```

то

```text
{H_chi[N],H_chi[M]} = s_chi^2 * D[beta].
```

Чтобы сохранить тот же `D` и ту же нормировку lapse, нужно

```text
s_chi^2 = 1.
```

Остаются только

```text
s=+1
s=-1.
```

Но

```text
-H[N] = H[-N].
```

То есть `s=-1` просто разворачивает нормальное/time направление эволюции. Статическая ньютоновская сила от этого не становится отталкивающей.

Итак:

```text
H -> -H
```

**не является антигравитацией**.

---

## Урок 29. Почему отрицательный Newton sign почти сразу создаёт ghost

Второй соблазн — сделать коэффициент при Einstein-Hilbert term отрицательным:

```text
S_grav ~ F_chi * integral sqrt(-g) R
```

и выбрать

```text
F_- < 0.
```

Тогда формально можно говорить об opposite sign effective Newton coupling.

Но этот же `F_chi` стоит перед кинетическим членом линейного spin-2 поля.

Поэтому

```text
F_- < 0
```

означает отрицательную kinetic energy гравитона относительно обычной positive-energy matter.

Получается ghost instability.

То есть прямой metric-sign flip даёт нужный знак ценой разрушения устойчивости.

Это второй no-go.

---

## Урок 30. Здоровая альтернатива: orientation order parameter sigma и mediator phi

Вместо отрицательной массы оставим Einstein gravity обычной и положительно-энергетической.

Микроскопически у нас уже есть

```text
chi = sign(Q)
Q ~ Y_L.
```

На coarse scale используем pseudoscalar order parameter

```text
sigma ~ +v  -> chi=+1
sigma ~ -v  -> chi=-1.
```

И второй pseudoscalar `phi`, который будет mediator.

Mirror действует как

```text
(phi,sigma) -> (-phi,-sigma).
```

Возьмём action

```text
S = integral sqrt(-g) [
      Mpl^2 R/2
    - (d phi)^2/2
    - (d sigma)^2/2
    - U(phi,sigma)
]
```

с

```text
U = mu^2 phi^2/2
  + lambda phi^4/4
  + g phi sigma
  + kappa (sigma^2-v^2)^2/4.
```

При

```text
mu^2 > 0
lambda > 0
kappa > 0
```

оба kinetic terms положительные, а quartic potential bounded at large field.

То есть две mirror branches существуют без отрицательной excitation energy.

### И теперь главный тест: HDA

Для canonical two-field sector

```text
H[N] = integral N [
    p_phi^2/2 + (phi')^2/2
  + p_sigma^2/2 + (sigma')^2/2
  + U(phi,sigma)
].
```

В antisymmetric bracket все local potential terms сокращаются.

Остаётся

```text
{H[N],H[M]}
 = D[N M' - M N'].
```

Новый spectral gate при `L=512` получает

```text
H-H bracket = 0.0007282211771021175
D target    = 0.0007282211771021123
abs error   = 5.204170427930421e-18
rel error   = 7.146414566848946e-15
```

То есть новый mirror/orientation matter sector воспроизводит canonical HDA principal identity практически до machine precision.

Это уже принципиальное отличие от произвольного `lambda_chi Y_L` term: мы построили continuum extension, у которой constraint structure заранее правильного типа.

---

## Урок 31. Как из orientation bit возникает реальное отталкивание

Пусть coarse object несёт orientation charge

```text
Q_chi = eta * m * chi.
```

Важно:

```text
m > 0
```

для обеих ветвей.

Меняется знак **charge**, а не знак энергии.

Exchange mediator `phi` даёт potential

```text
U_12(r)
 = -G_T m1 m2/r
   -alpha G_T m1 m2 chi1 chi2 exp(-m_phi r)/r.
```

Для одинаковых orientations

```text
chi1*chi2 = +1
```

новый канал добавляет притяжение.

Для противоположных

```text
chi1*chi2 = -1
```

он становится отталкивающим.

Положим

```text
x = m_phi*r.
```

Тогда magnitude orientation force относительно bare tensor gravity:

```text
F_chi/F_T = alpha*(1+x)*exp(-x).
```

Полное screening происходит ровно при

```text
alpha_crit(x) = exp(x)/(1+x).
```

А repulsion начинается при

```text
alpha > alpha_crit(x).
```

Несколько точных порогов:

```text
m_phi*r = 0.0 -> alpha_crit = 1.0000000000
m_phi*r = 0.1 -> alpha_crit = 1.0047008346
m_phi*r = 0.5 -> alpha_crit = 1.0991475138
m_phi*r = 1.0 -> alpha_crit = 1.3591409142
m_phi*r = 2.0 -> alpha_crit = 2.4630186996
m_phi*r = 5.0 -> alpha_crit = 24.7355265171
```

В long-range limit:

```text
m_phi*r << 1
```

получаем очень простой результат:

```text
alpha = 1 -> complete screening
alpha > 1 -> opposite-chi repulsion.
```

Тестовый пример:

```text
alpha = 2
m_phi*r = 0.1
```

даёт

```text
orientation repulsion / bare tensor gravity
= 1.990642319679...
```

поэтому после вычитания обычного притяжения остаётся

```text
net outward force / bare tensor gravity
= 0.990642319679...
```

Это уже настоящий математический **antigravity-like cross-sector effect**, но без negative mass и без ghost graviton.

---

## Урок 32. Почему это пока называется mirror-force, а не доказанная антигравитация

Новая конструкция не меняет непосредственно

```text
g_00 -> -g_00.
```

Она даёт

```text
tensor gravity
+
orientation-dependent fifth force.
```

Поэтому её правильное название на текущем этапе:

```text
healthy mirror-force candidate
```

или

```text
antigravity-like cross-sector repulsion.
```

Это даже сильнее с точки зрения устойчивости: spin-2 sector остаётся positive-energy, а repulsion несёт отдельный canonical mediator.

Но теперь возникает последний неудобный вопрос:

> откуда в микроскопической теории вообще берётся поле `sigma`? Не вставили ли мы его руками?

Оказалось, что minimal 16-cell даёт на этот вопрос очень конкретный ответ.

---

## Урок 33. Как сам 16-cell создаёт зеркальный order parameter sigma

![Mirror order Q4](figures/11-mirror-order-q4.svg)


У 16-cell есть четыре antipodal пары вершин. Каждый из его 16 тетраэдров выбирает по одной вершине из каждой пары.

Поэтому каждый тетраэдр естественно кодируется **четырьмя битами**:

```text
b = b1 b2 b3 b4
```

и всего таких клеток ровно

```text
2^4 = 16.
```

Два тетраэдра имеют общую треугольную грань тогда и только тогда, когда их строки отличаются ровно одним битом.

Значит dual graph этих 16 клеток — не произвольный граф, а точно

```text
Q4 = 4D hypercube
```

с

```text
16 vertices
32 edges
degree 4.
```

Это снова неожиданная встреча бинарности и геометрии.

### Почему локальные ориентации чередуются, но глобальный mirror bit всё равно существует

Правильная склейка соседних тетраэдров требует противоположных outward orientations. В logical geometry-qubit это соответствует

```text
Y_v * Y_w = -1
```

на каждом dual edge.

`Q4` — bipartite graph. Определим

```text
eta_v = (-1)^popcount(v).
```

На каждом ребре

```text
eta_v * eta_w = -1.
```

Теперь введём staggered variable

```text
sigma_v = eta_v * Y_v.
```

Тогда геометрическое условие

```text
Y_v * Y_w = -1
```

превращается в обычный uniform order:

```text
sigma_v * sigma_w = +1.
```

И появляется настоящий block order parameter

```text
Sigma = (1/16) * sum_v eta_v * Y_v.
```

У него есть две точные mirror-вакуумные ветви:

```text
Y_v = +eta_v  -> Sigma = +1
Y_v = -eta_v  -> Sigma = -1.
```

Mirror меняет их местами:

```text
Sigma -> -Sigma.
```

То есть цепочка

```text
Y_L / Q
 -> Sigma=+/-1
 -> coarse sigma(x)
```

теперь имеет **конкретный finite microscopic bridge**.

### Сколько энергии стоит испортить mirror order

Для orientation part gluing Hamiltonian

```text
H_Y = J * sum_<vw> Y_v Y_w
```

идеальная ветвь имеет

```text
E0 = -32J.
```

Если перевернуть один local orientation bit, портятся четыре соседних bond. Каждый стоит `2J`, поэтому

```text
Delta E_single = 8J.
```

Если сделать domain wall, перевернув половину hypercube, пересекаются восемь связей:

```text
Delta E_wall = 16J.
```

Новый gate получает эти числа **точно**.

### А переживает ли mirror order квантовые флуктуации?

Проверяется полный 16-qubit transverse-field Hamiltonian

```text
H = -J * sum_<vw> sigma_v sigma_w
    -h * sum_v X_v.
```

Его Hilbert space уже имеет

```text
2^16 = 65536
```

состояний.

Sparse Lanczos diagonalization при

```text
h/J = 0.2
```

даёт

```text
<Sigma^2>  = 0.9976539474
<|Sigma|>  = 0.9987478253
```

две нижние mirror-состояния образуют практически вырожденный doublet со splitting меньше `1e-12 J`, а gap до следующего уровня равен примерно

```text
7.970087877 J.
```

То есть на минимальном 16-cell block зеркальный order не только существует классически — он **очень устойчив квантово** при слабом flip field.

Для negative control увеличиваем quantum flipping до

```text
h/J = 4
```

и получаем

```text
<Sigma^2> = 0.1462741385.
```

Order разрушается, как и должен.

Поэтому `sigma` больше не является полностью придуманным continuum полем. Естественный кандидат:

```text
sigma(x)
 ~ coarse block average of eta_v Y_v.
```

Более того, soft fluctuations самого `Sigma/sigma` могут оказаться mediator mode, поэтому отдельный `phi` в консервативной двухполевой конструкции, возможно, не фундаментален.

Теперь главный bottleneck уже не

```text
откуда взять sigma?
```

а

```text
какова физическая normalization sigma
и как она связывается с matter charge?
```

Именно из этого должно выйти главное число

```text
alpha
```

без ручной подгонки.

---

# Феерический итог: вся теория в одной лестнице

Теперь весь рассказ можно собрать в одну цепочку:

```text
BITS
 -> q=2 binary selector
 -> C4 route shell
 -> octahedral S2 link
 -> minimal flag 16-cell globalization
 -> recursive PL S3
 -> d_H ~ 3
 -> d_s(slice) ~ 3
 -> z ~ 1
 -> d_s(history) ~ 4
 -> smooth IR metric
 -> SU(2) Peter-Weyl quantum geometry
 -> 4 qubits contain a unique j=2 sector
 -> helicity +2/-2 logical graviton qubit
 -> H_E
 -> K=[V,H_E]
 -> C(V), C(K)
 -> H_E + (1+beta^2) H_L
 -> route-normal generator
 -> sharp(N dM - M dN)
 -> HDA composition certificate
 -> admissible joint Jmax(epsilon) continuum paths
 -> quantum metric fluctuations
 -> candidate hyperuniform foam spectrum
 -> possible GW-driven squeezing/resonance
 -> mirror K: i -> -i
 -> orientation bit chi=+/-1
 -> Q -> -Q but metric -> same metric
 -> mirror-conjugate anomaly sign cancellation
 -> current metric antigravity falsifier: f_current = 0
 -> H->-H no-go: only time orientation reversal
 -> negative Einstein-Hilbert sign no-go: graviton ghost
 -> 16-cell tetrahedra form dual Q4
 -> staggered Sigma=(1/16)sum eta_v Y_v
 -> two exact mirror sectors Sigma=+/-1
 -> finite 16-qubit ordered mirror doublet
 -> coarse sigma(x)
 -> positive-kinetic mediator sector
 -> canonical mirror-sector HDA identity
 -> exact alpha_crit = exp(m_phi r)/(1+m_phi r)
 -> healthy cross-sector repulsion if alpha > alpha_crit
 -> physical normalization and first-principles alpha are the next killer gate
```

Самая короткая HDA-кульминация остаётся:

```text
C_cross = O(1)
C_GG    = O(epsilon)
D       = O(epsilon^-1)
```

поэтому

```text
C_cross / D = O(epsilon)
C_GG    / D = O(epsilon^2)
```

и

```text
Delta_full
 <= Delta_route
  + C_cross*epsilon
  + C_GG*epsilon^2
 -> 0
```

А если одновременно увеличивать spin cutoff достаточно медленно,

```text
Jmax = o(epsilon^-2/13)
```

то эта иерархия сохраняется вдоль допустимого joint-limit family.

Mirror calculation добавляет:

```text
M X M^-1 = +X
M Z M^-1 = +Z
M Y M^-1 = -Y
```

но

```text
metric(+chi) = metric(-chi).
```

Microscopic mirror-order gate добавляет:

```text
Sigma = (1/16) sum eta_v Y_v
Sigma = +/-1
Delta E_single = 8J
Delta E_wall   = 16J.
```

А healthy-force branch добавляет:

```text
{H_mirror[N],H_mirror[M]}
 = D_mirror[N dM - M dN]
```

и

```text
opposite-chi repulsion
<=>
alpha*(1+m_phi*r)*exp(-m_phi*r) > 1.
```

То есть **ориентация бинарна, стандартная метрика mirror-even, глобальный mirror order уже возникает на минимальном 16-cell block, а дополнительная положительно-энергетическая mirror-charge dynamics математически может давать отталкивание.**

---

# Что именно утверждает проект

## Статус: кандидатная теория / candidate theory

Репозиторий содержит несколько разных классов результатов, и их нельзя смешивать.

### Exact / algebraic

Например:

- бинарный selector `q=2`;
- SU(2) representation decompositions;
- `4 spin-1/2 qubits -> unique j=2 irrep`;
- operator-ordering identities;
- fixed-cutoff composition scaling theorem в заявленном habitat;
- mirror algebra `X,Z even; Y,Q odd`;
- conjugate-representation identity `d(Rbar) = -d(R)`;
- 16-cell tetrahedron dual graph `Q4`;
- exact staggered mirror vacua `Sigma=+/-1`;
- exact defect costs `8J` и `16J`;
- rescaling no-go `H_chi=sH_GR => same HDA only if s^2=1`;
- Yukawa screening threshold `alpha_crit=exp(x)/(1+x)`.

### Tested finite

Например:

- recursive PL manifold gates;
- held-out dimensional scaling;
- Peter-Weyl covariance/support;
- two-node geometry x route HDA regression;
- graviton-helicity qubit projector gate;
- Mathieu/Floquet resonance gate;
- 256-tetrahedron mirror metric/orientation gate;
- mirror-state swap and oriented-phase gate;
- numerical conjugate-anomaly cancellation stress-test;
- full 16-qubit `2^16=65536` mirror-order sparse diagonalization;
- ordered mirror doublet with `<Sigma^2>=0.9976539474` at `h/J=0.2`;
- disordered control `<Sigma^2>=0.1462741385` at `h/J=4`;
- two-pseudoscalar canonical HDA gate with relative defect `~7.15e-15`;
- mirror-force screening/repulsion threshold regression.

### Conditional physical interpretations

Например:

- отождествление старого `delta g` smoothing exponent с quantum RMS;
- из него следует кандидатный `P_foam(k) ~ k^1.003414`;
- существование физического TT infoton mode;
- ненулевой microscopic coupling `xi` к gravitational waves;
- возможный `Y_L J5^0` bridge к fermion chirality;
- continuum normalization и propagation law derived `Sigma/sigma` mode;
- coupling derived `sigma` к physical matter charge;
- существование отдельного physical mediator `phi` либо доказательство, что soft `sigma` mode сам является mediator;
- first-principles value `alpha`;
- microscopic Peter-Weyl x route x orientation HDA closure.

### Что расчёт уже опровергает внутри текущей архитектуры

Текущая mirror-even geometry **не** даёт

```text
chi -> -chi
   =>
g_00 -> -g_00.
```

Также простое

```text
H -> -H
```

не даёт статическую антигравитацию, а negative Einstein-Hilbert coefficient даёт ghost tensor sector.

То есть здоровая ветка должна быть сложнее прямого sign flip.

### Что теперь впервые построено положительно

Построен finite microscopic mirror-order bridge

```text
Y_L/Q -> staggered Sigma -> coarse sigma(x)
```

и continuum canonical mirror sector с positive kinetic energy и exact HDA principal identity. Он допускает orientation-charge Yukawa force, которая для opposite `chi` становится repulsive и может превысить tensor attraction при

```text
alpha > exp(m_phi*r)/(1+m_phi*r).
```

Это **кандидатный стабильный mirror-force mechanism**, а не экспериментально обнаруженная антигравитация.

**Проект не утверждает, что фундаментальная физика природы экспериментально доказана.**

Корректная формулировка результата:

> **В репозитории построена математически и вычислительно проверяемая кандидатная схема перехода от бинарной дискретной микроструктуры к 3D/4D-like гладкому пределу, SU(2)/Peter-Weyl квантовой геометрии и канонической HDA-гравитационной архитектуре. Дополнительно построены exact finite мост от четырёх microscopic qubits к spin-2 сектору, условная модель квантовой пены и GW-squeezing, mirror/chirality gate, finite 16-cell mirror-order bridge и continuum healthy mirror-force extension. Последняя сохраняет положительные kinetic terms, проходит canonical matter-HDA gate и допускает cross-sector repulsion выше точного порога `alpha_crit`; однако physical normalization и microscopic force coupling ещё не выведены. Это кандидатная теория, а не установленный экспериментальный закон природы.**

---

# Что ещё действительно остаётся открытым

![Девять физических ворот](figures/12-nine-physical-gates.svg)


После последних результатов список стал уже, но не исчез.

Остаются, среди прочего:

- настоящий uniform joint theorem для максимально общего `Jmax -> infinity`, `epsilon -> 0`, а не только admissible diagonal family;
- вывод Lorentzian quantum history measure и глобальной unitarity;
- вывод microscopic TT information-mode action и coupling `xi`;
- прямой quantum two-point function вакуума для проверки или опровержения `P(k) ~ k^1.003414`;
- реалистичная matter sector: gauge group, masses, generations, chirality и Yukawa structure;
- полная anomaly cancellation, включая global anomalies, не только exact `R/Rbar` cubic pair identity;
- объяснение, почему mirror sector скрыт/decoupled, если он существует;
- physical normalization, kinetic coefficient и propagation speed derived `Sigma/sigma` mode;
- derivation coupling of `Sigma/sigma` к matter и charge-to-mass coefficient `eta`;
- first-principles prediction of `alpha` и mediator mass/range;
- full Peter-Weyl x route x orientation quantum HDA regression;
- проверка fifth-force/equivalence-principle/cosmological bounds после появления physical scale;
- literal metric antigravity branch, если она всё ещё нужна после healthy mirror-force construction;
- physical scale и Newton constant;
- blind predictions в физических единицах;
- независимая внешняя репликация другой группой/кодом/экспериментом.

Отдельно: произвольная bare-causal nonflag globalization всё ещё может быть неединственной; доказанная выше уникальность относится к minimal 8-vertex flag semantics.

---

# Где смотреть строгие доказательства и проверки

- `THEORY_STATUS.md` — канонический человеческий ledger;
- `theory_gates.json` — machine-readable proof ledger;
- `FINAL_CORE_ARCHITECTURE_CERTIFICATE.md` — fixed-cutoff HDA composition theorem;
- `GRAVITON_INFOTON_FOAM_BRIDGE.md` — spin-2/helicity, foam и GW-resonance bridge;
- `MIRROR_CHIRALITY_GRAVITY.md` — mirror orientation, chirality, anomaly sign и antigravity falsifier;
- `MICROSCOPIC_MIRROR_ORDER.md` — Q4 dual 16-cell, staggered Sigma, defect/domain-wall energies и 16-qubit ordered phase;
- `ORIENTATION_ODD_HDA_CONSTRUCTION.md` — healthy mirror-force construction, HDA control и exact repulsion threshold;
- `SPATIAL_QUBIT_GEOMETRY_BRIDGE.md` — exact four-qubit geometry decomposition;
- `GLOBAL_MANIFOLD_Q2_COMPLETION.md` — q=2 PL globalization;
- `bcqg_bit_to_gravity_final.py` — executable canonical status aggregator;
- `bcqg_observer_smoothing_unified.py` — bit -> scaling -> smoothing;
- `bcqg_global_manifold_gate.py` — recursive PL-manifold gate;
- `scripts/path_normal_hda_gate.py` — route-normal HDA;
- `scripts/peter_weyl_two_node_euclidean_joint_gate.py` — two-node Peter-Weyl x route regression;
- `scripts/lorentzian_hit_depth_bound.py` — Lorentzian support wall;
- `scripts/graviton_infoton_foam_gate.py` — exact j=2 projector + Floquet resonance checks;
- `scripts/mirror_chirality_gravity_gate.py` — mirror/chirality/phase/anomaly finite gate;
- `scripts/mirror_order_16cell_gate.py` — microscopic 16-cell staggered mirror-order gate;
- `scripts/orientation_odd_hda_gate.py` — continuum mirror-sector HDA + Yukawa screening/repulsion gate.

Запуск канонического статуса:

```bash
python bcqg_bit_to_gravity_final.py --strict
```

Spin-2 / foam bridge:

```bash
python scripts/graviton_infoton_foam_gate.py
```

Mirror/chirality gate:

```bash
python scripts/mirror_chirality_gravity_gate.py \
  --trials 256 \
  --output verification_results/MIRROR_CHIRALITY_GRAVITY.json
```

Microscopic mirror-order gate:

```bash
python scripts/mirror_order_16cell_gate.py \
  --output verification_results/MIRROR_ORDER_16CELL.json
```

Healthy orientation-force gate:

```bash
python scripts/orientation_odd_hda_gate.py \
  --L 512 \
  --output verification_results/ORIENTATION_ODD_HDA.json
```

Если canonical aggregator возвращает

```text
core_candidate_architecture_closed: true
```

это означает закрытие заявленной **candidate architecture в её доказанном scope**, а не экспериментальное подтверждение окончательной теории природы.
