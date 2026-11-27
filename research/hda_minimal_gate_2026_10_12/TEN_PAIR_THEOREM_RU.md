# BQG: десять пар вместо бесконечного сканирования lapse (12 октября 2026)

**Статус: алгебраические теоремы — PROVED; физический HDA — OPEN.** Дополнение к [предыдущему 32D lower bound](RESULT_RU.md), без изменения fail-closed статусов проекта. Выполненная локальная Python-реплика и отчёт находятся в исследовательском пакете текущего диалога.

## Точная теорема K5

Пусть H[N]=sum_i N_i H_i, i=0,...,4, где каждый H_i=H_i^E+lambda_L H_i^L задан на одном общем инвариантном домене. Пусть независимо выбранный D[sharp omega] линеен по omega, а sharp не подбирается после измерения HH.

Для ориентированного ребра ij введём omega_ij=N_i M_j-N_j M_i. Тогда

```
[H[N],H[M]]/(i*hbar)-D[sharp omega]
  = sum_(i<j) omega_ij { [H_i,H_j]/(i*hbar)-D[sharp e_ij] } .
```

Таким образом, **10 операторных равенств** на общем домене необходимы и достаточны для HDA на **всех пятиузловых** функциях lapse. Они не доказывают HDA при переходе к более крупным графам или континууму. Одно уравнение на паре (0,1) недостаточно.

Для любого фиксированного psi обозначим r_ij=([H_i,H_j]/(i*hbar)-D_ij)psi и G_(ij)(kl)=<r_ij,r_kl>. Тогда

```
||sum omega_ij r_ij||^2 = omega^T G omega
||sum omega_ij r_ij|| <= sqrt(lambda_max(G)) * ||omega||,
||omega||^2 = ||N||^2 ||M||^2 - (N dot M)^2 .
```

Достаточно десяти *полных* выходных векторов, их попарных скалярных произведений и диагонализации 10x10 Gram. Это верхняя оценка для одной psi; для операторной теоремы требуется общий домен и uniform estimate по физическим состояниям и уровню уточнения. Допустимые omega являются decomposable bivectors, удовлетворяют 5 Plücker identities; spectral bound на всей R^10 несколько консервативен.

## Разложение полного Lorentzian commutator

```
[H_i,H_j] = [E_i,E_j]
 + lambda_L ([E_i,L_j]+[L_i,E_j])
 + lambda_L^2 [L_i,L_j] .
```

Это точная полиномиальная структура; готовые амплитуды только евклидова H_E нельзя считать результатом для полного H. Нужны все 24 epsilon terms и invariant domain. Upstream [LORENTZIAN_HH_REACHABLE_SPACE.md](../../LORENTZIAN_HH_REACHABLE_SPACE.md) сообщает **границу поддержки**, не амплитуды: 11 314 085 конечных spin-network states для одной пары в консервативном перечислении и transient Jmax >= 13/2.

## Жёсткая алгебраическая преграда для D внутри vertex-only C^5

В конечной коммутативной алгебре A=C^5 с поточечным произведением всякая derivation D:A->A с обычным Leibniz тождеством нулевая: e_i^2=e_i => D(e_i)=2e_i D(e_i); все координаты D(e_i) равны нулю. Это НЕ запрещает матричные/noncommutative derivations, метрические displacement-операторы на гладком habitat или конечные дифференциальные calculus.

Наоборот, d:A->C^{edges} с (df)_ij=f_j-f_i ненулевая, а **точное** midpoint product rule:

```
d(f*g) = average(f)*d(g) + average(g)*d(f).
```

Независимый HDA covector \(\omega=\bar N dM-\bar M dN\) при этом точно равен \(N_iM_j-N_jM_i\). Он **не определяет** сам physical D: отсутствуют общий full graph-changing habitat, квантовая sharp_{E,q} и operator ordering.

## Воспроизведение и результаты

```bash
python research/hda_minimal_gate_2026_10_12/ten_pair_gate.py
```

SymPy: 10 basis edge pairs, пять Plücker identities, exterior Gram identity, нулевая размерность pointwise Leibniz derivation algebra и точный midpoint differential — PASS.

Синтетический 6x6 Hermitian control с 250 независимыми N,M: разложение по парам и 10x10 Gram — PASS с относительными погрешностями порядка 1e-15. D в тесте синтетический и выбран **независимо** от H; это проверка вычислительного алгоритма, не результат HDA для BQG.

**Физический следующий gate:** построить полный H_E+lambda_L H_L и отдельно D с одним исходным и выходным basis, spin/link/route embeddings, flux sharp, common domain, затем 10 различных r_ij, Gram и масштабную convergence. Отсутствие D означает BLOCKED, не нулевую невязку и не PASS.
