# Архивный аудит BQG — 10 октября 2026

> **Исторический срез, не текущий статус.** Этот аудит фиксирует состояние на 10 октября 2026 года. Для актуальной машинной истины используйте `theory_gates.json`, `physicalization_gates.json` и `OPEN_PROBLEMS.md`; для более поздних специализированных результатов — `BQG_QFI_Q2_EMERGENT_GRAPH_AUDIT_20261118.md` и `W5_PRODUCTION_AUDIT_20261113.md`.

Базовый срез: a55c77b. Полного доказательства BQG или вывода физической гравитационной динамики нет. Структурная замкнутость и физическая завершённость false. Корректный реестр может содержать научный no_go; этот статус не засчитывается как закрытый core gate.

## Научные результаты

- Точная кинематическая изометрия не создаёт физического пространства автоматически.
- На сохранённых actual32 Euclidean H0 образах minimum master eigenvalue3.97703148; конечные окна у нуля исчезают. Полный kernel вне code не исключён.
- Для source column0 выполнены настоящие H0H1 и H1H0: по1004 выходных состояний; commutator norm2.87945381; 96.35% squared norm вне исходного code.74 reverse matrix elements согласованы до2.78e−16.
- Эти расчёты относятся к sine Euclidean ordering с zero-aware volume, не к полному Lorentzian constraint. Общий selfadjoint domain и все32 two-hit columns не проверены.
- Независимое microscopic D на том же habitat с одинаковым lapse interpolation отсутствует в проверенной реализации. Ненулевой commutator не устанавливает anomaly.
- Условная цепочка spectral/refinement→connected history→Einstein IR сформулирована; её microscopic предпосылки не выведены.
- Девять физических gates остаются открытыми. Depth6 требует свидетелей [32],[311],[221].

## Исправления управления статусом

Validator допускает научный no_go, но не включает его в closure. Объявленная structural closure=false. Итоговый агрегатор учитывает все девять обязательных физических условий. Ancilla использует общий occupation-order spin basis. CI сохраняет свидетельство ожидаемого no-go и различает exit2 от ошибок. Добавлены используемые networkx/torch; зависимостям всё ещё нужны закреплённые версии для полной численной воспроизводимости.

## Свидетельства и границы воспроизводимости

[Пять маршрутов](research/audit_2026_10_10/FIVE_ROUTES_AUDIT.md), [условное доказательство](research/audit_2026_10_10/CONDITIONAL_PROOF.md), [two-hit расчёт](research/audit_2026_10_10/TWO_HIT_FINDINGS.md).

Команды из корня:

```sh
python scripts/verify_theory_gates.py
python scripts/verify_physicalization_gates.py
python scripts/verify_depth6_frontier.py
python scripts/bqg_ancilla_peter_weyl_refinement_gate.py
python research/audit_2026_10_10/replay_two_hit.py
python research/audit_2026_10_10/check_conditional_chain.py
```

Two-hit replay проверяет суммирование сохранённых ветвей, не заново выводит SU(2) amplitudes. Pickle загружать только из доверенного checkout. Архивные гигабайтные вычисления целиком не переиграны. Исторические отчёты сохранены как исторические свидетельства; эта страница задаёт текущую интерпретацию.
