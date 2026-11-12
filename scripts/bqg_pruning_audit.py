"""Opt-in, semantics-preserving audit of numerical sparse-state pruning.

The collector is context-local and disabled by default. Instrumented operators
continue to execute their original threshold decisions; when a coefficient is
discarded, the audit records the exact complex values that led to the decision.
"""
from __future__ import annotations

from collections import defaultdict
from contextvars import ContextVar
import math
from typing import Any

_ACTIVE: ContextVar["PruningCollector | None"] = ContextVar(
    "bqg_pruning_collector", default=None
)


class PruningCollector:
    def __init__(self) -> None:
        self.events: list[dict[str, Any]] = []
        self._by_site: dict[str, dict[str, float | int]] = defaultdict(
            lambda: {"count": 0, "pathwise_l1": 0.0, "pathwise_l2_sq": 0.0,
                     "max_abs": 0.0}
        )

    def record(self, site_id: str, key: Any, value: complex, threshold: float,
               *, kind: str, prior: complex | None = None,
               incoming: complex | None = None) -> None:
        z = complex(value)
        # Keep exact-zero cancellations as events too: they carry no lost
        # norm, but prove that a previously accumulated key was deleted.
        row = {
            "site_id": str(site_id),
            "kind": str(kind),
            "key_repr": repr(key),
            "real": float(z.real),
            "imag": float(z.imag),
            "abs": float(abs(z)),
            "threshold": float(threshold),
        }
        if prior is not None:
            row["prior_real"] = float(complex(prior).real)
            row["prior_imag"] = float(complex(prior).imag)
        if incoming is not None:
            row["incoming_real"] = float(complex(incoming).real)
            row["incoming_imag"] = float(complex(incoming).imag)
        self.events.append(row)
        stats = self._by_site[str(site_id)]
        stats["count"] += 1
        stats["pathwise_l1"] += abs(z)
        stats["pathwise_l2_sq"] += abs(z) ** 2
        stats["max_abs"] = max(float(stats["max_abs"]), abs(z))

    def snapshot(self) -> dict[str, Any]:
        summary = {}
        for site, raw in sorted(self._by_site.items()):
            summary[site] = {
                "coefficient_count": int(raw["count"]),
                "pathwise_l1": float(raw["pathwise_l1"]),
                "pathwise_l2": math.sqrt(float(raw["pathwise_l2_sq"])),
                "max_abs": float(raw["max_abs"]),
            }
        return {
            "schema": "bqg-pruning-audit-v1",
            "event_count": len(self.events),
            "events": list(self.events),
            "by_site": summary,
            "semantics": (
                "Records actual discarded aggregate coefficients at instrumented "
                "sites. Pathwise totals are not coherent output norms and do not "
                "include uninstrumented sites."
            ),
        }


class capture:
    """Context manager enabling an isolated collector in the current context."""
    def __init__(self) -> None:
        self.collector = PruningCollector()
        self._token = None

    def __enter__(self) -> PruningCollector:
        self._token = _ACTIVE.set(self.collector)
        return self.collector

    def __exit__(self, exc_type, exc, tb) -> bool:
        if self._token is not None:
            _ACTIVE.reset(self._token)
        return False


def record_threshold(site_id: str, key: Any, value: complex, threshold: float) -> None:
    collector = _ACTIVE.get()
    if collector is not None and abs(value) <= threshold:
        collector.record(site_id, key, value, threshold, kind="coefficient_threshold")


def record_merge(site_id: str, key: Any, prior: complex, incoming: complex,
                 result: complex, threshold: float) -> None:
    collector = _ACTIVE.get()
    if collector is not None and abs(result) <= threshold:
        collector.record(
            site_id, key, result, threshold, kind="incremental_merge_result",
            prior=prior, incoming=incoming,
        )
