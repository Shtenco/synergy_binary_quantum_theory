"""Opt-in, semantics-preserving audit of numerical sparse-state pruning.

The collector is context-local and disabled by default. Instrumented operators
continue to execute their original threshold decisions; when a coefficient is
discarded, the audit records the exact complex values that led to the decision.
"""
from __future__ import annotations

from collections import defaultdict
from contextvars import ContextVar
import functools
import math
from typing import Any, Callable

_ACTIVE: ContextVar["PruningCollector | None"] = ContextVar(
    "bqg_pruning_collector", default=None
)
_BYPASS: ContextVar[bool] = ContextVar("bqg_pruning_bypass", default=False)


class PruningCollector:
    def __init__(self) -> None:
        self.events: list[dict[str, Any]] = []
        self._by_site: dict[str, dict[str, float | int]] = defaultdict(
            lambda: {"count": 0, "pathwise_l1": 0.0, "pathwise_l2_sq": 0.0,
                     "max_abs": 0.0}
        )
        self._by_site_kind: dict[tuple[str, str], dict[str, float | int]] = defaultdict(
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
        kind_stats = self._by_site_kind[(str(site_id), str(kind))]
        kind_stats["count"] += 1
        kind_stats["pathwise_l1"] += abs(z)
        kind_stats["pathwise_l2_sq"] += abs(z) ** 2
        kind_stats["max_abs"] = max(float(kind_stats["max_abs"]), abs(z))

    def replay(self, rows: list[dict[str, Any]] | tuple[dict[str, Any], ...]) -> None:
        """Replay a cached call's exact audit rows into this collector.

        Rows are already encoded, so replay preserves key_repr and event order
        without rerunning the numerical operator body.
        """
        for source in rows:
            row = dict(source)
            self.events.append(row)
            site = str(row["site_id"])
            kind = str(row["kind"])
            magnitude = float(row["abs"])
            for stats in (self._by_site[site], self._by_site_kind[(site, kind)]):
                stats["count"] += 1
                stats["pathwise_l1"] += magnitude
                stats["pathwise_l2_sq"] += magnitude ** 2
                stats["max_abs"] = max(float(stats["max_abs"]), magnitude)

    def snapshot(self) -> dict[str, Any]:
        def encode(raw):
            return {
                "coefficient_count": int(raw["count"]),
                "pathwise_l1": float(raw["pathwise_l1"]),
                "pathwise_l2": math.sqrt(float(raw["pathwise_l2_sq"])),
                "max_abs": float(raw["max_abs"]),
            }
        summary = {
            site: encode(raw) for site, raw in sorted(self._by_site.items())
        }
        by_kind = {
            f"{site}::{kind}": encode(raw)
            for (site, kind), raw in sorted(self._by_site_kind.items())
        }
        return {
            "schema": "bqg-pruning-audit-v1",
            "event_count": len(self.events),
            "events": list(self.events),
            "by_site": summary,
            "by_site_and_kind": by_kind,
            "semantics": (
                "Event kinds distinguish production threshold removals, shadow-mode "
                "would-prune candidates, and exact merge cancellations. Pathwise "
                "totals are not coherent output norms and do not include "
                "uninstrumented sites."
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


class bypass:
    """Temporarily bypass only explicitly instrumented numerical amplitude thresholds."""
    def __init__(self) -> None:
        self._token = None

    def __enter__(self):
        self._token = _BYPASS.set(True)
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        if self._token is not None:
            _BYPASS.reset(self._token)
        return False


def keep(value: complex, threshold: float) -> bool:
    """Preserve the production rule unless an explicit shadow bypass is active."""
    z = complex(value)
    return z != 0j if _BYPASS.get() else abs(z) > threshold


def filter_state(site_id: str, state: dict, threshold: float) -> dict:
    """Apply one sparse-state amplitude threshold, or retain all exact nonzeros in shadow mode."""
    kept = {}
    for key, value in state.items():
        record_threshold(site_id, key, value, threshold)
        if keep(value, threshold):
            kept[key] = value
    return kept


def record_threshold(site_id: str, key: Any, value: complex, threshold: float) -> None:
    collector = _ACTIVE.get()
    z = complex(value)
    # Structural zeros are not amplitude-pruning losses. They remain absent
    # even in the shadow mode and are accounted only for merge cancellations.
    if collector is not None and z != 0j and abs(z) <= threshold:
        collector.record(
            site_id, key, z, threshold,
            kind="would_prune_bypassed" if _BYPASS.get() else "coefficient_threshold",
        )


def record_merge(site_id: str, key: Any, prior: complex, incoming: complex,
                 result: complex, threshold: float) -> None:
    collector = _ACTIVE.get()
    z = complex(result)
    p, inc = complex(prior), complex(incoming)
    exact_cancellation = z == 0j and p != 0j and inc != 0j
    subthreshold_nonzero = z != 0j and abs(z) <= threshold
    if collector is not None and (exact_cancellation or subthreshold_nonzero):
        kind = (
            "exact_cancellation" if exact_cancellation else
            "would_prune_bypassed" if _BYPASS.get() else
            "incremental_merge_result"
        )
        collector.record(
            site_id, key, z, threshold, kind=kind,
            prior=p, incoming=inc,
        )


def audit_lru_cache(maxsize=None):
    """LRU-cache a numerical function while replaying its audit trace on hits.

    The cached value is paired with the exact ordered PRUNE event rows emitted
    by the function and all nested instrumented calls. A private collector is
    active even when the caller is not auditing, so a later warm-cache audited
    call receives the same trace as a cold call. Context-local bypass semantics
    are inherited while computing the cache entry.
    """
    def decorate(function: Callable):
        @functools.lru_cache(maxsize=maxsize)
        def cached(bypass_enabled, *args, **kwargs):
            # Bypass changes numerical semantics and event classification, so
            # it is part of the cache key rather than ambient hidden state.
            local = PruningCollector()
            token = _ACTIVE.set(local)
            try:
                value = function(*args, **kwargs)
            finally:
                _ACTIVE.reset(token)
            return value, tuple(dict(row) for row in local.events)

        @functools.wraps(function)
        def wrapped(*args, **kwargs):
            value, rows = cached(_BYPASS.get(), *args, **kwargs)
            active = _ACTIVE.get()
            if active is not None:
                active.replay(rows)
            return value

        wrapped.cache_clear = cached.cache_clear
        wrapped.cache_info = cached.cache_info
        wrapped.cache_parameters = cached.cache_parameters
        return wrapped
    return decorate
