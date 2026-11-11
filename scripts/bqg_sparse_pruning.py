#!/usr/bin/env python3
"""Deterministic coherent sparse-amplitude aggregation and pruning ledger.

This module is intentionally independent of the K5 producer. It supplies the
small algebraic primitive that the producer can adopt once its policy is
threaded through every operator/cache layer.

No contribution is dropped while accumulating. Contributions are grouped by
canonical key, sorted deterministically, and summed with math.fsum separately
on real and imaginary parts. Pruning is a separate explicit operation.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from hashlib import sha256
import json
import math
from typing import Hashable, Iterable, Mapping, Sequence


@dataclass(frozen=True)
class PruningPolicy:
    mode: str = "no_prune"
    threshold: float = 0.0
    stage: str = "unspecified"

    def __post_init__(self) -> None:
        if self.mode not in {"no_prune", "aggregate_then_prune"}:
            raise ValueError("mode must be no_prune or aggregate_then_prune")
        if not math.isfinite(self.threshold) or self.threshold < 0:
            raise ValueError("threshold must be finite and nonnegative")
        if self.mode == "no_prune" and self.threshold != 0:
            raise ValueError("no_prune requires threshold=0")


@dataclass(frozen=True)
class PruningLedgerEntry:
    stage: str
    mode: str
    threshold: float
    input_keys: int
    kept_keys: int
    discarded_keys: int
    max_discarded_amplitude: float
    discarded_l1: float
    discarded_l2: float
    discarded_vector_sha256: str
    coherent_aggregation_complete: bool


def _stable_key(key: Hashable) -> str:
    """Stable representation for supported canonical tuple/scalar keys."""
    return repr(key)


def _stable_complex_sum(values: Sequence[complex]) -> complex:
    ordered = sorted((complex(v) for v in values), key=lambda z: (z.real, z.imag))
    re = math.fsum(z.real for z in ordered)
    im = math.fsum(z.imag for z in ordered)
    return complex(re, im)


def coherent_aggregate(
    contributions: Iterable[tuple[Hashable, complex]],
) -> dict[Hashable, complex]:
    """Coherently sum every contribution per key without thresholding."""
    buckets: dict[Hashable, list[complex]] = {}
    for key, raw in contributions:
        value = complex(raw)
        if not (math.isfinite(value.real) and math.isfinite(value.imag)):
            raise ValueError(f"non-finite amplitude for key {key!r}")
        buckets.setdefault(key, []).append(value)
    out: dict[Hashable, complex] = {}
    for key in sorted(buckets, key=_stable_key):
        value = _stable_complex_sum(buckets[key])
        if value != 0j:
            out[key] = value
    return out


def aggregate_mappings(
    mappings: Sequence[Mapping[Hashable, complex]],
    scales: Sequence[complex] | None = None,
) -> dict[Hashable, complex]:
    """Coherently combine sparse maps, optionally applying one scale per map."""
    if scales is None:
        scales = [1.0 + 0.0j] * len(mappings)
    if len(scales) != len(mappings):
        raise ValueError("scales and mappings must have equal lengths")
    contributions = []
    for mapping, scale in zip(mappings, scales):
        s = complex(scale)
        for key, value in mapping.items():
            contributions.append((key, s * complex(value)))
    return coherent_aggregate(contributions)


def _discard_digest(discarded: Mapping[Hashable, complex]) -> str:
    records = [
        (_stable_key(k), float(complex(v).real).hex(), float(complex(v).imag).hex())
        for k, v in sorted(discarded.items(), key=lambda item: _stable_key(item[0]))
    ]
    payload = json.dumps(records, ensure_ascii=True, separators=(",", ":"))
    return sha256(payload.encode("utf-8")).hexdigest()


def apply_pruning(
    aggregated: Mapping[Hashable, complex],
    policy: PruningPolicy,
    *,
    coherent_aggregation_complete: bool = True,
) -> tuple[dict[Hashable, complex], PruningLedgerEntry]:
    """Apply pruning only to an already aggregated canonical-key mapping."""
    if not coherent_aggregation_complete:
        raise ValueError("pruning is forbidden before coherent key aggregation")
    source = {k: complex(v) for k, v in aggregated.items() if complex(v) != 0j}
    if policy.mode == "no_prune":
        kept, discarded = dict(source), {}
    elif policy.mode == "aggregate_then_prune":
        kept = {k: v for k, v in source.items() if abs(v) > policy.threshold}
        discarded = {k: v for k, v in source.items() if abs(v) <= policy.threshold}
    else:  # defensive; __post_init__ already validates
        raise ValueError(f"unsupported policy mode {policy.mode!r}")

    vals = list(discarded.values())
    entry = PruningLedgerEntry(
        stage=policy.stage,
        mode=policy.mode,
        threshold=policy.threshold,
        input_keys=len(source),
        kept_keys=len(kept),
        discarded_keys=len(discarded),
        max_discarded_amplitude=max((abs(v) for v in vals), default=0.0),
        discarded_l1=math.fsum(abs(v) for v in vals),
        discarded_l2=math.sqrt(math.fsum(abs(v) ** 2 for v in vals)),
        discarded_vector_sha256=_discard_digest(discarded),
        coherent_aggregation_complete=True,
    )
    return kept, entry


def verify_ledger_identity(
    unpruned: Mapping[Hashable, complex],
    pruned: Mapping[Hashable, complex],
    entry: PruningLedgerEntry,
    *,
    atol: float = 0.0,
) -> bool:
    """Verify that unpruned-pruned equals the discarded vector encoded by ledger.

    The digest is a compact identity check; this function additionally checks
    the norm/count invariants available from the two supplied maps.
    """
    keys = set(unpruned) | set(pruned)
    delta = {k: complex(unpruned.get(k, 0j)) - complex(pruned.get(k, 0j))
             for k in keys}
    delta = {k: v for k, v in delta.items() if abs(v) > atol}
    return (
        len(delta) == entry.discarded_keys
        and math.isclose(math.fsum(abs(v) for v in delta.values()),
                         entry.discarded_l1, rel_tol=1e-12, abs_tol=atol)
        and math.isclose(math.sqrt(math.fsum(abs(v)**2 for v in delta.values())),
                         entry.discarded_l2, rel_tol=1e-12, abs_tol=atol)
        and _discard_digest(delta) == entry.discarded_vector_sha256
    )


def telescoping_error_bound(
    discarded_norms: Sequence[float],
    downstream_operator_norm_products: Sequence[float],
) -> float:
    """Return sum_s ||A_{n-1}...A_{s+1}|| * ||delta_s||."""
    if len(discarded_norms) != len(downstream_operator_norm_products):
        raise ValueError("one downstream norm product is required per discarded vector")
    terms = []
    for d, a in zip(discarded_norms, downstream_operator_norm_products):
        if not (math.isfinite(d) and math.isfinite(a)) or d < 0 or a < 0:
            raise ValueError("norm inputs must be finite and nonnegative")
        terms.append(d * a)
    return math.fsum(terms)
