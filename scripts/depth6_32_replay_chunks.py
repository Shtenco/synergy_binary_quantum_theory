#!/usr/bin/env python3
"""Deterministically plan replay chunks over already-frozen [3,2] records.

This module never reconstructs structural identity or multiplicities.  It only
partitions frozen ledger records that are not already durable into smaller
recovery chunks using the persisted numerical proxy as a load measure.
"""
from __future__ import annotations


def plan_missing_chunks(records, completed_orbit_ids, *, target_proxy: int = 250000,
                        max_blocks: int | None = None) -> list[dict]:
    target_proxy = int(target_proxy)
    if target_proxy <= 0:
        raise RuntimeError('target_proxy must be positive')
    if max_blocks is not None and int(max_blocks) <= 0:
        raise RuntimeError('max_blocks must be positive when provided')
    if max_blocks is not None:
        max_blocks = int(max_blocks)

    by_id = {}
    for rec in records:
        oid = int(rec['orbit_id'])
        if oid in by_id:
            raise RuntimeError(f'duplicate orbit {oid}')
        proxy = int(rec.get('proxy', int(rec['m']) * int(rec['coord_dim'])))
        if proxy <= 0:
            raise RuntimeError(f'orbit {oid}: non-positive proxy')
        by_id[oid] = {
            'orbit_id': oid,
            'proxy': proxy,
            'm': int(rec['m']),
            'coord_dim': int(rec['coord_dim']),
        }

    completed = {int(x) for x in completed_orbit_ids}
    unknown = completed.difference(by_id)
    if unknown:
        raise RuntimeError(f'completed set contains unknown orbits: {sorted(unknown)[:10]}')

    missing = [v for oid, v in by_id.items() if oid not in completed]
    # First-fit decreasing with deterministic tie-break.  This keeps heavy
    # records isolated while making input ordering irrelevant.
    missing.sort(key=lambda r: (-r['proxy'], r['orbit_id']))

    bins: list[dict] = []
    for rec in missing:
        placed = False
        if rec['proxy'] <= target_proxy:
            for b in bins:
                if ((max_blocks is None or len(b['orbit_ids']) < max_blocks)
                    and b['proxy'] + rec['proxy'] <= target_proxy):
                    b['orbit_ids'].append(rec['orbit_id'])
                    b['proxy'] += rec['proxy']
                    b['columns'] += rec['m']
                    placed = True
                    break
        if not placed:
            bins.append({
                'orbit_ids': [rec['orbit_id']],
                'proxy': rec['proxy'],
                'columns': rec['m'],
            })

    for idx, b in enumerate(bins):
        b['chunk'] = idx
        b['blocks'] = len(b['orbit_ids'])
        b['orbit_ids'].sort()
        b['oversize_singleton'] = bool(b['proxy'] > target_proxy)
        if b['oversize_singleton'] and b['blocks'] != 1:
            raise RuntimeError('oversize recovery chunk must be singleton')

    planned = [oid for b in bins for oid in b['orbit_ids']]
    expected = sorted(r['orbit_id'] for r in missing)
    if sorted(planned) != expected or len(planned) != len(set(planned)):
        raise RuntimeError('recovery chunk plan is not an exact partition of missing orbits')
    return bins
