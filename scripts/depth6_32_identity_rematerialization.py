#!/usr/bin/env python3
"""Fail-closed acceptance gate for depth-6 [3,2] assignment identity re-materialization.

This module does not enumerate the shell, construct S5 orbits, compute
multiplicities, build Jucys selectors, or evaluate master maps.  It validates
an identity list produced by an explicitly quarantined deterministic replay of
the frozen Stage-A assignment constructor.

Passing this gate certifies assignment identity only.  It is not actual-q
coverage, rank evidence, a master-kernel certificate, or numerical closure.
"""
from __future__ import annotations

import hashlib
import json


def _identity(record: dict) -> dict:
    return {
        'orbit_id': int(record['orbit_id']),
        'rep': [int(x) for x in record['rep']],
        'm': int(record['m']),
        'coord_dim': int(record['coord_dim']),
        'proxy': int(record.get('proxy', int(record['coord_dim']) * int(record['m']))),
    }


def cost_balanced(records, nshards: int):
    """Reproduce the frozen Stage-A greedy cost balancing from explicit records."""
    if int(nshards) <= 0:
        raise RuntimeError('nshards must be positive')
    bins = [[] for _ in range(int(nshards))]
    loads = [0] * int(nshards)
    normalized = [_identity(r) for r in records]
    for rec in sorted(normalized, key=lambda x: (x['proxy'], x['m'], x['orbit_id']), reverse=True):
        j = min(range(int(nshards)), key=lambda b: (loads[b], b))
        bins[j].append(rec)
        loads[j] += rec['proxy']
    return bins, loads


def _hash_records(records: list[dict]) -> str:
    raw = json.dumps(records, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode('utf-8')
    return hashlib.sha256(raw).hexdigest()


def validate_rematerialized_identity(
    records,
    *,
    persisted_records,
    historical_batches,
    expected_shell_states: int,
    observed_shell_states: int,
    expected_s5_orbits: int,
    observed_s5_orbits: int,
    target_blocks: int,
    target_columns: int,
    nshards: int,
) -> dict:
    """Accept only an exact identity reconstruction consistent with frozen evidence."""
    if int(observed_shell_states) != int(expected_shell_states) or int(observed_s5_orbits) != int(expected_s5_orbits):
        raise RuntimeError(
            'frozen shell/orbit count mismatch: '
            f'shell={observed_shell_states}/{expected_shell_states} '
            f'orbits={observed_s5_orbits}/{expected_s5_orbits}'
        )

    normalized = sorted((_identity(r) for r in records), key=lambda x: x['orbit_id'])
    seen = set()
    for rec in normalized:
        oid = rec['orbit_id']
        if oid in seen:
            raise RuntimeError(f'duplicate rematerialized orbit_id {oid}')
        seen.add(oid)
        if len(rec['rep']) != 10:
            raise RuntimeError(f'orbit {oid}: rep must contain 10 doubled-spin labels')
        if rec['m'] <= 0 or rec['coord_dim'] <= 0 or rec['proxy'] <= 0:
            raise RuntimeError(f'orbit {oid}: non-positive identity dimensions')

    blocks = len(normalized)
    columns = sum(r['m'] for r in normalized)
    if blocks != int(target_blocks) or columns != int(target_columns):
        raise RuntimeError(
            f'rematerialized target mismatch blocks={blocks}/{target_blocks} '
            f'columns={columns}/{target_columns}'
        )

    by_id = {r['orbit_id']: r for r in normalized}
    persisted_matches = 0
    for raw in persisted_records:
        p = _identity(raw)
        got = by_id.get(p['orbit_id'])
        if got is None:
            raise RuntimeError(f"persisted identity mismatch: missing orbit {p['orbit_id']}")
        for key in ('rep', 'm', 'coord_dim'):
            if got[key] != p[key]:
                raise RuntimeError(
                    f"persisted identity mismatch orbit={p['orbit_id']} field={key}: "
                    f"{got[key]} != {p[key]}"
                )
        persisted_matches += 1

    bins, loads = cost_balanced(normalized, int(nshards))
    historical_matches = 0
    for shard0, observed in historical_batches.items():
        shard = int(shard0)
        if shard < 0 or shard >= int(nshards):
            raise RuntimeError(f'historical batch mismatch: invalid shard {shard}')
        expected = {
            'assigned_blocks': len(bins[shard]),
            'assigned_columns': sum(r['m'] for r in bins[shard]),
            'proxy': loads[shard],
        }
        for key in ('assigned_blocks', 'assigned_columns'):
            if int(observed[key]) != int(expected[key]):
                raise RuntimeError(
                    f'historical batch mismatch shard={shard} field={key}: '
                    f"{observed[key]} != {expected[key]}"
                )
        if 'proxy' in observed and observed['proxy'] is not None:
            if abs(float(observed['proxy']) - float(expected['proxy'])) > max(1e-9, 1e-12 * abs(float(expected['proxy']))):
                raise RuntimeError(
                    f'historical batch mismatch shard={shard} field=proxy: '
                    f"{observed['proxy']} != {expected['proxy']}"
                )
        historical_matches += 1

    return {
        'schema_version': 1,
        'kind': 'BQG_DEPTH6_32_IDENTITY_REMATERIALIZATION_GATE',
        'status': 'PASS_IDENTITY_REMATERIALIZATION_ONLY',
        'blocks': blocks,
        'columns': columns,
        'shell_states': int(observed_shell_states),
        's5_orbits': int(observed_s5_orbits),
        'persisted_identity_matches': persisted_matches,
        'historical_batch_matches': historical_matches,
        'nshards': int(nshards),
        'identity_sha256': _hash_records(normalized),
        'rank_certified': False,
        'actual_q_coverage_certified': False,
        'numerical_closure_claimed': False,
        'claim_boundary': (
            'Deterministic assignment identity reconstruction only. Passing this gate does not '
            'certify actual-q coverage, numerical rank, a master kernel, or closure of [3,2].'
        ),
    }
