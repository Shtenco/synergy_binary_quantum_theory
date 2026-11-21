#!/usr/bin/env python3
"""Build and verify immutable depth-6 persisted assignment ledgers.

The module accepts explicit persisted assignment records only.  It contains no
shell, orbit, multiplicity or Jucys reconstruction path.
"""
from __future__ import annotations

import hashlib
import json

REQUIRED_PROVENANCE = (
    'source_shard',
    'source_keymeta_payload_sha256',
    'source_raw_run_id',
    'source_raw_artifact_id',
    'source_raw_artifact_digest',
    'source_tar_sha256',
    'engine_git_blob_sha',
)


def _missing_value(record: dict, key: str) -> bool:
    if key not in record:
        return True
    value = record[key]
    return value is None or (isinstance(value, str) and value == '')


def normalize_assignment_record(record: dict) -> dict:
    missing = [k for k in REQUIRED_PROVENANCE if _missing_value(record, k)]
    if missing:
        raise RuntimeError(f'missing provenance fields: {missing}')
    if 'orbit_id' not in record or 'rep' not in record or 'm' not in record:
        raise RuntimeError('assignment record requires orbit_id, rep, m')
    rep = [int(x) for x in record['rep']]
    if len(rep) != 10:
        raise RuntimeError(f"orbit {record['orbit_id']}: rep must contain 10 doubled-spin labels")
    m = int(record['m'])
    if m <= 0:
        raise RuntimeError(f"orbit {record['orbit_id']}: m must be positive")
    source_shard = int(record['source_shard'])
    if source_shard < 0:
        raise RuntimeError(f"orbit {record['orbit_id']}: source_shard must be non-negative")
    out = {
        'orbit_id': int(record['orbit_id']),
        'rep': rep,
        'm': m,
        'coord_dim': int(record.get('coord_dim', 0)),
        'source_shard': source_shard,
        'source_keymeta_payload_sha256': str(record['source_keymeta_payload_sha256']),
        'source_raw_run_id': str(record['source_raw_run_id']),
        'source_raw_artifact_id': str(record['source_raw_artifact_id']),
        'source_raw_artifact_digest': str(record['source_raw_artifact_digest']),
        'source_tar_sha256': str(record['source_tar_sha256']),
        'engine_git_blob_sha': str(record['engine_git_blob_sha']),
    }
    for optional, caster in (
        ('observed_sec', float),
        ('observed_rows', int),
        ('observed_row_blocks', int),
    ):
        if optional in record and record[optional] is not None:
            out[optional] = caster(record[optional])
    return out


def _hash_material(records: list[dict], target_blocks: int, target_columns: int) -> str:
    payload = {
        'irrep': '[3,2]',
        'target_blocks': int(target_blocks),
        'target_columns': int(target_columns),
        'records': records,
    }
    raw = json.dumps(payload, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode('utf-8')
    return hashlib.sha256(raw).hexdigest()


def build_assignment_ledger(records, *, target_blocks: int, target_columns: int) -> dict:
    normalized = sorted((normalize_assignment_record(r) for r in records), key=lambda r: r['orbit_id'])
    seen: set[int] = set()
    for r in normalized:
        oid = r['orbit_id']
        if oid in seen:
            raise RuntimeError(f'duplicate orbit_id {oid}')
        seen.add(oid)

    persisted_blocks = len(normalized)
    persisted_columns = sum(r['m'] for r in normalized)
    if persisted_blocks > target_blocks:
        raise RuntimeError(f'persisted block count exceeds target: {persisted_blocks}>{target_blocks}')
    if persisted_columns > target_columns:
        raise RuntimeError(f'persisted column count exceeds target: {persisted_columns}>{target_columns}')
    if persisted_blocks == target_blocks and persisted_columns != target_columns:
        raise RuntimeError(
            'complete block count but target column mismatch '
            f'{persisted_columns}!={target_columns}'
        )

    complete = persisted_blocks == target_blocks and persisted_columns == target_columns
    return {
        'schema_version': 1,
        'kind': 'BQG_DEPTH6_32_PERSISTED_ASSIGNMENT_LEDGER',
        'irrep': '[3,2]',
        'target_blocks': int(target_blocks),
        'target_columns': int(target_columns),
        'records': normalized,
        'coverage': {
            'status': 'COMPLETE_PERSISTED_ASSIGNMENT_LEDGER' if complete else 'INCOMPLETE_PERSISTED_ASSIGNMENT_LEDGER',
            'persisted_blocks': persisted_blocks,
            'persisted_columns': persisted_columns,
            'target_blocks': int(target_blocks),
            'target_columns': int(target_columns),
            'missing_blocks': int(target_blocks) - persisted_blocks,
            'missing_columns': int(target_columns) - persisted_columns,
            'retry_allowed': bool(complete),
            'rank_certified': False,
        },
        'ledger_sha256': _hash_material(normalized, target_blocks, target_columns),
        'claim_boundary': (
            'Persisted frozen assignment identity only. retry_allowed means assignment completeness, '
            'not actual-q coverage, rank certification, or numerical closure.'
        ),
    }


def verify_assignment_ledger(ledger: dict) -> dict:
    rebuilt = build_assignment_ledger(
        ledger.get('records', []),
        target_blocks=int(ledger.get('target_blocks', -1)),
        target_columns=int(ledger.get('target_columns', -1)),
    )
    if ledger.get('kind') != rebuilt['kind'] or ledger.get('irrep') != '[3,2]':
        raise RuntimeError('wrong assignment ledger kind/irrep')
    if ledger.get('ledger_sha256') != rebuilt['ledger_sha256']:
        raise RuntimeError('assignment ledger SHA256 mismatch')
    supplied = ledger.get('coverage', {})
    expected = rebuilt['coverage']
    for key in ('status','persisted_blocks','persisted_columns','target_blocks','target_columns','missing_blocks','missing_columns','retry_allowed','rank_certified'):
        if supplied.get(key) != expected.get(key):
            raise RuntimeError(f'assignment ledger coverage mismatch for {key}')
    return rebuilt
