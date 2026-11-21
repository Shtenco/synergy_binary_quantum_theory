#!/usr/bin/env python3
"""Integrity checks for the complete depth-6 [3,2] assignment identity snapshot.

The snapshot is a deterministic identity re-materialization from the immutable
Stage-A engine.  It is useful for scheduling already-frozen orbit identities,
but it is deliberately not actual-q evidence, a numerical rank witness, or a
license to reconstruct the structural/Jucys proof.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

LEDGER_NAME = 'BQG_DEPTH6_32_REMATERIALIZED_ASSIGNMENT_LEDGER_2026-10-02.json'
GATE_NAME = 'BQG_DEPTH6_32_IDENTITY_REMATERIALIZATION_GATE_2026-10-02.json'


def _require(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(msg)


def _records_sha(records: list[dict]) -> str:
    raw = json.dumps(records, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode('utf-8')
    return hashlib.sha256(raw).hexdigest()


def verify_snapshot(root: Path, frontier: dict) -> dict:
    meta = frontier.get('rematerialized_assignment_ledger', {})
    _require(meta.get('path') == LEDGER_NAME, 'rematerialized identity ledger path mismatch')
    _require(meta.get('gate_path') == GATE_NAME, 'rematerialized identity gate path mismatch')
    lp, gp = root / LEDGER_NAME, root / GATE_NAME
    _require(lp.is_file(), f'missing rematerialized identity ledger: {LEDGER_NAME}')
    _require(gp.is_file(), f'missing rematerialized identity gate: {GATE_NAME}')
    ledger = json.loads(lp.read_text())
    gate = json.loads(gp.read_text())

    _require(ledger.get('kind') == 'BQG_DEPTH6_32_REMATERIALIZED_ASSIGNMENT_LEDGER', 'wrong rematerialized identity kind')
    _require(ledger.get('status') == 'COMPLETE_IDENTITY_REMATERIALIZATION_ONLY', 'rematerialized identity status mismatch')
    _require(ledger.get('target_blocks') == 2755, 'rematerialized identity block target mismatch')
    _require(ledger.get('target_columns') == 130903, 'rematerialized identity column target mismatch')
    _require(ledger.get('shell_states') == 264962, 'rematerialized shell count mismatch')
    _require(ledger.get('s5_orbits') == 2757, 'rematerialized S5 orbit count mismatch')
    records = ledger.get('records', [])
    _require(len(records) == 2755, 'rematerialized identity must contain 2755 records')
    ids = [int(r['orbit_id']) for r in records]
    _require(len(set(ids)) == 2755, 'rematerialized identity contains duplicate orbit ids')
    _require(sum(int(r['m']) for r in records) == 130903, 'rematerialized identity column sum mismatch')
    for r in records:
        _require(len(r.get('rep', [])) == 10, f"rematerialized identity orbit {r.get('orbit_id')} has invalid rep")
        _require(int(r.get('m', 0)) > 0 and int(r.get('coord_dim', 0)) > 0, 'rematerialized identity contains non-positive dimensions')
    sha = _records_sha(records)
    _require(ledger.get('identity_sha256') == sha, 'rematerialized identity SHA mismatch')

    _require(gate.get('kind') == 'BQG_DEPTH6_32_IDENTITY_REMATERIALIZATION_GATE', 'wrong rematerialization gate kind')
    _require(gate.get('status') == 'PASS_IDENTITY_REMATERIALIZATION_ONLY', 'rematerialization gate did not pass')
    _require(gate.get('blocks') == 2755 and gate.get('columns') == 130903, 'rematerialization gate target mismatch')
    _require(gate.get('shell_states') == 264962 and gate.get('s5_orbits') == 2757, 'rematerialization gate shell/orbit mismatch')
    _require(gate.get('persisted_identity_matches') == 15, 'rematerialization persisted-identity cross-check drifted')
    _require(gate.get('historical_batch_matches') == 102, 'rematerialization historical-batch cross-check drifted')
    for obj, label in ((ledger, 'ledger'), (gate, 'gate')):
        _require(obj.get('rank_certified') is False, f'rematerialized identity {label} must not claim rank')
        _require(obj.get('actual_q_coverage_certified') is False, f'rematerialized identity {label} must not claim actual-q coverage')
        _require(obj.get('numerical_closure_claimed') is False, f'rematerialized identity {label} must not claim numerical closure')

    _require(meta.get('status') == 'COMPLETE_IDENTITY_REMATERIALIZATION_ONLY', 'frontier rematerialized status mismatch')
    _require(meta.get('blocks') == 2755 and meta.get('columns') == 130903, 'frontier rematerialized target mismatch')
    _require(meta.get('shell_states') == 264962 and meta.get('s5_orbits') == 2757, 'frontier rematerialized shell/orbit mismatch')
    _require(meta.get('persisted_identity_matches') == 15 and meta.get('historical_batch_matches') == 102, 'frontier rematerialized cross-check mismatch')
    _require(meta.get('identity_sha256') == sha, 'frontier rematerialized identity SHA mismatch')
    _require(meta.get('rank_certified') is False and meta.get('actual_q_coverage_certified') is False, 'frontier rematerialized identity crossed numerical claim boundary')
    _require(meta.get('canonical_numerical_compute_allowed') is False, 'identity completeness cannot bypass replay-readiness gate')

    return {
        'assignment_identity_complete': True,
        'rematerialized_identity_blocks': 2755,
        'rematerialized_identity_columns': 130903,
        'rematerialized_identity_sha256': sha,
        'rematerialized_rank_certified': False,
        'rematerialized_actual_q_certified': False,
    }
