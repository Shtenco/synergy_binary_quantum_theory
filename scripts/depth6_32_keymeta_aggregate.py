#!/usr/bin/env python3
"""Fail-closed serialization helpers for the depth-6 [3,2] KEYMETA union gate."""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any

QKEY_ENCODING = 'BQG_QKEY_JSON_V1'


def _tag_q(value: Any):
    if value is None:
        return ['none', None]
    if isinstance(value, bool):
        return ['bool', value]
    if isinstance(value, int):
        return ['int', value]
    if isinstance(value, float):
        if value != value or value in (float('inf'), float('-inf')):
            raise RuntimeError('non-finite q-key float is not serializable')
        return ['float', value]
    if isinstance(value, str):
        return ['str', value]
    if isinstance(value, tuple):
        return ['tuple', [_tag_q(x) for x in value]]
    if isinstance(value, list):
        return ['list', [_tag_q(x) for x in value]]
    raise RuntimeError(f'unsupported q-key type: {type(value).__name__}')


def q_key_to_text(q: Any) -> str:
    """Canonical, collision-resistant JSON text for one Python q key."""
    return json.dumps(_tag_q(q), sort_keys=False, separators=(',', ':'), ensure_ascii=True)


def _encode_count_map(mapping: dict) -> dict[str, int]:
    if not isinstance(mapping, dict):
        raise RuntimeError('q count map must be a dict')
    out: dict[str, int] = {}
    for q, n0 in mapping.items():
        n = int(n0)
        if n < 0:
            raise RuntimeError(f'negative q count for {q!r}')
        key = q_key_to_text(q)
        if key in out:
            raise RuntimeError(f'q-key encoding collision for {q!r}')
        out[key] = n
    return out


def json_safe_gate(gate: dict) -> dict:
    if not isinstance(gate, dict):
        raise RuntimeError('gate must be a dict')
    out = dict(gate)
    out['q_key_encoding'] = QKEY_ENCODING
    out['q_block_occupancy'] = _encode_count_map(gate.get('q_block_occupancy', {}))
    out['q_row_totals'] = _encode_count_map(gate.get('q_row_totals', {}))
    if out.get('rank_certified') is True:
        raise RuntimeError('KEYMETA union serialization must not claim rank')
    if out.get('numerical_closure_claimed') is True:
        raise RuntimeError('KEYMETA union serialization must not claim numerical closure')
    return out


def write_json_gate(path: Path, gate: dict) -> dict:
    """Atomically persist a JSON-safe gate and return the serialized payload."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json_safe_gate(gate)
    raw = (json.dumps(payload, indent=2, sort_keys=True) + '\n').encode('utf-8')
    fd, tmp = tempfile.mkstemp(prefix=path.name + '.', suffix='.tmp', dir=str(path.parent))
    try:
        with os.fdopen(fd, 'wb') as f:
            f.write(raw)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass
        raise
    return payload
