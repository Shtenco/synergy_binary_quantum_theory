#!/usr/bin/env python3
"""KEYMETA-first numerical replay for the frozen depth-6 [3,2] assignments.

Scientific/claim boundary
-------------------------
This module consumes an already-frozen assignment record ``(orbit_id, rep, m,
coord_dim, provenance...)``.  It does not enumerate the depth-6 shell, rebuild
S5 orbits, recompute multiplicities, or construct a fresh structural/Jucys
assignment identity.

The only numerical operations performed here are:

1. construct a numerical coordinate basis of the *frozen* multiplicity ``m``
   through ``depth6_32_frozen_m_basis.build_frozen_m_basis``;
2. call the immutable Stage-A ``master_map`` on that basis;
3. record the actual non-zero q-row occupancy and dimensions;
4. persist that small KEYMETA payload before any optional heavy matrix payload.

KEYMETA completeness is only a prerequisite for global peeling/SVD.  Neither a
single block nor a complete KEYMETA union is a rank certificate by itself.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import pickle
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

IRREP_KEY = '32'
IRREP_LABEL = '[3,2]'
FROZEN_ENGINE_BLOB_SHA = 'cbbaa80f6b6b61d353458db1fa1f03f824730532'


def _json_hash(payload: dict) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode('utf-8')
    return hashlib.sha256(raw).hexdigest()


def _frozen_identity(record: dict) -> dict:
    """Return the immutable numerical identity used by this replay.

    Deliberately no structural reconstruction is available here.  Optional
    provenance is copied only when it was already present in the input record.
    """
    if not isinstance(record, dict):
        raise RuntimeError('frozen assignment record must be an object')
    missing = [k for k in ('orbit_id', 'rep', 'm', 'coord_dim') if k not in record]
    if missing:
        raise RuntimeError(f'frozen assignment record missing fields: {missing}')

    rep = tuple(int(x) for x in record['rep'])
    if len(rep) != 10:
        raise RuntimeError('frozen assignment rep must contain 10 doubled-spin labels')
    m = int(record['m'])
    d = int(record['coord_dim'])
    if m <= 0:
        raise RuntimeError('frozen assignment m must be positive')
    if d <= 0:
        raise RuntimeError('frozen assignment coord_dim must be positive')

    out = {
        'orbit_id': int(record['orbit_id']),
        'rep': list(rep),
        'm': m,
        'coord_dim': d,
    }
    for key in (
        'source_shard',
        'source_keymeta_payload_sha256',
        'source_raw_run_id',
        'source_raw_artifact_id',
        'source_raw_artifact_digest',
        'source_tar_sha256',
        'engine_git_blob_sha',
    ):
        if key in record and record[key] is not None:
            out[key] = record[key]
    return out


def _validate_basis_contract(W: np.ndarray, identity: dict, basis_cert: dict) -> np.ndarray:
    W = np.asarray(W, dtype=complex)
    expected = (int(identity['coord_dim']), int(identity['m']))
    if W.ndim != 2 or tuple(W.shape) != expected:
        raise RuntimeError(f'bad frozen basis shape: got={tuple(W.shape)} expected={expected}')
    if not isinstance(basis_cert, dict):
        raise RuntimeError('frozen basis certificate must be an object')
    if basis_cert.get('multiplicity_recomputed') is True:
        raise RuntimeError('frozen basis contract violated: multiplicity was recomputed')
    for forbidden in ('shell_recomputed', 'orbit_list_recomputed', 'structural_closure_claimed'):
        if basis_cert.get(forbidden) is True:
            raise RuntimeError(f'frozen basis contract violated: {forbidden}=true')
    return W


def _compute_block_material(record: dict, selector_backend, master_backend, basis_backend,
                            *, tol: float = 1e-8):
    identity = _frozen_identity(record)
    rep = tuple(identity['rep'])
    m = int(identity['m'])
    d = int(identity['coord_dim'])

    W, basis_cert = basis_backend.build_frozen_m_basis(
        rep, IRREP_KEY, m, d, selector_backend, tol=tol
    )
    W = _validate_basis_contract(W, identity, basis_cert)

    result = master_backend.master_map(rep, IRREP_KEY, W)
    if not isinstance(result, tuple) or len(result) != 2:
        raise RuntimeError('master_map must return (map, basis)')
    master_map, returned_W = result
    if not isinstance(master_map, dict):
        raise RuntimeError('master_map payload must be a dict')
    returned_W = np.asarray(returned_W, dtype=complex)
    if returned_W.ndim != 2 or returned_W.shape[1] != m:
        raise RuntimeError(
            f'bad master-map returned basis shape: got={tuple(returned_W.shape)} expected_columns={m}'
        )

    q_rows: dict[Any, int] = {}
    rows = 0
    for q, matrix in master_map.items():
        M = np.asarray(matrix)
        if M.ndim != 2 or M.shape[1] != m:
            raise RuntimeError(
                f'bad master-map shape for orbit {identity["orbit_id"]}: '
                f'q={q!r} shape={tuple(M.shape)} expected_columns={m}'
            )
        nrows = int(M.shape[0])
        if nrows <= 0:
            raise RuntimeError(
                f'bad master-map shape for orbit {identity["orbit_id"]}: '
                f'q={q!r} has non-positive row count {nrows}'
            )
        q_rows[q] = nrows
        rows += nrows

    basis_summary = {
        'status': basis_cert.get('status'),
        'max_residual': basis_cert.get('max_residual'),
        'multiplicity_recomputed': bool(basis_cert.get('multiplicity_recomputed', False)),
        'shell_recomputed': bool(basis_cert.get('shell_recomputed', False)),
        'orbit_list_recomputed': bool(basis_cert.get('orbit_list_recomputed', False)),
    }
    metadata = {
        'schema_version': 1,
        'kind': 'BQG_DEPTH6_32_KEYMETA_FIRST_BLOCK',
        'irrep': IRREP_LABEL,
        'irrep_key': IRREP_KEY,
        **identity,
        'frozen_assignment_sha256': _json_hash(identity),
        'q_rows': q_rows,
        'q_key_count': len(q_rows),
        'rows': rows,
        'basis_witness': basis_summary,
        'rank_certified': False,
        'global_uniqueness_certified': False,
        'numerical_closure_claimed': False,
        'structural_recompute_performed': False,
        'claim_boundary': (
            'Actual non-zero master-map q occupancy for one frozen [3,2] assignment. '
            'No global uniqueness, peeling, rank, kernel, or numerical closure claim.'
        ),
    }
    return metadata, master_map, W


def compute_block_keymeta(record: dict, selector_backend, master_backend, basis_backend,
                          *, tol: float = 1e-8) -> dict:
    """Compute only the small per-orbit KEYMETA view.

    The master matrices are used transiently to discover their actual q keys and
    row counts, then discarded by this public API.
    """
    metadata, _, _ = _compute_block_material(
        record, selector_backend, master_backend, basis_backend, tol=tol
    )
    return metadata


def validate_keymeta_union(blocks, *, target_blocks: int, target_columns: int) -> dict:
    """Validate an in-memory union of per-orbit KEYMETA records fail-closed.

    Exact target coverage authorizes the *next* global q-occupancy/peeling step;
    it does not certify that any column peels or that the master operator has
    full rank.
    """
    target_blocks = int(target_blocks)
    target_columns = int(target_columns)
    if target_blocks <= 0 or target_columns <= 0:
        raise RuntimeError('target_blocks and target_columns must be positive')

    seen: set[int] = set()
    columns = 0
    q_block_occupancy: dict[Any, int] = {}
    q_row_totals: dict[Any, int] = {}
    for block in blocks:
        oid = int(block['orbit_id'])
        if oid in seen:
            raise RuntimeError(f'duplicate orbit {oid}')
        seen.add(oid)
        m = int(block['m'])
        if m <= 0:
            raise RuntimeError(f'orbit {oid}: non-positive m')
        columns += m
        q_rows = block.get('q_rows', {})
        if not isinstance(q_rows, dict):
            raise RuntimeError(f'orbit {oid}: q_rows must be a dict')
        for q, n0 in q_rows.items():
            n = int(n0)
            if n <= 0:
                raise RuntimeError(f'orbit {oid}: q={q!r} has non-positive rows')
            q_block_occupancy[q] = q_block_occupancy.get(q, 0) + 1
            q_row_totals[q] = q_row_totals.get(q, 0) + n

    block_count = len(seen)
    if block_count > target_blocks:
        raise RuntimeError(f'KEYMETA block count exceeds target: {block_count}>{target_blocks}')
    if columns > target_columns:
        raise RuntimeError(f'KEYMETA column count exceeds target: {columns}>{target_columns}')
    if block_count == target_blocks and columns != target_columns:
        raise RuntimeError(
            f'complete KEYMETA block count but target column mismatch {columns}!={target_columns}'
        )

    complete = block_count == target_blocks and columns == target_columns
    return {
        'schema_version': 1,
        'kind': 'BQG_DEPTH6_32_KEYMETA_UNION_GATE',
        'irrep': IRREP_LABEL,
        'blocks': block_count,
        'columns': columns,
        'target_blocks': target_blocks,
        'target_columns': target_columns,
        'missing_blocks': target_blocks - block_count,
        'missing_columns': target_columns - columns,
        'q_key_count': len(q_block_occupancy),
        'q_block_occupancy': q_block_occupancy,
        'q_row_totals': q_row_totals,
        'exact_keymeta_coverage': bool(complete),
        'global_peeling_allowed': bool(complete),
        'rank_certified': False,
        'numerical_closure_claimed': False,
        'claim_boundary': (
            'Complete KEYMETA coverage only unlocks global actual-q peeling/SVD. '
            'It is not itself a rank or closure certificate.'
        ),
    }


def _atomic_pickle_gzip(path: Path, payload) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = gzip.compress(pickle.dumps(payload, protocol=5), compresslevel=6)
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


def compute_and_persist_block(record: dict, selector_backend, master_backend, basis_backend,
                              *, keymeta_path: Path, raw_path: Path | None = None,
                              tol: float = 1e-8) -> dict:
    """Compute one frozen block and enforce KEYMETA-before-RAW persistence order."""
    metadata, master_map, W = _compute_block_material(
        record, selector_backend, master_backend, basis_backend, tol=tol
    )

    # Scientific/recovery invariant: the compact actual-q identity becomes
    # durable before any heavy matrix artifact is attempted.
    _atomic_pickle_gzip(Path(keymeta_path), metadata)

    if raw_path is not None:
        raw_payload = {
            'schema_version': 1,
            'kind': 'BQG_DEPTH6_32_FROZEN_MASTER_BLOCK_RAW',
            'irrep': IRREP_LABEL,
            'orbit_id': metadata['orbit_id'],
            'm': metadata['m'],
            'coord_dim': metadata['coord_dim'],
            'frozen_assignment_sha256': metadata['frozen_assignment_sha256'],
            'master_map': master_map,
            'basis': W,
            'rank_certified': False,
            'numerical_closure_claimed': False,
        }
        _atomic_pickle_gzip(Path(raw_path), raw_payload)
    return metadata


def _load_cli_record(path: Path) -> dict:
    payload = json.loads(Path(path).read_text())
    # Production CLI accepts a single frozen record, or a one-record wrapper.
    # It deliberately does not derive missing assignments from structural data.
    if isinstance(payload, dict) and 'records' in payload:
        records = payload['records']
        if not isinstance(records, list) or len(records) != 1:
            raise RuntimeError('assignment input wrapper must contain exactly one frozen record')
        payload = records[0]
    if not isinstance(payload, dict):
        raise RuntimeError('assignment input must be a frozen assignment object')
    return payload


def _require_production_provenance(record: dict) -> None:
    required = (
        'source_shard',
        'source_keymeta_payload_sha256',
        'source_raw_run_id',
        'source_raw_artifact_id',
        'source_raw_artifact_digest',
        'source_tar_sha256',
        'engine_git_blob_sha',
    )
    missing = [k for k in required if record.get(k) in (None, '')]
    if missing:
        raise RuntimeError(f'production frozen assignment missing provenance: {missing}')
    if str(record['engine_git_blob_sha']) != FROZEN_ENGINE_BLOB_SHA:
        raise RuntimeError(
            'frozen engine blob mismatch: '
            f"record={record['engine_git_blob_sha']} expected={FROZEN_ENGINE_BLOB_SHA}"
        )


def main() -> int:
    ap = argparse.ArgumentParser(description='Depth-6 [3,2] frozen-assignment KEYMETA-first replay')
    ap.add_argument('--assignment', type=Path, required=True,
                    help='JSON containing exactly one frozen assignment record')
    ap.add_argument('--keymeta-out', type=Path, required=True,
                    help='small gzip-pickle KEYMETA output; always written first')
    ap.add_argument('--raw-out', type=Path,
                    help='optional heavy gzip-pickle master-map/basis output written only after KEYMETA')
    ap.add_argument('--tol', type=float, default=1e-8)
    args = ap.parse_args()

    record = _load_cli_record(args.assignment)
    _require_production_provenance(record)

    # These imports are intentionally local.  CI/runtime places the immutable
    # frozen Stage-A bundle on PYTHONPATH.  Importing this module alone therefore
    # never triggers or requires structural reconstruction machinery.
    import bqg_depth6_generic_jucys_selector as selector_backend
    import bqg_depth6_generic_master_row as master_backend
    import depth6_32_frozen_m_basis as basis_backend

    metadata = compute_and_persist_block(
        record,
        selector_backend,
        master_backend,
        basis_backend,
        keymeta_path=args.keymeta_out,
        raw_path=args.raw_out,
        tol=args.tol,
    )
    print(
        json.dumps(
            {
                'status': 'KEYMETA_PERSISTED_BEFORE_RAW',
                'orbit_id': metadata['orbit_id'],
                'm': metadata['m'],
                'q_key_count': metadata['q_key_count'],
                'rows': metadata['rows'],
                'keymeta_out': str(args.keymeta_out),
                'raw_out': str(args.raw_out) if args.raw_out is not None else None,
                'rank_certified': False,
                'numerical_closure_claimed': False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
