#!/usr/bin/env python3
"""Extract tiny support metadata from already-computed mixed master maps.

The output contains actual persisted numerical master-map output keys and row
counts, but no matrix coefficients. It is therefore suitable for constructing
a candidate numerical peeling schedule without recomputing shell/orbit/Jucys
or geometric structural support.

This metadata by itself is NOT a rank certificate and MUST NOT be interpreted
as numerical closure.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pickle
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'proof_tools'))
from safe_archive import safe_extract_tar
from safe_pickle import safe_load_numpy_pickle

import numpy as np

KIND = 'BQG_MIXED_MASTER_KEY_METADATA'
SCHEMA_VERSION = 1


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def extract_tar(src: Path, out: Path) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    with tarfile.open(src, 'r:gz') as tf:
        safe_extract_tar(tf, out)
    summaries = list(out.rglob('summary.json'))
    if len(summaries) != 1:
        raise RuntimeError(f'expected one summary.json, got {len(summaries)}')
    return summaries[0]


def extract_metadata(src: Path, work: Path, out: Path, provenance: dict) -> dict:
    summary_path = extract_tar(src, work)
    summary = json.loads(summary_path.read_text())
    blocks = {}
    total_q = 0
    total_rows = 0
    for rec in summary.get('records', []):
        if not rec.get('ok'):
            raise RuntimeError(f'non-ok numerical block in source shard: {rec}')
        i = int(rec['orbit_index'])
        m = int(rec['m'])
        matches = list(work.rglob(f'b{i:05d}.pkl'))
        if len(matches) != 1:
            raise RuntimeError(f'orbit {i}: expected one pickle, got {len(matches)}')
        with matches[0].open('rb') as f:
            mp = safe_load_numpy_pickle(f)
        q_rows = {}
        for q, a0 in mp.items():
            a = np.asarray(a0)
            if a.ndim != 2 or a.shape[1] != m:
                raise RuntimeError(f'orbit {i}, q={q!r}: bad map shape {a.shape}, m={m}')
            # The pickle output deliberately retains q as its exact Python tuple
            # rather than lossy stringification.
            q_rows[q] = int(a.shape[0])
            total_q += 1
            total_rows += int(a.shape[0])
        blocks[i] = {
            'm': m,
            'record': rec,
            'q_rows': q_rows,
        }
    payload = {
        'schema_version': SCHEMA_VERSION,
        'kind': KIND,
        'irrep': summary.get('irrep'),
        'shard': summary.get('shard'),
        'shards': summary.get('shards'),
        'active_blocks': summary.get('active_blocks'),
        'target_columns': summary.get('target_columns'),
        'assigned_blocks': summary.get('assigned_blocks'),
        'assigned_columns': summary.get('assigned_columns'),
        'source_tar': src.name,
        'source_tar_sha256': sha256(src),
        **provenance,
        'blocks': blocks,
        'stats': {
            'blocks': len(blocks),
            'columns': sum(int(x['m']) for x in blocks.values()),
            'q_keys': total_q,
            'rows': total_rows,
        },
        'claim_boundary': (
            'Actual output-key metadata extracted from persisted numerical master maps only. '
            'This is not a structural proof and not a numerical rank certificate.'
        ),
    }
    with gzip.open(out, 'wb', compresslevel=9) as f:
        pickle.dump(payload, f, protocol=5)
    payload['stats']['output_bytes'] = out.stat().st_size
    return payload


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', type=Path, required=True)
    ap.add_argument('--work', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--json', type=Path, required=True)
    ap.add_argument('--source-run-id')
    ap.add_argument('--source-artifact-id')
    ap.add_argument('--source-artifact-digest')
    ap.add_argument('--engine-git-blob-sha')
    a = ap.parse_args()
    x = extract_metadata(a.src, a.work, a.out, {
        'source_run_id': a.source_run_id,
        'source_artifact_id': a.source_artifact_id,
        'source_artifact_digest': a.source_artifact_digest,
        'engine_git_blob_sha': a.engine_git_blob_sha,
    })
    audit = {k: x[k] for k in (
        'schema_version','kind','irrep','shard','shards','active_blocks',
        'target_columns','assigned_blocks','assigned_columns','source_tar',
        'source_tar_sha256','source_run_id','source_artifact_id',
        'source_artifact_digest','engine_git_blob_sha','stats','claim_boundary')}
    a.json.write_text(json.dumps(audit, indent=2, sort_keys=True) + '\n')
    print(json.dumps(x['stats'], indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
