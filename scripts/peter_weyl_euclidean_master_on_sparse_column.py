#!/usr/bin/env python3
"""Apply the fixed-cutoff Euclidean normal master to one serialized sparse state.

For the explicitly Hermitian frozen Euclidean node operator,

    M_E = sum_v H_v^dag H_v = sum_v H_v^2.

This producer is intended for direct block-Lanczos continuation after Q1 is
available. It does not call compose_on_sparse and performs no additional
post-second-hit tolerance pruning beyond the thresholds internal to the frozen
Peter-Weyl reference implementation.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k5_peter_weyl_safe_hda_column as PW


def decode(rows):
    out = {}
    for row in rows:
        key = (tuple(int(x) for x in row['spins']), tuple(int(x) for x in row['K_labels']))
        a = row['amp']
        out[key] = complex(float(a[0]), float(a[1]))
    return out


def encode(state):
    return [
        {'spins': [int(x) for x in spins], 'K_labels': [int(x) for x in Ks],
         'amp': [float(complex(z).real), float(complex(z).imag)]}
        for (spins, Ks), z in sorted(state.items(), key=lambda kv: repr(kv[0]))
    ]


def add_exact(dst, src):
    for k, z in src.items():
        dst[k] = dst.get(k, 0.0j) + complex(z)
    for k in [k for k, z in dst.items() if z == 0.0j]:
        del dst[k]


def norm(state):
    return math.sqrt(sum(abs(z) ** 2 for z in state.values()))


def max_spin(state):
    return max((max(k[0]) / 2.0 for k in state), default=0.0)


def run(input_path: Path, jmax2: int):
    payload = json.loads(input_path.read_text(encoding='utf-8'))
    if 'state' not in payload:
        raise ValueError('input payload has no sparse state')
    state = decode(payload['state'])
    if not state:
        raise ValueError('empty input state')

    total = {}
    nodes = []
    finite = True
    PW.T_cached.cache_clear()
    for v in range(5):
        first = PW.apply_H_cached_state(state, v, jmax2)
        second = PW.apply_H_cached_state(first, v, jmax2)
        add_exact(total, second)
        finite = finite and all(np.isfinite(complex(z).real) and np.isfinite(complex(z).imag) for z in first.values())
        finite = finite and all(np.isfinite(complex(z).real) and np.isfinite(complex(z).imag) for z in second.values())
        nodes.append({
            'node': v,
            'first_support': len(first),
            'first_norm': norm(first),
            'second_support': len(second),
            'second_norm': norm(second),
            'second_max_spin': max_spin(second),
        })
        PW.T_cached.cache_clear()
        gc.collect()

    rows = encode(total)
    source_sha = hashlib.sha256(Path(PW.__file__).read_bytes()).hexdigest()
    checks = {
        'five_node_contributions': len(nodes) == 5,
        'finite_amplitudes': bool(finite),
        'fixed_cutoff_respected': max_spin(total) <= jmax2 / 2.0 + 1e-12,
        'no_compose_on_sparse': True,
        'no_additional_post_second_hit_prune': True,
    }
    return {
        'schema': 'BQG_EUCLIDEAN_MASTER_ON_SPARSE_COLUMN_V1',
        'passed': bool(all(checks.values())),
        'source_schema': payload.get('schema'),
        'source_label': payload.get('q1_index', payload.get('input_index')),
        'Jmax': jmax2 / 2.0,
        'operator_source_sha256': source_sha,
        'input_support': len(state),
        'input_norm': norm(state),
        'output_support': len(total),
        'output_norm': norm(total),
        'output_max_spin': max_spin(total),
        'node_contributions': nodes,
        'checks': checks,
        'state': rows,
        'claim_boundary': 'One fixed-regulator Euclidean master action on a supplied sparse Lanczos state; no spectral termination or physical projector claim.',
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--jmax2', type=int, default=5)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    out = run(a.input, a.jmax2)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, separators=(',', ':')) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in out.items() if k != 'state'}, indent=2))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
