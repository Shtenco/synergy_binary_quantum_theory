#!/usr/bin/env python3
"""Aggregate 32 raw Euclidean master-image columns into mu2 and Lanczos R1.

Consumes a current Euclidean 5x32 one-hit packet with its pruning perturbation
bound and all 32 raw-reference Y_i=M_E|b_i> master-image columns. Produces
mu0, mu1=V0^dag Y, mu2=Y^dag Y and R1=mu2-mu1^dag mu1. It can also serialize
the actual orthonormal first master-Lanczos block Q1 and edge matrix B1.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import k5_peter_weyl_safe_hda_column as PW
import bqg_spectral_history_graph_gate as SPECTRAL


def decode_state(rows):
    out = {}
    for row in rows:
        key = (tuple(int(x) for x in row['spins']), tuple(int(x) for x in row['K_labels']))
        amp = row['amp']
        out[key] = complex(float(amp[0]), float(amp[1]))
    return out


def inner(a, b):
    if len(a) > len(b):
        return np.conj(inner(b, a))
    return sum((np.conj(z) * b.get(k, 0.0j) for k, z in a.items()), 0.0j)


def gram(states):
    n = len(states)
    G = np.zeros((n, n), complex)
    for i in range(n):
        for j in range(i, n):
            z = inner(states[i], states[j])
            G[i, j] = z
            G[j, i] = np.conj(z)
    return G


def opnorm(a):
    return float(np.linalg.norm(a, 2)) if a.size else 0.0


def add_scaled(dst, src, scale):
    for k, z in src.items():
        dst[k] = dst.get(k, 0.0j) + scale * z
    for k in [k for k, z in dst.items() if z == 0.0j]:
        del dst[k]


def encode_state(state):
    return [
        {'spins': [int(x) for x in spins], 'K_labels': [int(x) for x in Ks],
         'amp': [float(complex(z).real), float(complex(z).imag)]}
        for (spins, Ks), z in sorted(state.items(), key=lambda kv: repr(kv[0]))
    ]


def load_master_images(root: Path):
    files = sorted(root.rglob('master_image_*.json'))
    rows = {}
    operator_sha = None
    for p in files:
        payload = json.loads(p.read_text(encoding='utf-8'))
        if payload.get('schema') != 'BQG_EUCLIDEAN_MASTER_IMAGE_COLUMN_V1':
            continue
        i = int(payload['input_index'])
        if i in rows:
            raise RuntimeError(f'duplicate master image input {i}')
        if not payload.get('passed', False):
            raise RuntimeError(f'master image input {i} did not pass integrity checks')
        sha = str(payload.get('operator_source_sha256', ''))
        if operator_sha is None:
            operator_sha = sha
        elif sha != operator_sha:
            raise RuntimeError('operator source SHA mismatch across master-image columns')
        rows[i] = (payload, decode_state(payload['state']))
    if set(rows) != set(range(32)):
        raise RuntimeError(f'expected all 32 master images, got {sorted(rows)}')
    return [rows[i][1] for i in range(32)], operator_sha


def load_euclidean_retained(packet_root: Path):
    manifest_path = packet_root / 'euclidean_packet_manifest.json'
    if not manifest_path.exists():
        found = list(packet_root.rglob('euclidean_packet_manifest.json'))
        if len(found) != 1:
            raise RuntimeError('cannot uniquely locate euclidean_packet_manifest.json')
        manifest_path = found[0]
    root = manifest_path.parent
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    if not manifest.get('passed', False):
        raise RuntimeError('Euclidean packet did not pass')
    cert = manifest.get('pruning_error_certificate')
    if not isinstance(cert, dict) or 'unpruned_minus_retained_M_EE_operator_norm_upper_bound' not in cert:
        raise RuntimeError('Euclidean packet lacks required pruning perturbation certificate')
    prune_bound = float(cert['unpruned_minus_retained_M_EE_operator_norm_upper_bound'])

    by_node = {v: [None] * 32 for v in range(5)}
    for row in manifest.get('columns', []):
        if row.get('family') != 'E':
            continue
        v = int(row['node'])
        i = int(row['input_index'])
        p = root / row['path']
        payload = json.loads(p.read_text(encoding='utf-8'))
        by_node[v][i] = decode_state(payload['complete_gauss_outgoing_column']['state'])
    if any(x is None for v in range(5) for x in by_node[v]):
        raise RuntimeError('Euclidean packet is missing retained one-hit columns')

    M = np.zeros((32, 32), complex)
    for v in range(5):
        M += gram(by_node[v])
    return manifest, 0.5 * (M + M.conj().T), prune_bound


def run(master_root: Path, e_packet_root: Path):
    Y, operator_sha = load_master_images(master_root)
    e_manifest, mu1_retained, prune_bound = load_euclidean_retained(e_packet_root)
    basis = PW.basis_full_jhalf()

    mu0 = np.eye(32, dtype=complex)
    mu1_raw = np.zeros((32, 32), complex)
    for i, state in enumerate(Y):
        for j, key in enumerate(basis):
            mu1_raw[j, i] = state.get(key, 0.0j)
    mu1_herm_error = opnorm(mu1_raw - mu1_raw.conj().T)
    mu1_raw = 0.5 * (mu1_raw + mu1_raw.conj().T)
    mu2_raw = gram(Y)
    mu2 = 0.5 * (mu2_raw + mu2_raw.conj().T)

    retained_consistency = opnorm(mu1_raw - mu1_retained)
    numerical_slack = 2.0e-8 * max(1.0, opnorm(mu1_retained))
    retained_consistency_pass = retained_consistency <= prune_bound + numerical_slack

    rr = mu2 - mu1_raw.conj().T @ mu1_raw
    R1 = 0.5 * (rr + rr.conj().T)
    r_ev = np.linalg.eigvalsh(R1)
    r_scale = max(1.0, float(np.max(np.abs(r_ev))))
    r_tol = 3.0e-9 * r_scale
    r_rank = int(np.sum(r_ev > r_tol))
    r_psd = bool(float(np.min(r_ev)) >= -r_tol)
    residual_opnorm = max(0.0, float(np.max(r_ev)))

    # Actual sparse residual columns R = M V0 - V0 A0.
    residual_states = []
    for i, yi in enumerate(Y):
        ri = dict(yi)
        for j, key in enumerate(basis):
            z = mu1_raw[j, i]
            if z != 0.0j:
                ri[key] = ri.get(key, 0.0j) - z
                if ri[key] == 0.0j:
                    del ri[key]
        residual_states.append(ri)
    direct_R1_raw = gram(residual_states)
    direct_R1 = 0.5 * (direct_R1_raw + direct_R1_raw.conj().T)
    residual_gram_identity_error = opnorm(direct_R1 - R1)

    # Rank-revealing first block factorization R=Q1 B1.
    rev, rU = np.linalg.eigh(R1)
    support = rev > r_tol
    q1_states = []
    if np.any(support):
        Ur = rU[:, support]
        lr = rev[support]
        for a in range(len(lr)):
            st = {}
            for i, ri in enumerate(residual_states):
                c = Ur[i, a] / np.sqrt(lr[a])
                if c != 0.0j:
                    add_scaled(st, ri, c)
            q1_states.append(st)
        B1 = np.sqrt(lr)[:, None] * Ur.conj().T
        q1G = gram(q1_states)
        q1_orth_error = opnorm(q1G - np.eye(len(q1_states)))
        recurrence_error = 0.0
        for i, yi in enumerate(Y):
            rec = {}
            for j, key in enumerate(basis):
                if mu1_raw[j, i] != 0.0j:
                    rec[key] = rec.get(key, 0.0j) + mu1_raw[j, i]
            for a, qa in enumerate(q1_states):
                add_scaled(rec, qa, B1[a, i])
            diff = dict(yi)
            add_scaled(diff, rec, -1.0)
            recurrence_error = max(recurrence_error, np.sqrt(max(float(inner(diff, diff).real), 0.0)))
    else:
        B1 = np.zeros((0, 32), complex)
        q1_orth_error = 0.0
        recurrence_error = max(np.sqrt(max(float(inner(r, r).real), 0.0)) for r in residual_states)

    m1_ev = np.linalg.eigvalsh(mu1_raw)
    m2_ev = np.linalg.eigvalsh(mu2)
    checks = {
        'all_32_raw_master_images_present': len(Y) == 32,
        'operator_source_sha_present': bool(operator_sha),
        'euclidean_pruning_certificate_present': np.isfinite(prune_bound),
        'raw_boundary_overlap_hermitian': mu1_herm_error <= numerical_slack,
        'raw_vs_retained_mu1_within_pruning_bound': retained_consistency_pass,
        'mu1_positive_semidefinite': float(np.min(m1_ev)) >= -3e-9 * max(1.0, float(np.max(np.abs(m1_ev)))),
        'mu2_positive_semidefinite': float(np.min(m2_ev)) >= -3e-9 * max(1.0, float(np.max(np.abs(m2_ev)))),
        'R1_positive_semidefinite': r_psd,
        'direct_sparse_residual_Gram_matches_moment_R1': residual_gram_identity_error <= max(numerical_slack, 1e-8),
        'Q1_is_orthonormal_on_selected_support': q1_orth_error <= 5e-8,
        'first_block_Lanczos_recurrence_reconstructs_master_images': recurrence_error <= 5e-8,
    }

    result = {
        'schema': 'BQG_EUCLIDEAN_PROJECTED_HEAT_MU2_V1',
        'status': 'actual Euclidean master moments mu0/mu1/mu2 and first spectral-history edge Gram',
        'operator_source_sha256': operator_sha,
        'euclidean_packet_sha256': e_manifest.get('packet_sha256'),
        'mu1_retained_pruning_operator_bound': prune_bound,
        'mu1_raw_vs_retained_operator_error': retained_consistency,
        'numerical_slack': numerical_slack,
        'mu1_raw_hermiticity_error_before_symmetrization': mu1_herm_error,
        'mu1': {
            'eigenvalue_min': float(np.min(m1_ev)),
            'eigenvalue_max': float(np.max(m1_ev)),
            'trace': float(np.trace(mu1_raw).real),
            'frobenius_norm': float(np.linalg.norm(mu1_raw)),
            'matrix': SPECTRAL.encode_matrix(mu1_raw),
        },
        'mu2': {
            'eigenvalue_min': float(np.min(m2_ev)),
            'eigenvalue_max': float(np.max(m2_ev)),
            'trace': float(np.trace(mu2).real),
            'frobenius_norm': float(np.linalg.norm(mu2)),
            'matrix': SPECTRAL.encode_matrix(mu2),
        },
        'first_master_lanczos_residual_gram': {
            'definition': 'R1=mu2-mu1^dagger mu1 = B1^dagger B1',
            'rank': r_rank,
            'rank_tolerance': r_tol,
            'eigenvalue_min': float(np.min(r_ev)),
            'eigenvalue_max': float(np.max(r_ev)),
            'trace': float(np.trace(R1).real),
            'frobenius_norm': float(np.linalg.norm(R1)),
            'matrix': SPECTRAL.encode_matrix(R1),
            'B1_matrix': SPECTRAL.encode_matrix(B1),
            'Q1_rank': len(q1_states),
            'Q1_orthonormality_error': q1_orth_error,
            'direct_sparse_residual_Gram_identity_error': residual_gram_identity_error,
            'Lanczos_recurrence_max_sparse_norm_error': recurrence_error,
        },
        'integrity_checks': checks,
        'passed': bool(all(checks.values())),
        'claim_boundary': (
            'Actual finite Euclidean normal-master moments and first master-Krylov edge only. '
            'No finite termination is certified here; Lorentzian and HDA/Dtarget production data remain required for physical BQG history.'
        ),
    }

    spectral_packet = {
        'schema': SPECTRAL.SCHEMA,
        'depth': 0,
        'seed_label': 'actual_q2_boundary_Euclidean_master_history',
        'moments': {
            '0': SPECTRAL.encode_matrix(mu0),
            '1': SPECTRAL.encode_matrix(mu1_raw),
            '2': SPECTRAL.encode_matrix(mu2),
        },
        'termination_certificate': {
            'mode': 'direct_block_residual',
            'residual_norm': residual_opnorm ** 0.5,
            'certified': False,
            'reason': 'R1 is measured but a propagated numerical termination bound for the raw two-hit operator has not yet been certified',
        },
        'physical_preconditions': {
            'domain_complete': False,
            'master_constraint_certified': False,
            'quantum_hda_or_explicit_dtarget_certified': False,
            'source_seed_complete_for_claim': False,
        },
        'tolerances': {
            'hermiticity': max(numerical_slack, 1e-10),
            'mu0_identity': 1e-12,
            'hankel_rank_rtol': 1e-10,
            'psd': 1e-8,
            'moment_reproduction': 1e-8,
            'termination_residual': 1e-10,
            'zero_eigenvalue': 1e-10,
        },
        'heat_sigma': [0.0, 0.01, 0.1, 1.0, 10.0],
        'zeta_s': [0.5, 1.0, 2.0],
        'provenance': {
            'operator_source_sha256': operator_sha,
            'euclidean_packet_sha256': e_manifest.get('packet_sha256'),
            'mu2_schema': result['schema'],
        },
    }
    return result, spectral_packet, q1_states


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--master-images', type=Path, required=True)
    ap.add_argument('--euclidean-packet', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--spectral-packet', type=Path, required=True)
    ap.add_argument('--q1-dir', type=Path)
    a = ap.parse_args()
    result, packet, q1_states = run(a.master_images, a.euclidean_packet)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    a.spectral_packet.parent.mkdir(parents=True, exist_ok=True)
    a.spectral_packet.write_text(json.dumps(packet, indent=2) + '\n', encoding='utf-8')
    if a.q1_dir is not None:
        a.q1_dir.mkdir(parents=True, exist_ok=True)
        meta = {
            'schema': 'BQG_EUCLIDEAN_MASTER_LANCZOS_Q1_V1',
            'rank': len(q1_states),
            'B1_matrix': result['first_master_lanczos_residual_gram']['B1_matrix'],
            'operator_source_sha256': result['operator_source_sha256'],
            'claim_boundary': 'First Euclidean master-Lanczos block only; no full history/projector claim.'
        }
        (a.q1_dir / 'q1_manifest.json').write_text(json.dumps(meta, indent=2) + '\n', encoding='utf-8')
        for i, st in enumerate(q1_states):
            payload = {
                'schema': 'BQG_EUCLIDEAN_MASTER_LANCZOS_Q1_COLUMN_V1',
                'q1_index': i,
                'support': len(st),
                'state': encode_state(st),
            }
            (a.q1_dir / f'q1_{i:02d}.json').write_text(json.dumps(payload, separators=(',', ':')) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k not in ('mu1', 'mu2', 'first_master_lanczos_residual_gram')}, indent=2))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
