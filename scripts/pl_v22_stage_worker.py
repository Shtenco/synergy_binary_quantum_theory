#!/usr/bin/env python3
"""BQG V20: ORIGINAL PL16 Lorentzian three-stage SU(2) path.

Requires the full production repository checkout. No synthetic physics is
used by run(). Unit-test mocks are NOT physical amplitude measurements.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import traceback
from pathlib import Path

import numpy as np
from v21_checkpoint_guard import validate_physical_stage
from v21_stage_logic import predecessor_ok

JMAX2 = 7


def ordered_operations(targets, mode, aux):
    a, b, c = targets
    i, j, k = aux
    if mode == "forward":
        return [("CV", c, k, i), ("CK", b, j, k), ("CK", a, i, j)]
    if mode == "adjoint":
        return [("CK", a, j, i), ("CK", b, k, j), ("CV", c, i, k)]
    raise ValueError("INVALID_MODE")


def scalar_filter(s):
    """Original filter is applied after leg 2, not earlier."""
    return {q: z for q, z in s.items() if q[2] in (0, 2)}


def do_leg(state, operation, source_node, cv, ck):
    name, target, left, right = operation
    if name == "CV":
        out, leak = cv(state, source_node, target, left, right, JMAX2)
        diag = {"CV_complete_basis_leakage": float(leak)}
    else:
        out, diag = ck(state, source_node, target, left, right, JMAX2)
        diag = {str(k): float(v) for k, v in diag.items()
                if isinstance(v, (float, int, np.floating, np.integer))}
    if any(not (math.isfinite(z.real) and math.isfinite(z.imag))
           for z in out.values()):
        raise ValueError("NONFINITE_AMPLITUDE")
    if any(not math.isfinite(v) for v in diag.values()):
        raise ValueError("NONFINITE_LEAKAGE")
    return out, diag


def save_charged(path, state, n_edges=32, n_vertices=16):
    """Preserve ALL complex amplitudes, no tolerance pruning."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(state.items())
    spins = np.array([q[0] for q, z in rows], dtype=np.int16).reshape(-1, n_edges)
    Ks = np.array([q[1] for q, z in rows], dtype=np.int16).reshape(-1, n_vertices)
    charges = np.array([q[2:6] for q, z in rows], dtype=np.int16).reshape(-1, 4)
    amp = np.array([z for q, z in rows], dtype=np.complex128)
    np.savez_compressed(path, spins=spins, K_other=Ks, charges=charges, amp=amp)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_charged(path, digest, n_edges=32, n_vertices=16):
    path = Path(path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError("CHECKPOINT_DIGEST_MISMATCH")
    with np.load(path, allow_pickle=False) as z:
        if set(z.files) != {"spins", "K_other", "charges", "amp"}:
            raise ValueError("WRONG_COVARIANT_CHECKPOINT_KEYS")
        s, k, c, a = [z[key] for key in ("spins", "K_other", "charges", "amp")]
        if not (s.shape == (len(a), n_edges)
                and k.shape == (len(a), n_vertices)
                and c.shape == (len(a), 4)
                and a.ndim == 1 and a.dtype.kind == "c"):
            raise ValueError("BAD_CHECKPOINT_SHAPE")
        if not np.isfinite(a).all():
            raise ValueError("BAD_CHECKPOINT_AMPLITUDES")
        state = {}
        for spins, K, charges, z in zip(s, k, c, a):
            key = (tuple(map(int, spins)), tuple(map(int, K)), *map(int, charges))
            if key in state:
                raise ValueError("DUPLICATE_CHECKPOINT_KEY")
            state[key] = complex(z)
    return state


def validate_predecessor(prev, expected, stage):
    if prev.get("stage") != stage - 1 or not prev.get("passed"):
        raise ValueError("PREVIOUS_STAGE_NOT_PASSED")
    for key in ("source_commit", "source_node", "omit", "cycle", "mode",
                "order", "aux", "Jmax2", "volume_backend", "volume_backend_sha256"):
        if prev.get(key) != expected.get(key):
            raise ValueError("PREDECESSOR_MISMATCH_" + key)
    if not re.fullmatch("[0-9a-f]{64}", prev.get("npz_sha256", "")):
        raise ValueError("BAD_PREDECESSOR_DIGEST")
    if not re.fullmatch("[0-9a-f]{64}", prev.get("input_state_sha256", "")):
        raise ValueError("MISSING_INPUT_STATE_PROVENANCE")


def run(opts):
    """Only production SU(2) primitives, original order and scalar filter."""
    if not re.fullmatch("[0-9a-f]{40}", opts.commit):
        raise ValueError("NEED_FROZEN_SOURCE_COMMIT")
    if not opts.volume_backend:
        raise ValueError("NEED_VOLUME_BACKEND_IDENTITY")
    try:
        actual_head = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, timeout=10
        ).strip()
    except (OSError, subprocess.SubprocessError) as exc:
        raise ValueError("REPOSITORY_COMMIT_NOT_VERIFIABLE") from exc
    if actual_head != opts.commit:
        raise ValueError("CHECKOUT_SOURCE_COMMIT_MISMATCH")
    if opts.stage not in (1, 2, 3):
        raise ValueError("BAD_STAGE")
    if len(opts.aux) != 3 or any(x not in (0, 1) for x in opts.aux):
        raise ValueError("BAD_AUX")

    from pl_dual_complex import DualComplex, seed_16cell_boundary
    from pl_peter_weyl_euclidean import PLPeterWeylEuclidean
    from pl_lorentzian_physical_pair_worker_v4 import pair_indices
    from pl_covariant_backend import install_pl_graph
    from tetrahedral_volume_backend import install_tetrahedral_volume_backend
    import peter_weyl_lorentzian_logical_projection_gate as LP

    D = DualComplex(seed_16cell_boundary())
    G = PLPeterWeylEuclidean(D)
    p, q, _, _ = pair_indices(D, opts.source_node, opts.omit, opts.cycle)
    backend_file = Path(__file__).with_name("tetrahedral_volume_backend.py")
    backend_sha = hashlib.sha256(backend_file.read_bytes()).hexdigest()
    perm = p if opts.order == "p" else q
    targets = tuple(D.neighbor[(opts.source_node, t)] for t in perm)
    allops = ordered_operations(targets, opts.mode, opts.aux)

    meta = {
        "source_commit": opts.commit, "source_node": opts.source_node,
        "omit": opts.omit, "cycle": opts.cycle,
        "mode": opts.mode, "order": opts.order,
        "aux": opts.aux, "Jmax2": JMAX2,
        "volume_backend": opts.volume_backend,
        "volume_backend_sha256": backend_sha,
        "graph": "PL-16cell", "stage": opts.stage,
        "operator_signature": list(allops[opts.stage - 1]),
        "post_second_scalar_filter": opts.stage >= 2
    }
    if opts.stage == 1:
        seed = ((1,) * len(G.EDGES), (0,) * D.n_tets)
        meta["input_state_sha256"] = hashlib.sha256(repr(seed).encode("utf-8")).hexdigest()
    else:
        if not opts.previous_prefix:
            raise ValueError("MISSING_PREVIOUS_STAGE")
        prevpath = opts.previous_prefix.with_suffix(".json")
        previous = json.loads(prevpath.read_text("utf-8"))
        predecessor_ok(previous, meta, opts.stage, allops[0], allops[1])
        meta["input_state_sha256"] = previous["input_state_sha256"]
        meta["predecessor_npz_sha256"] = previous["npz_sha256"]
        meta["predecessor_metadata_sha256"] = hashlib.sha256(prevpath.read_bytes()).hexdigest()
        state = load_charged(
            opts.previous_prefix.with_suffix(".npz"),
            previous["npz_sha256"], len(G.EDGES), D.n_tets
        )

    with install_tetrahedral_volume_backend():
        with install_pl_graph(G):
            oldJ = LP.JMAX2
            LP.JMAX2 = JMAX2
            restore, caches = LP.install_sine_cached_stack()
            try:
                if opts.stage == 1:
                    state = LP.CV.gauss_to_covariant({seed: 1. + 0j}, opts.source_node)
                after, diag = do_leg(
                    state, allops[opts.stage - 1], opts.source_node,
                    LP.RAW.COMP.C_volume_component,
                    LP.RAW.KCOMP.C_K_component
                )
            finally:
                restore()
                LP.JMAX2 = oldJ

    if opts.stage == 2:
        after = scalar_filter(after)
    guard = validate_physical_stage(diag, after, opts.stage)
    digest = save_charged(
        opts.output_prefix.with_suffix(".npz"),
        after, len(G.EDGES), D.n_tets
    )
    meta.update({
        "passed": True, "status": "REAL_ORIGINAL_STAGE_ONLY",
        "npz_sha256": digest, "state_support": len(after),
        "state_norm": math.sqrt(sum(abs(z)**2 for z in after.values())),
        "diagnostics": diag, "guard": guard,
        "physical_full_pair_completed": False,
        "physical_HDA_proved": False
    })
    opts.output_prefix.with_suffix(".json").write_text(
        json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return meta


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--stage", type=int, choices=[1, 2, 3], required=True)
    p.add_argument("--commit", default=os.environ.get("GITHUB_SHA", ""))
    p.add_argument("--volume-backend", required=True)
    p.add_argument("--source-node", type=int, default=0)
    p.add_argument("--omit", type=int, choices=range(4), default=0)
    p.add_argument("--cycle", type=int, choices=range(3), default=0)
    p.add_argument("--mode", choices=["forward", "adjoint"], required=True)
    p.add_argument("--order", choices=["p", "q"], required=True)
    p.add_argument("--aux", type=int, nargs=3, required=True)
    p.add_argument("--previous-prefix", type=Path)
    p.add_argument("--output-prefix", type=Path, required=True)
    a = p.parse_args()
    try:
        print(json.dumps(run(a), indent=2))
        return 0
    except Exception as e:
        report = {
            "status": "STAGE_FAILED_NOT_PHYSICS", "passed": False,
            "error": str(e), "error_type": type(e).__name__,
            "traceback": traceback.format_exc()
        }
        a.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        a.output_prefix.with_suffix(".json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(report, indent=2))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
