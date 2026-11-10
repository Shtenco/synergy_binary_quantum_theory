#!/usr/bin/env python3
"""Fail-closed L1 amplitude package producer.

This first adapter ingests the existing Euclidean full-E depth-two shards only.
It does NOT claim to produce the joint E+S+R physical package.  Its purpose is
to establish a provenance-bound, order-stable amplitude package and regress it
against the historical collector before additional operators are connected.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Mapping

import numpy as np

SCHEMA = "L1_AMPLITUDE_PACKAGE"
SCHEMA_VERSION = 1
REPS = (0, 1, 5)
NODES = tuple(range(24))
EXPECTED_BACKEND = "LocalPLPeterWeylEuclidean exact active-cone"
SOURCE_BLOB_SHA = {
    "worker": "9b1845aa263bacf23a0f79d5140bb82ceb56c8a1",
    "collector": "3ba2514ecd553fc54a59d3b56c1c43ce24aae09f",
}
LEGACY_SUM_TOL = 1e-10
GLOBAL_DROP_TOL = 1e-10
REGRESSION_TOL = 3e-7
SOURCE_TOL = 1e-12


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def canonical_json_sha(obj) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    return hashlib.sha256(raw).hexdigest()

def git_blob_sha(path: Path) -> str:
    data = Path(path).read_bytes()
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def state_norm(state: Mapping[bytes, complex]) -> float:
    return math.sqrt(math.fsum(abs(v) ** 2 for v in state.values()))


def state_relerr(a, b) -> float:
    keys = set(a) | set(b)
    num = math.sqrt(math.fsum(abs(a.get(k, 0j) - b.get(k, 0j)) ** 2 for k in keys))
    den = max(state_norm(b), 1e-300)
    return num / den


def load_state(path: Path, n_edges: int, n_vertices: int):
    with np.load(path, allow_pickle=False) as z:
        required = {"spins", "Ks", "amp"}
        if set(z.files) != required:
            raise ValueError(f"{path}: NPZ members must be exactly {sorted(required)}, got {z.files}")
        spins, Ks, amp = z["spins"], z["Ks"], z["amp"]
    if spins.dtype != np.dtype("int16") or Ks.dtype != np.dtype("int16"):
        raise ValueError(f"{path}: spin and intertwiner arrays must be int16")
    if amp.dtype != np.dtype("complex128"):
        raise ValueError(f"{path}: amp must be complex128")
    if spins.ndim != 2 or spins.shape[1] != n_edges:
        raise ValueError(f"{path}: invalid spins shape {spins.shape}")
    if Ks.ndim != 2 or Ks.shape[1] != n_vertices:
        raise ValueError(f"{path}: invalid Ks shape {Ks.shape}")
    if amp.ndim != 1 or len(amp) != len(spins) or len(amp) != len(Ks):
        raise ValueError(f"{path}: row count mismatch")
    if not np.isfinite(amp.real).all() or not np.isfinite(amp.imag).all():
        raise ValueError(f"{path}: non-finite complex amplitude")
    if np.any(spins < 0) or np.any(Ks < 0):
        raise ValueError(f"{path}: negative spin/intertwiner label")
    out = {}
    for i, a in enumerate(amp):
        key = spins[i].tobytes() + Ks[i].tobytes()
        if key in out:
            raise ValueError(f"{path}: duplicate state key at row {i}")
        out[key] = complex(a)
    return out, spins.copy(), Ks.copy()


def stable_sum_states(states):
    """Order-independent-by-contract compensated sum; no intermediate pruning."""
    buckets = defaultdict(list)
    for state in states:
        for key, value in state.items():
            buckets[key].append(complex(value))
    out = {}
    discarded_sq = []
    total_sq = []
    for key in sorted(buckets):
        vals = buckets[key]
        # math.fsum is correctly rounded for each real component; deterministic
        # sorting also makes the serialized result reproducible.
        vals.sort(key=lambda z: (z.real, z.imag))
        z = complex(math.fsum(v.real for v in vals), math.fsum(v.imag for v in vals))
        total_sq.append(abs(z) ** 2)
        if abs(z) > GLOBAL_DROP_TOL:
            out[key] = z
        else:
            discarded_sq.append(abs(z) ** 2)
    discarded_norm = math.sqrt(math.fsum(discarded_sq))
    total_norm = math.sqrt(math.fsum(total_sq))
    discarded_relative_norm = discarded_norm / max(total_norm, 1e-300)
    return out, discarded_norm, discarded_relative_norm


def legacy_add(dst, src, tol=LEGACY_SUM_TOL):
    """Exact historical collector behavior, retained only as regression oracle."""
    for key, value in src.items():
        z = dst.get(key, 0j) + value
        if abs(z) > tol:
            dst[key] = z
        elif key in dst:
            del dst[key]


def inner(a, b):
    # Compensated complex inner product; retain the full complex result.
    if len(a) <= len(b):
        terms = (np.conj(v) * b.get(k, 0j) for k, v in a.items())
    else:
        terms = (np.conj(a.get(k, 0j)) * v for k, v in b.items())
    terms = list(terms)
    return complex(math.fsum(z.real for z in terms),
                   math.fsum(z.imag for z in terms))


def matrix(left, right):
    return np.asarray([[inner(left[i], right[j]) for j in range(3)] for i in range(3)], complex)


def orbit(M):
    n = max(float(np.linalg.norm(M)), 1e-30)
    diag = [M[i, i].real for i in range(3)]
    adj = [M[0, 1].real, M[1, 0].real, M[1, 2].real, M[2, 1].real]
    opp = [M[0, 2].real, M[2, 0].real]
    a, b, c = float(np.mean(diag)), float(np.mean(adj)), float(np.mean(opp))
    fit = np.array([[a, b, c], [b, a, b], [c, b, a]])
    return {
        "a_same": a, "b_adjacent": b, "c_opposite": c,
        "hermiticity_relative_defect": float(np.linalg.norm(M - M.conj().T) / n),
        "max_imaginary_entry": float(np.max(np.abs(M.imag))),
        "diagonal_spread": float(max(diag) - min(diag)),
        "adjacent_spread": float(max(adj) - min(adj)),
        "opposite_spread": float(max(opp) - min(opp)),
        "three_representative_orbit_residual": float(np.linalg.norm(M.real - fit) / max(np.linalg.norm(M.real), 1e-30)),
        "lambda_A1": a + 4*b + c,
        "lambda_E": a - 2*b + c,
        "lambda_T2": a - c,
        "Delta_ET": 2*(c-b),
    }


def metrics(U, V):
    K, A, B = matrix(U, U), matrix(U, V), matrix(V, V)
    ko, ao, bo = orbit(K), orbit(A), orbit(B)
    dyn = {}
    for ir in ("A1", "E", "T2"):
        kval, aval, bval = ko[f"lambda_{ir}"], ao[f"lambda_{ir}"], bo[f"lambda_{ir}"]
        if kval <= 1e-12:
            raise ValueError(f"nonpositive K_{ir}={kval}")
        h1, h2 = aval/kval, bval/kval
        dyn[ir] = {"K": kval, "A": aval, "B": bval,
                   "h1_normalized": h1, "h2_normalized": h2,
                   "Sigma2_depth2": h2-h1*h1}
    return {
        "K": K, "A": A, "B": B, "K_orbit": ko, "A_orbit": ao, "B_orbit": bo,
        "dynamic_irreps": dyn,
        "Delta_ET_h2": dyn["E"]["h2_normalized"]-dyn["T2"]["h2_normalized"],
        "Delta_ET_Sigma2": dyn["E"]["Sigma2_depth2"]-dyn["T2"]["Sigma2_depth2"],
    }


def pack_state(state, n_edges, n_vertices):
    rows = sorted(state.items(), key=lambda kv: kv[0])
    if not rows:
        return (np.zeros((0, n_edges), np.int16),
                np.zeros((0, n_vertices), np.int16),
                np.zeros((0,), np.complex128))
    spins = np.stack([np.frombuffer(k[:2*n_edges], dtype=np.int16) for k, _ in rows])
    Ks = np.stack([np.frombuffer(k[2*n_edges:], dtype=np.int16) for k, _ in rows])
    amp = np.asarray([v for _, v in rows], np.complex128)
    return spins, Ks, amp


def run(root: Path, out_dir: Path, source_commit: str):
    root, out_dir = Path(root), Path(out_dir)
    if (not source_commit or len(source_commit) != 40
            or any(ch not in "0123456789abcdef" for ch in source_commit.lower())):
        raise ValueError("--source-commit must be a full 40-character hexadecimal Git SHA")
    scripts_dir = Path(__file__).resolve().parent
    actual_worker_blob = git_blob_sha(scripts_dir/"collective_l1_metric_edge_depth2_shard_worker.py")
    actual_collector_blob = git_blob_sha(scripts_dir/"collective_l1_metric_edge_depth2_shard_collect.py")
    if actual_worker_blob != SOURCE_BLOB_SHA["worker"]:
        raise ValueError(f"worker source blob mismatch: {actual_worker_blob}")
    if actual_collector_blob != SOURCE_BLOB_SHA["collector"]:
        raise ValueError(f"collector source blob mismatch: {actual_collector_blob}")
    producer_blob = git_blob_sha(Path(__file__))

    U, V, V_LEGACY, U_SHARDS, V_SHARDS, shard_hashes, shard_metas = [], [], [], [], [], {}, {}
    source_checks, shard_stats, aggregation = {}, {}, {}
    for edge in REPS:
        rows, sources, shard_states, meta_rows = [], [], [], []
        for node in NODES:
            stem = f"edge_{edge}_node_{node}"
            mp, up, vp = root/f"{stem}.json", root/f"u_{edge}_node_{node}.npz", root/f"v_{edge}_node_{node}.npz"
            if not all(p.is_file() for p in (mp, up, vp)):
                raise FileNotFoundError(f"missing shard triplet for {stem}")
            m = json.loads(mp.read_text(encoding="utf-8"))
            if m.get("passed") is not True or m.get("science_status") != "L1_METRIC_EDGE_DEPTH2_FULL_E_SHARD":
                raise ValueError(f"failed or wrong-science shard {stem}")
            if m.get("backend") != EXPECTED_BACKEND:
                raise ValueError(f"unexpected backend in {stem}: {m.get('backend')}")
            required_checks = {
                "L1_nodes_384", "L1_dual_links_768", "parent_has_24_fine_nodes",
                "edge_has_four_barycentric_chambers", "first_edge_state_nonzero",
                "one_node_depth2_state_nonzero", "finite_amplitudes",
                "depth2_spin_wall_j_le_5_over_2", "depth2_seed_parity_restored",
                "node_index_valid",
            }
            shard_checks = m.get("checks", {})
            if not required_checks.issubset(shard_checks) or not all(shard_checks.values()):
                raise ValueError(f"incomplete or failed worker checks in {stem}")
            expected_edge = {0: [0, 1], 1: [0, 2], 5: [2, 3]}[edge]
            if m.get("edge") != expected_edge:
                raise ValueError(f"edge label mismatch in {stem}: {m.get('edge')}")
            if int(m.get("edge_index", -1)) != edge or int(m.get("node_index", -1)) != node:
                raise ValueError(f"metadata index mismatch in {stem}")
            if int(m.get("parent_coarse_tetra", -1)) != 0:
                raise ValueError(f"unexpected parent tetra in {stem}")
            if not math.isclose(float(m.get("Jmax", -1)), 2.5, rel_tol=0.0, abs_tol=1e-12):
                raise ValueError(f"unexpected Jmax in {stem}")
            u, us, uk = load_state(up, 768, 384)
            v, vs, vk = load_state(vp, 768, 384)
            rows.append(m); sources.append(u); shard_states.append(v); meta_rows.append(m)
            for p in (mp, up, vp):
                shard_hashes[str(p.name)] = sha256_file(p)
        u0 = sources[0]
        mismatch = max((state_relerr(s, u0) for s in sources[1:]), default=0.0)
        if mismatch >= SOURCE_TOL:
            raise ValueError(f"source state mismatch for edge {edge}: {mismatch}")
        U.append(u0)
        U_SHARDS.extend(sources)
        V_SHARDS.extend(shard_states)
        v_stable, discarded, discarded_rel = stable_sum_states(shard_states)
        v_reverse, _, _ = stable_sum_states(list(reversed(shard_states)))
        permutation_error = state_relerr(v_stable, v_reverse)
        if permutation_error != 0.0:
            raise ValueError(f"stable aggregation is order-dependent for edge {edge}: {permutation_error}")
        v_legacy = {}
        for state in shard_states:
            legacy_add(v_legacy, state)
        mismatch_agg = state_relerr(v_stable, v_legacy)
        if mismatch_agg > REGRESSION_TOL:
            raise ValueError(f"new aggregation does not reproduce legacy collector for edge {edge}: {mismatch_agg}")
        V.append(v_stable)
        V_LEGACY.append(v_legacy)
        source_checks[str(edge)] = {
            "support": len(u0), "norm": state_norm(u0),
            "max_relative_source_mismatch_across_24_shards": mismatch,
        }
        shard_stats[str(edge)] = {
            "shard_support_min": min(len(s) for s in shard_states),
            "shard_support_max": max(len(s) for s in shard_states),
            "reconstructed_block_support": len(v_stable),
            "reconstructed_block_norm": state_norm(v_stable),
            "sum_of_shard_runtimes_seconds": math.fsum(float(m["second_seconds"]) for m in rows),
        }
        aggregation[str(edge)] = {
            "legacy_collector_relative_difference": mismatch_agg,
            "discarded_global_l2_norm": discarded,
            "discarded_global_l2_relative_norm": discarded_rel,
            "stable_sum_order": "state-key sorted; per-key real/imag math.fsum",
            "intermediate_pruning": False,
            "reverse_order_relative_difference": permutation_error,
        }
        shard_metas[str(edge)] = rows

    result = metrics(U, V)
    legacy_result = metrics(U, V_LEGACY)
    moment_regression = {
        name: float(np.linalg.norm(result[name] - legacy_result[name]) / max(np.linalg.norm(legacy_result[name]), 1e-300))
        for name in ("K", "A", "B")
    }
    finite = all(np.isfinite(M.real).all() and np.isfinite(M.imag).all()
                 for M in (result["K"], result["A"], result["B"]))
    symmetry = all(result[k]["hermiticity_relative_defect"] < 3e-7
                   and result[k]["three_representative_orbit_residual"] < 3e-7
                   for k in ("K_orbit", "A_orbit", "B_orbit"))
    variance = min(x["Sigma2_depth2"] for x in result["dynamic_irreps"].values()) > -3e-6
    checks = {
        "all_72_shards_present_and_passed": True,
        "source_states_identical": True,
        "finite_amplitudes_and_moments": bool(finite),
        "legacy_collector_state_regression": all(x["legacy_collector_relative_difference"] <= REGRESSION_TOL for x in aggregation.values()),
        "legacy_collector_moment_regression": all(x <= REGRESSION_TOL for x in moment_regression.values()),
        "global_threshold_loss_within_frozen_tolerance": all(
            x["discarded_global_l2_relative_norm"] <= REGRESSION_TOL
            for x in aggregation.values()),
        "S4_three_orbit_consistency": bool(symmetry),
        "nonnegative_depth2_variance_with_frozen_tolerance": bool(variance),
    }
    # This is intentionally a precursor package, not the complete L1 package.
    package = {
        "schema": SCHEMA, "schema_version": SCHEMA_VERSION,
        "package_status": "PASS_E_ONLY_REGRESSION_PRECURSOR" if all(checks.values()) else "FAIL_E_ONLY_REGRESSION",
        "science_status": "L1_EUCLIDEAN_E_ONLY_PRECURSOR_NOT_FULL_PHYSICS",
        "complete_physical_package": False,
        "missing_required_operators": ["S", "R"],
        "missing_production_products": ["joint_E_S_R_Krylov_basis", "metric_calibration_on_same_W",
                                        "full_Lorentzian_Schur_blocks", "direct_HDA_generators_and_residual"],
        "level": "L1", "source_commit": source_commit,
        "source_commit_attestation": "caller-supplied; must equal git rev-parse HEAD at production time; not independently queried by this script",
        "operator_registry": {
            "E": {"status": "PRESENT_EUCLIDEAN_PRECURSOR", "definition": "existing full physical-sine Euclidean L1 depth-two shard worker",
                  "backend": EXPECTED_BACKEND, "source_blob_sha": SOURCE_BLOB_SHA["worker"]},
            "S": {"status": "MISSING"}, "R": {"status": "MISSING"},
        },
        "carrier": {"nodes": 384, "dual_links": 768, "parent_coarse_tetra": 0,
                    "parent_fine_nodes": 24, "edge_representatives": list(REPS),
                    "Jmax": 2.5, "basis_status": "NOT_YET_JOINT_E_S_R_W"},
        "provenance": {
            "source_worker_blob_sha": SOURCE_BLOB_SHA["worker"],
            "source_collector_blob_sha": SOURCE_BLOB_SHA["collector"],
            "producer_blob_sha": producer_blob,
            "shard_file_sha256": shard_hashes,
            "shard_manifest_sha256": canonical_json_sha(shard_hashes),
            "source_commit": source_commit,
            "producer_source_sha256": sha256_file(Path(__file__).resolve()),
            "worker_identity_limitation": "worker/collector source blobs are pinned by this schema, but legacy shard metadata does not self-attest its generating source blob",
        },
        "source_consistency": source_checks, "shard_diagnostics": shard_stats,
        "aggregation_diagnostics": aggregation,
        "legacy_moment_regression_relative_errors": moment_regression,
        "moments": {
            "K_matrix": [[[float(z.real), float(z.imag)] for z in row] for row in result["K"]],
            "A_matrix": [[[float(z.real), float(z.imag)] for z in row] for row in result["A"]],
            "B_matrix": [[[float(z.real), float(z.imag)] for z in row] for row in result["B"]],
            "K_orbit": result["K_orbit"], "A_orbit": result["A_orbit"], "B_orbit": result["B_orbit"],
            "dynamic_irreps": result["dynamic_irreps"],
            "Delta_ET_h2": result["Delta_ET_h2"], "Delta_ET_Sigma2": result["Delta_ET_Sigma2"],
        },
        "checks": checks,
        "scope_guard": "No S/R amplitudes, common E+S+R W, metric calibration, Lorentzian Schur reduction, or production HDA are claimed.",
    }
    arrays = {}
    for i, state in enumerate(U):
        s, k, a = pack_state(state, 768, 384)
        arrays[f"u_{REPS[i]}_spins"], arrays[f"u_{REPS[i]}_Ks"], arrays[f"u_{REPS[i]}_amp"] = s, k, a
    for i, state in enumerate(V):
        s, k, a = pack_state(state, 768, 384)
        arrays[f"v_{REPS[i]}_spins"], arrays[f"v_{REPS[i]}_Ks"], arrays[f"v_{REPS[i]}_amp"] = s, k, a
    # Preserve every input shard amplitude, not merely its hash or aggregate.
    for idx, state in enumerate(U_SHARDS):
        edge, node = REPS[idx // 24], idx % 24
        s, k, a = pack_state(state, 768, 384)
        arrays[f"u_shard_{edge}_{node}_spins"], arrays[f"u_shard_{edge}_{node}_Ks"], arrays[f"u_shard_{edge}_{node}_amp"] = s, k, a
    for idx, state in enumerate(V_SHARDS):
        edge, node = REPS[idx // 24], idx % 24
        s, k, a = pack_state(state, 768, 384)
        arrays[f"v_shard_{edge}_{node}_spins"], arrays[f"v_shard_{edge}_{node}_Ks"], arrays[f"v_shard_{edge}_{node}_amp"] = s, k, a
    out_dir.mkdir(parents=True, exist_ok=True)
    npz_path = out_dir/"L1_AMPLITUDE_PACKAGE.npz"
    json_path = out_dir/"L1_AMPLITUDE_PACKAGE.json"
    with tempfile.NamedTemporaryFile(dir=out_dir, suffix=".npz", delete=False) as tmp:
        tmp_name = tmp.name
    try:
        np.savez_compressed(tmp_name, **arrays)
        os.replace(tmp_name, npz_path)
    finally:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
    package["amplitude_archive_sha256"] = sha256_file(npz_path)
    package["package_sha256"] = canonical_json_sha(package)
    manifest_text = json.dumps(package, indent=2, sort_keys=True)+"\n"
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=out_dir,
                                     suffix=".json", delete=False) as tmp:
        tmp.write(manifest_text)
        manifest_tmp = tmp.name
    try:
        os.replace(manifest_tmp, json_path)
    finally:
        if os.path.exists(manifest_tmp):
            os.unlink(manifest_tmp)
    return package


def verify_package(manifest_path: Path):
    manifest_path = Path(manifest_path)
    package = json.loads(manifest_path.read_text(encoding="utf-8"))
    if package.get("schema") != SCHEMA or package.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("unknown L1 amplitude package schema/version")
    recorded = package.get("package_sha256")
    unsigned = dict(package)
    unsigned.pop("package_sha256", None)
    if not recorded or canonical_json_sha(unsigned) != recorded:
        raise ValueError("package manifest SHA256 mismatch")
    archive = manifest_path.with_name("L1_AMPLITUDE_PACKAGE.npz")
    if not archive.is_file() or sha256_file(archive) != package.get("amplitude_archive_sha256"):
        raise ValueError("amplitude archive missing or SHA256 mismatch")
    if package.get("complete_physical_package") is not False:
        raise ValueError("v1 E-only verifier refuses an unrecognized completeness claim")
    if package.get("science_status") != "L1_EUCLIDEAN_E_ONLY_PRECURSOR_NOT_FULL_PHYSICS":
        raise ValueError("unexpected science status")
    provenance = package.get("provenance", {})
    if provenance.get("source_worker_blob_sha") != SOURCE_BLOB_SHA["worker"]:
        raise ValueError("worker source provenance mismatch")
    if provenance.get("source_collector_blob_sha") != SOURCE_BLOB_SHA["collector"]:
        raise ValueError("collector source provenance mismatch")
    if len(provenance.get("producer_blob_sha", "")) != 40:
        raise ValueError("missing producer Git blob identity")
    shard_hashes = provenance.get("shard_file_sha256", {})
    expected_shard_files = {
        name
        for edge in REPS for node in NODES
        for name in (f"edge_{edge}_node_{node}.json",
                     f"u_{edge}_node_{node}.npz",
                     f"v_{edge}_node_{node}.npz")
    }
    if set(shard_hashes) != expected_shard_files or any(
            len(v) != 64 or any(ch not in "0123456789abcdef" for ch in v.lower())
            for v in shard_hashes.values()):
        raise ValueError("incomplete or malformed shard SHA256 manifest")
    if canonical_json_sha(shard_hashes) != provenance.get("shard_manifest_sha256"):
        raise ValueError("shard manifest digest mismatch")
    expected = set()
    for edge in REPS:
        for prefix in ("u", "v"):
            expected.update({f"{prefix}_{edge}_spins", f"{prefix}_{edge}_Ks", f"{prefix}_{edge}_amp"})
        for node in NODES:
            for prefix in ("u_shard", "v_shard"):
                expected.update({f"{prefix}_{edge}_{node}_spins", f"{prefix}_{edge}_{node}_Ks",
                                 f"{prefix}_{edge}_{node}_amp"})
    groups = [f"{p}_{edge}" for edge in REPS for p in ("u", "v")]
    groups += [f"{p}_shard_{edge}_{node}" for edge in REPS for node in NODES
               for p in ("u", "v")]
    with np.load(archive, allow_pickle=False) as z:
        if set(z.files) != expected:
            raise ValueError("amplitude archive member set mismatch")
        for group in groups:
            spins, Ks, amp = z[group+"_spins"], z[group+"_Ks"], z[group+"_amp"]
            if spins.dtype != np.dtype("int16") or spins.ndim != 2 or spins.shape[1] != 768:
                raise ValueError(f"invalid spin array: {group}")
            if Ks.dtype != np.dtype("int16") or Ks.ndim != 2 or Ks.shape[1] != 384:
                raise ValueError(f"invalid intertwiner array: {group}")
            if amp.dtype != np.dtype("complex128") or amp.ndim != 1:
                raise ValueError(f"invalid amplitude array: {group}")
            if len(spins) != len(Ks) or len(spins) != len(amp):
                raise ValueError(f"row count mismatch: {group}")
            if not np.isfinite(amp.real).all() or not np.isfinite(amp.imag).all():
                raise ValueError(f"non-finite amplitudes: {group}")
            seen = set()
            for i in range(len(amp)):
                key = spins[i].tobytes() + Ks[i].tobytes()
                if key in seen:
                    raise ValueError(f"duplicate state key: {group}")
                seen.add(key)
    return {"passed": True, "schema": SCHEMA, "package_sha256": recorded,
            "amplitude_archive_sha256": package["amplitude_archive_sha256"],
            "array_count": len(expected), "package_status": package.get("package_status")}


def self_test():
    # Demonstrate the exact failure mode of intermediate thresholding:
    # a small term is discarded before a later term makes the total material.
    key = b"k"
    a, b = {key: 0.6e-10}, {key: 0.6e-10}
    legacy = {}
    legacy_add(legacy, a); legacy_add(legacy, b)
    stable, _, _ = stable_sum_states([a, b])
    check1 = key not in legacy and abs(stable[key]-1.2e-10) < 1e-25
    # Stable aggregation is invariant under shard permutation.
    x = {b"x": 1.25+2.5j, b"y": -0.5j}
    y = {b"x": -0.25-1.5j, b"z": 3.0}
    s1, _, _ = stable_sum_states([x, y])
    s2, _, _ = stable_sum_states([y, x])
    check2 = state_relerr(s1, s2) == 0.0
    # Complex inner product must preserve conjugation and imaginary parts.
    ip = inner({b"z": 1.0+2.0j}, {b"z": 3.0-1.0j})
    check3 = abs(ip - complex(1.0-2.0j) * complex(3.0-1.0j)) < 1e-14
    checks = {"intermediate_threshold_counterexample_detected": bool(check1),
              "stable_aggregation_permutation_invariant": bool(check2),
              "complex_inner_product_preserves_phase": bool(check3)}
    return {"passed": all(checks.values()), "checks": checks}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path)
    ap.add_argument("--output-dir", type=Path, default=Path("verification_results/l1_amplitude_package"))
    ap.add_argument("--source-commit", default=os.environ.get("GITHUB_SHA", ""))
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--verify-manifest", type=Path)
    args = ap.parse_args()
    if args.verify_manifest:
        try:
            out = verify_package(args.verify_manifest)
        except Exception as exc:
            out = {"passed": False, "error_type": type(exc).__name__, "error": str(exc)}
        print(json.dumps(out, indent=2))
        return 0 if out.get("passed") else 1
    if args.self_test:
        out = self_test()
        print(json.dumps(out, indent=2))
        return 0 if out["passed"] else 1
    if args.root is None:
        ap.error("--root is required unless --self-test is used")
    try:
        out = run(args.root, args.output_dir, args.source_commit)
    except Exception as exc:
        fail = {"schema": SCHEMA, "schema_version": SCHEMA_VERSION,
                "package_status": "FAIL_E_ONLY_REGRESSION",
                "science_status": "L1_AMPLITUDE_PACKAGE_INPUT_FAIL",
                "error_type": type(exc).__name__, "error": str(exc)}
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir/"L1_AMPLITUDE_PACKAGE.json").write_text(json.dumps(fail, indent=2)+"\n")
        print(json.dumps(fail, indent=2))
        return 1
    print(json.dumps({k: v for k, v in out.items() if k not in ("provenance", "moments", "shard_diagnostics")}, indent=2))
    return 0 if out["package_status"] == "PASS_E_ONLY_REGRESSION_PRECURSOR" else 1


if __name__ == "__main__":
    raise SystemExit(main())
