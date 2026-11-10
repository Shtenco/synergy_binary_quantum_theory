#!/usr/bin/env python3
"""Validate an L1_AMPLITUDE_PACKAGE archive and its immutable manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

SCHEMA = "L1_AMPLITUDE_PACKAGE"
VERSION = 1
EDGES = (0, 1, 5)
NODES = range(24)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def git_blob_sha(path: Path) -> str:
    data = Path(path).read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def canonical_sha(obj) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    return hashlib.sha256(raw).hexdigest()


def expected_names():
    names = set()
    for e in EDGES:
        for prefix in ("u", "v"):
            names.update({f"{prefix}_{e}_spins", f"{prefix}_{e}_Ks", f"{prefix}_{e}_amp"})
        for n in NODES:
            for prefix in ("u_shard", "v_shard"):
                names.update({f"{prefix}_{e}_{n}_spins",
                              f"{prefix}_{e}_{n}_Ks",
                              f"{prefix}_{e}_{n}_amp"})
    return names


def _state_from_npz(path):
    with np.load(path, allow_pickle=False) as data:
        if set(data.files) != {"spins", "Ks", "amp"}:
            raise ValueError(f"source shard {path.name}: unexpected NPZ members")
        spins, Ks, amp = data["spins"], data["Ks"], data["amp"]
        if (spins.dtype != np.dtype("int16") or spins.ndim != 2 or spins.shape[1] != 768
                or Ks.dtype != np.dtype("int16") or Ks.ndim != 2 or Ks.shape[1] != 384
                or amp.dtype != np.dtype("complex128") or amp.ndim != 1
                or len(spins) != len(Ks) or len(spins) != len(amp)):
            raise ValueError(f"source shard {path.name}: invalid array dtype/shape")
        if (not np.isfinite(amp.real).all() or not np.isfinite(amp.imag).all()
                or np.any(spins < 0) or np.any(Ks < 0)):
            raise ValueError(f"source shard {path.name}: invalid labels or amplitudes")
        state = {}
        for i, value in enumerate(amp):
            key = spins[i].tobytes() + Ks[i].tobytes()
            if key in state:
                raise ValueError(f"source shard {path.name}: duplicate state key")
            state[key] = complex(value)
        return state


def _validate_shard_metadata(path, edge, node):
    meta = json.loads(path.read_text(encoding="utf-8"))
    expected_edges = {0: [0, 1], 1: [0, 2], 5: [2, 3]}
    required = {
        "L1_nodes_384", "L1_dual_links_768", "parent_has_24_fine_nodes",
        "edge_has_four_barycentric_chambers", "first_edge_state_nonzero",
        "one_node_depth2_state_nonzero", "finite_amplitudes",
        "depth2_spin_wall_j_le_5_over_2", "depth2_seed_parity_restored",
        "node_index_valid",
    }
    if meta.get("passed") is not True or meta.get("science_status") != "L1_METRIC_EDGE_DEPTH2_FULL_E_SHARD":
        raise ValueError(f"source shard metadata failed or wrong science status: {path.name}")
    if meta.get("backend") != "LocalPLPeterWeylEuclidean exact active-cone":
        raise ValueError(f"source shard metadata backend mismatch: {path.name}")
    if type(meta.get("edge_index")) is not int or meta["edge_index"] != edge:
        raise ValueError(f"source shard edge index mismatch: {path.name}")
    if type(meta.get("node_index")) is not int or meta["node_index"] != node:
        raise ValueError(f"source shard node index mismatch: {path.name}")
    if type(meta.get("parent_coarse_tetra")) is not int or meta["parent_coarse_tetra"] != 0:
        raise ValueError(f"source shard parent tetra mismatch: {path.name}")
    if meta.get("edge") != expected_edges[edge] or meta.get("Jmax") != 2.5:
        raise ValueError(f"source shard edge/cutoff mismatch: {path.name}")
    if meta.get("definition") != "u_e=(1/2) sum_{4 chambers->e} H_c|Omega>; v_e,w=H_w u_e; collector sums w=0..23":
        raise ValueError(f"source shard operator definition mismatch: {path.name}")
    checks = meta.get("checks", {})
    if not required.issubset(checks) or any(value is not True for value in checks.values()):
        raise ValueError(f"source shard checks missing or failed: {path.name}")
    return meta


def _state_from_arrays(z, stem):
    spins, Ks, amp = z[stem+"_spins"], z[stem+"_Ks"], z[stem+"_amp"]
    return {
        spins[i].tobytes() + Ks[i].tobytes(): complex(amp[i])
        for i in range(len(amp))
    }


def _state_norm(state):
    return math.sqrt(math.fsum(abs(value) ** 2 for value in state.values()))


def _state_relerr(actual, expected):
    keys = set(actual) | set(expected)
    numerator = math.sqrt(math.fsum(
        abs(actual.get(key, 0j) - expected.get(key, 0j)) ** 2 for key in keys
    ))
    return numerator / max(_state_norm(expected), 1e-300)


def _stable_sum_states(states):
    contributions = {}
    for state in states:
        for key, value in state.items():
            contributions.setdefault(key, []).append(complex(value))
    result = {}
    for key in sorted(contributions):
        values = sorted(contributions[key], key=lambda value: (value.real, value.imag))
        value = complex(
            math.fsum(item.real for item in values),
            math.fsum(item.imag for item in values),
        )
        if abs(value) > 1e-10:
            result[key] = value
    return result


def validate(manifest_path: Path, archive_path: Path, producer_source_path: Path, shard_root: Path):
    manifest_path, archive_path = Path(manifest_path), Path(archive_path)
    p = json.loads(manifest_path.read_text(encoding="utf-8"))
    if p.get("schema") != SCHEMA or p.get("schema_version") != VERSION:
        raise ValueError("unknown package schema/version")
    if p.get("complete_physical_package") is not False:
        raise ValueError("v1 E-only package cannot claim complete physical status")
    if p.get("package_status") != "PASS_E_ONLY_REGRESSION_PRECURSOR":
        raise ValueError("package is not a passed E-only regression precursor")
    if p.get("science_status") != "L1_EUCLIDEAN_E_ONLY_PRECURSOR_NOT_FULL_PHYSICS":
        raise ValueError("unexpected science status")
    if p.get("missing_required_operators") != ["S", "R"]:
        raise ValueError("S/R absence must remain explicit")
    if (not isinstance(p.get("source_commit"), str) or len(p["source_commit"]) != 40
            or any(ch not in "0123456789abcdef" for ch in p["source_commit"].lower())):
        raise ValueError("missing full hexadecimal source commit")
    if not p.get("source_commit_attestation", "").startswith("caller-supplied;"):
        raise ValueError("source commit attestation limitation must remain explicit")
    claimed_package_sha = p.get("package_sha256")
    body = dict(p)
    body.pop("package_sha256", None)
    if claimed_package_sha != canonical_sha(body):
        raise ValueError("manifest package_sha256 mismatch")
    if p.get("amplitude_archive_sha256") != sha256(archive_path):
        raise ValueError("amplitude archive SHA256 mismatch")
    provenance = p.get("provenance", {})
    hashes = provenance.get("shard_file_sha256", {})
    expected_shard_files = {
        name
        for e in EDGES for n in NODES
        for name in (f"edge_{e}_node_{n}.json",
                     f"u_{e}_node_{n}.npz",
                     f"v_{e}_node_{n}.npz")
    }
    if set(hashes) != expected_shard_files or any(
            len(v) != 64 or any(ch not in "0123456789abcdef" for ch in v.lower())
            for v in hashes.values()
    ):
        raise ValueError("expected exact hashes for all 72 shard JSON/NPZ triplets")
    if provenance.get("shard_manifest_sha256") != canonical_sha(hashes):
        raise ValueError("shard manifest digest mismatch")
    producer_source_path, shard_root = Path(producer_source_path), Path(shard_root)
    if provenance.get("producer_source_sha256") != sha256(producer_source_path):
        raise ValueError("producer source SHA256 mismatch")
    if provenance.get("producer_blob_sha") != git_blob_sha(producer_source_path):
        raise ValueError("producer Git blob SHA mismatch")
    for name, expected_hash in hashes.items():
        candidate = shard_root / name
        if not candidate.is_file() or sha256(candidate) != expected_hash:
            raise ValueError(f"missing or modified source shard file: {name}")
    if "worker_identity_limitation" not in provenance:
        raise ValueError("worker source self-attestation limitation must remain explicit")
    registry = p.get("operator_registry", {})
    if registry.get("E", {}).get("source_blob_sha") != "9b1845aa263bacf23a0f79d5140bb82ceb56c8a1":
        raise ValueError("unrecognized E worker source blob")
    if provenance.get("source_collector_blob_sha") != "3ba2514ecd553fc54a59d3b56c1c43ce24aae09f":
        raise ValueError("unrecognized historical collector source blob")
    if registry.get("E", {}).get("backend") != "LocalPLPeterWeylEuclidean exact active-cone":
        raise ValueError("unexpected E backend")
    if p.get("carrier") != {
        "nodes": 384, "dual_links": 768, "parent_coarse_tetra": 0,
        "parent_fine_nodes": 24, "edge_representatives": [0, 1, 5],
        "Jmax": 2.5, "basis_status": "NOT_YET_JOINT_E_S_R_W"
    }:
        raise ValueError("unexpected or mutated carrier contract")
    required_checks = {
        "all_72_shards_present_and_passed",
        "source_states_identical",
        "finite_amplitudes_and_moments",
        "legacy_collector_state_regression",
        "legacy_collector_moment_regression",
        "global_threshold_loss_within_frozen_tolerance",
        "S4_three_orbit_consistency",
        "nonnegative_depth2_variance_with_frozen_tolerance",
    }
    checks = p.get("checks", {})
    if set(checks) != required_checks or not all(value is True for value in checks.values()):
        raise ValueError("missing, unexpected, or failed producer checks")
    if provenance.get("source_worker_blob_sha") != "9b1845aa263bacf23a0f79d5140bb82ceb56c8a1":
        raise ValueError("unrecognized worker source provenance")
    if p.get("operator_registry", {}).get("S", {}).get("status") != "MISSING" or p.get("operator_registry", {}).get("R", {}).get("status") != "MISSING":
        raise ValueError("S/R must be explicitly missing in v1")

    expected = expected_names()
    with np.load(archive_path, allow_pickle=False) as z:
        actual = set(z.files)
        if actual != expected:
            raise ValueError(f"archive member set mismatch: missing={len(expected-actual)} extra={len(actual-expected)}")
        checked_states = 0
        aggregation_reconstruction = {}
        source_reconstruction = {}
        for e in EDGES:
            for prefix in ("u", "v"):
                checked_states += _check_state(z, f"{prefix}_{e}", 768, 384)
            for n in NODES:
                for prefix in ("u_shard", "v_shard"):
                    checked_states += _check_state(z, f"{prefix}_{e}_{n}", 768, 384)
            u_shards = []
            v_shards = []
            for n in NODES:
                stem = f"edge_{e}_node_{n}"
                _validate_shard_metadata(shard_root / f"{stem}.json", e, n)
                for prefix, archive_prefix, target in (
                    ("u", "u_shard", u_shards), ("v", "v_shard", v_shards)
                ):
                    source_state = _state_from_npz(shard_root / f"{prefix}_{e}_node_{n}.npz")
                    archived_state = _state_from_arrays(z, f"{archive_prefix}_{e}_{n}")
                    if source_state != archived_state:
                        raise ValueError(
                            f"source NPZ and archived raw amplitudes differ: {prefix}_{e}_node_{n}"
                        )
                    target.append(source_state)
            u_reference = u_shards[0]
            u_error = max((_state_relerr(state, u_reference) for state in u_shards[1:]), default=0.0)
            archived_u_error = _state_relerr(_state_from_arrays(z, f"u_{e}"), u_reference)
            if u_error > 1e-12 or archived_u_error > 1e-12:
                raise ValueError(f"source amplitude reconstruction mismatch for edge {e}: shards={u_error}, aggregate={archived_u_error}")
            v_reconstructed = _stable_sum_states(v_shards)
            v_error = _state_relerr(_state_from_arrays(z, f"v_{e}"), v_reconstructed)
            if v_error > 1e-14:
                raise ValueError(f"depth-two aggregate does not equal raw shard sum for edge {e}: {v_error}")
            source_reconstruction[str(e)] = {
                "max_shard_relative_error": u_error,
                "archived_source_relative_error": archived_u_error,
            }
            aggregation_reconstruction[str(e)] = {
                "archived_response_relative_error": v_error,
                "reconstructed_support": len(v_reconstructed),
            }
    return {
        "passed": True,
        "science_status": "PACKAGE_INTEGRITY_VALID_E_ONLY_PRECURSOR",
        "checked_sparse_state_groups": checked_states,
        "verified_source_shard_files": len(hashes),
        "independently_reconstructed_source_states": source_reconstruction,
        "independently_reconstructed_depth2_aggregates": aggregation_reconstruction,
        "archive_sha256": p["amplitude_archive_sha256"],
        "package_sha256": claimed_package_sha,
        "physical_completeness": False,
    }


def _check_state(z, stem, n_edges, n_vertices):
    spins, Ks, amp = z[f"{stem}_spins"], z[f"{stem}_Ks"], z[f"{stem}_amp"]
    if spins.dtype != np.dtype("int16") or spins.ndim != 2 or spins.shape[1] != n_edges:
        raise ValueError(f"{stem}: invalid spins dtype/shape")
    if Ks.dtype != np.dtype("int16") or Ks.ndim != 2 or Ks.shape[1] != n_vertices:
        raise ValueError(f"{stem}: invalid Ks dtype/shape")
    if amp.dtype != np.dtype("complex128") or amp.ndim != 1 or len(amp) != len(spins) or len(amp) != len(Ks):
        raise ValueError(f"{stem}: invalid amplitude dtype/shape")
    if not np.isfinite(amp.real).all() or not np.isfinite(amp.imag).all():
        raise ValueError(f"{stem}: non-finite amplitude")
    if np.any(spins < 0) or np.any(Ks < 0):
        raise ValueError(f"{stem}: negative spin/intertwiner label")
    keys = [spins[i].tobytes()+Ks[i].tobytes() for i in range(len(amp))]
    if len(set(keys)) != len(keys):
        raise ValueError(f"{stem}: duplicate state key")
    if keys != sorted(keys):
        raise ValueError(f"{stem}: sparse rows are not in canonical key order")
    return 1


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--producer-source", type=Path, required=True)
    ap.add_argument("--shard-root", type=Path, required=True)
    a = ap.parse_args()
    try:
        result = validate(a.manifest, a.archive, a.producer_source, a.shard_root)
    except Exception as exc:
        print(json.dumps({"passed": False, "science_status": "PACKAGE_INTEGRITY_FAIL",
                          "error_type": type(exc).__name__, "error": str(exc)}, indent=2))
        return 1
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
