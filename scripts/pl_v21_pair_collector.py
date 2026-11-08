#!/usr/bin/env python3
"""Collect genuine complete 32-path V4 and compare with four independent V2."""
from __future__ import annotations
import argparse, hashlib, json, math, re, subprocess
from itertools import product
from pathlib import Path
import numpy as np
from v21_stage_logic import operations, predecessor_ok
from v21_checkpoint_guard import BOUND, validate_physical_stage
from pl_v20_stage_worker import load_charged


def read_stage(prefix, attestation):
    prefix = Path(prefix)
    meta_path = prefix.with_suffix(".json")
    npz_path = prefix.with_suffix(".npz")
    if not isinstance(attestation, dict):
        raise ValueError("MISSING_EXTERNAL_STAGE_ATTESTATION")
    metadata_pin = attestation.get("metadata_sha256")
    npz_pin = attestation.get("npz_sha256")
    if not isinstance(metadata_pin, str) or not isinstance(npz_pin, str):
        raise ValueError("INCOMPLETE_EXTERNAL_STAGE_ATTESTATION")
    if hashlib.sha256(meta_path.read_bytes()).hexdigest() != metadata_pin:
        raise ValueError("STAGE_METADATA_ATTESTATION_MISMATCH")
    if hashlib.sha256(npz_path.read_bytes()).hexdigest() != npz_pin:
        raise ValueError("STAGE_NPZ_ATTESTATION_MISMATCH")
    meta = json.loads(meta_path.read_text("utf-8"))
    if meta.get("status") != "REAL_ORIGINAL_STAGE_ONLY" or meta.get("passed") is not True:
        raise ValueError("MISSING_SUCCESSFUL_ORIGINAL_STAGE")
    if not isinstance(meta.get("diagnostics"), dict):
        raise ValueError("MISSING_STAGE_DIAGNOSTICS")
    signature = meta.get("operator_signature")
    state = load_charged(npz_path, meta.get("npz_sha256", ""))
    output_sha = meta.get("output_state_sha256")
    if output_sha != meta.get("npz_sha256"):
        raise ValueError("OUTPUT_CHECKPOINT_DIGEST_METADATA_MISMATCH")
    if len(state) != meta.get("state_support"):
        raise ValueError("STATE_SUPPORT_MISMATCH")
    recomputed = validate_physical_stage(
        meta["diagnostics"], state, meta.get("stage"), signature, BOUND)
    if meta.get("guard") != recomputed:
        raise ValueError("STORED_GUARD_RECOMPUTATION_MISMATCH")
    return meta, state


def certify_triple(item, shared, p_targets, q_targets):
    attestations = item.get("stage_attestations")
    if not isinstance(attestations, dict):
        raise ValueError("MISSING_TRIPLE_STAGE_ATTESTATIONS")
    m1, _ = read_stage(item["stage1"], attestations.get("stage1"))
    m2, _ = read_stage(item["stage2"], attestations.get("stage2"))
    m3, s3 = read_stage(item["stage3"], attestations.get("stage3"))
    mode, order, aux = m3["mode"], m3["order"], m3["aux"]
    if (mode not in ("forward", "adjoint") or order not in ("p", "q")
            or not isinstance(aux, list) or len(aux) != 3
            or any(x not in (0, 1) for x in aux)):
        raise ValueError("INVALID_PATH_INDEX")
    expected_ops = operations(p_targets if order == "p" else q_targets, mode, aux)
    for meta, stage, op in ((m1, 1, expected_ops[0]),
                            (m2, 2, expected_ops[1]),
                            (m3, 3, expected_ops[2])):
        if meta.get("stage") != stage or tuple(meta.get("operator_signature", ())) != tuple(op):
            raise ValueError("MISSING_STAGE_OR_WRONG_OPERATOR_SIGNATURE")
        for key in ("source_commit", "source_node", "omit", "cycle", "Jmax2",
                    "volume_backend", "volume_backend_sha256", "seed_state_sha256"):
            if meta.get(key) != shared[key]:
                raise ValueError("MIXED_PHYSICAL_SOURCE_OR_BACKEND_" + key)
    if m1.get("input_state_sha256") != m1.get("seed_state_sha256"):
        raise ValueError("STAGE1_INPUT_IS_NOT_DECLARED_SEED")
    predecessor_ok(m1, m2, 2, expected_ops[0], expected_ops[1])
    predecessor_ok(m2, m3, 3, expected_ops[0], expected_ops[1])
    return (mode, order, tuple(aux)), s3


def read_v2(spec, *, expected_source=None, expected_mode=None, expected_index=None,
            expected_source_commit=None, expected_backend_sha=None):
    file = Path(spec["npz"])
    if hashlib.sha256(file.read_bytes()).hexdigest() != spec["sha256"]:
        raise ValueError("INDEPENDENT_V2_SHA256_MISMATCH")
    if expected_source is not None:
        meta_path = Path(spec["metadata_json"])
        if hashlib.sha256(meta_path.read_bytes()).hexdigest() != spec["metadata_sha256"]:
            raise ValueError("V2_METADATA_SHA256_MISMATCH")
        meta = json.loads(meta_path.read_text("utf-8"))
        if (meta.get("passed") is not True or meta.get("mode") != expected_mode
                or meta.get("operator_version") != "tetrahedral-charged-volume-v2"):
            raise ValueError("V2_PRIMITIVE_FAILED_OR_WRONG_MODE_OR_VERSION")
        for field in ("index", "source_node"):
            value = meta.get(field)
            if not isinstance(value, int) or isinstance(value, bool):
                raise ValueError("V2_PRIMITIVE_BAD_INTEGER_" + field.upper())
        if meta.get("index") != expected_index or meta.get("source_node") != expected_source:
            raise ValueError("V2_PRIMITIVE_WRONG_INPUT_OR_ORIENTATION")
        if (expected_source_commit is not None
                and meta.get("source_commit") != expected_source_commit):
            raise ValueError("V2_SOURCE_COMMIT_MISMATCH")
        if (expected_backend_sha is not None
                and meta.get("volume_backend_sha256") != expected_backend_sha):
            raise ValueError("V2_VOLUME_BACKEND_MISMATCH")
        if meta.get("volume_backend") != "tetrahedral-volume-v1":
            raise ValueError("V2_WRONG_VOLUME_BACKEND")
        if float(meta.get("Jmax", -1)) != 3.5:
            raise ValueError("V2_PRIMITIVE_WRONG_CUTOFF")
        checks = meta.get("checks")
        if (not isinstance(checks, dict) or not checks
                or not all(value is True for value in checks.values())):
            raise ValueError("V2_PRIMITIVE_FAILED_ACCEPTANCE")
    with np.load(file, allow_pickle=False) as data:
        if set(data.files) != {"spins", "Ks", "amp"}:
            raise ValueError("INVALID_V2_FORMAT")
        spins, Ks, amp = data["spins"], data["Ks"], data["amp"]
        if (spins.shape != (len(amp), 32) or Ks.shape != (len(amp), 16)
                or spins.dtype.kind not in "iu" or Ks.dtype.kind not in "iu"
                or amp.ndim != 1 or amp.dtype.kind != "c" or not np.isfinite(amp).all()):
            raise ValueError("V2_BAD_SHAPE_DTYPE_OR_NONFINITE")
        result = {}
        for sp, ks, z in zip(spins, Ks, amp):
            key = tuple(map(int, sp)), tuple(map(int, ks))
            if key in result:
                raise ValueError("DUPLICATE_V2_KEY")
            result[key] = complex(z)
    return result


def _l2_norm(values):
    norm = 0.0
    for value in values:
        norm = math.hypot(norm, abs(value))
    return norm


def projection_rejected_norm(state, source, projection):
    norm = 0.0
    for key, amplitude in state.items():
        if projection.covariant_scalar_to_gauss_key(key, source) is None:
            norm = math.hypot(norm, abs(amplitude))
    return norm


def compare_states(a, b):
    keys = set(a) | set(b)
    delta = _l2_norm(a.get(key, 0j) - b.get(key, 0j) for key in keys)
    norm_a = _l2_norm(a.values())
    norm_b = _l2_norm(b.values())
    denominator = max(norm_a, norm_b)
    if not all(math.isfinite(value) for value in (delta, norm_a, norm_b)):
        relative = math.inf
    else:
        relative = (delta / denominator if denominator > 0.0
                    else (0.0 if delta == 0.0 else math.inf))
    return {"relative_delta": relative, "absolute_delta": delta,
            "norm_v4": norm_a, "norm_v2": norm_b,
            "support_v4": len(a), "support_v2": len(b)}


def _stable_complex_sum(values):
    """Accurately sum finite complex contributions using separate fsum channels."""
    real = math.fsum(float(complex(value).real) for value in values)
    imag = math.fsum(float(complex(value).imag) for value in values)
    return complex(real, imag)


def sum_states_exact(states):
    """Sum complex state maps without tolerance pruning or order-sensitive loss."""
    buckets = {}
    for state in states:
        for key, amplitude in state.items():
            buckets.setdefault(key, []).append(amplitude)
    total = {key: _stable_complex_sum(values)
             for key, values in buckets.items()}
    return {key: value for key, value in total.items() if value != 0j}


def sum_states_legacy_pruned(states, threshold):
    """Reproduce the historical ordered add(..., tol=threshold) semantics."""
    total = {}
    for state in states:
        for key, amplitude in state.items():
            value = total.get(key, 0j) + amplitude
            if abs(value) > threshold:
                total[key] = value
            elif key in total:
                del total[key]
    return total


def assemble_pair_exact(results, coefficients):
    """Assemble signed pair using compensated real/imaginary summation."""
    buckets = {}
    for name in sorted(coefficients):
        coefficient = coefficients[name]
        for key, amplitude in results[name].items():
            buckets.setdefault(key, []).append(coefficient * amplitude)
    total = {key: _stable_complex_sum(values)
             for key, values in buckets.items()}
    total = {key: value for key, value in total.items() if value != 0j}
    norm = _l2_norm(total.values())
    pair_status = ("VERIFIED_EXACT_ZERO_PAIR" if not total
                   else "VERIFIED_NONZERO_PAIR")
    return total, norm, pair_status


def collect(manifest):
    source_commit = manifest.get("source_commit")
    if not isinstance(source_commit, str) or not re.fullmatch(r"[0-9a-f]{40}", source_commit):
        raise ValueError("BAD_MANIFEST_SOURCE_COMMIT")
    try:
        actual_head = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, timeout=10
        ).strip()
    except (OSError, subprocess.SubprocessError) as exc:
        raise ValueError("COLLECTOR_COMMIT_NOT_VERIFIABLE") from exc
    if actual_head != source_commit:
        raise ValueError("COLLECTOR_SOURCE_COMMIT_MISMATCH")
    dirty = subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        text=True, timeout=10
    )
    if dirty.strip():
        raise ValueError("COLLECTOR_TRACKED_WORKTREE_NOT_CLEAN")
    if not _sha(manifest.get("seed_state_sha256")):
        raise ValueError("BAD_MANIFEST_SEED_DIGEST")

    backend_file = Path(__file__).with_name("tetrahedral_volume_backend.py")
    actual_backend_sha = hashlib.sha256(backend_file.read_bytes()).hexdigest()
    if actual_backend_sha != manifest.get("volume_backend_sha256"):
        raise ValueError("COLLECTOR_VOLUME_BACKEND_SHA_MISMATCH")

    from pl_dual_complex import DualComplex, seed_16cell_boundary
    from pl_lorentzian_physical_pair_worker_v4 import pair_indices
    import pl_lorentzian_triple_worker as W
    import peter_weyl_lorentzian_gauss_action_gate as LGA
    import peter_weyl_lorentzian_logical_projection_gate as LP

    D = DualComplex(seed_16cell_boundary())
    source, omit, cycle = (int(manifest[k]) for k in ("source_node", "omit", "cycle"))
    p, q, ip, iq = pair_indices(D, source, omit, cycle)
    targets_p = tuple(D.neighbor[(source, t)] for t in p)
    targets_q = tuple(D.neighbor[(source, t)] for t in q)
    shared = {k: manifest[k] for k in ("source_commit", "source_node", "omit", "cycle",
              "Jmax2", "volume_backend", "volume_backend_sha256", "seed_state_sha256")}
    if shared["volume_backend"] != "tetrahedral-volume-v1":
        raise ValueError("WRONG_MANIFEST_VOLUME_BACKEND_ID")
    if shared["Jmax2"] != 7:
        raise ValueError("WRONG_JMAX")
    paths = {}
    for item in manifest.get("paths", ()):
        key, state = certify_triple(item, shared, targets_p, targets_q)
        if key in paths:
            raise ValueError("DUPLICATE_AUX_PATH")
        paths[key] = state
    allkeys = set(product(("forward", "adjoint"), ("p", "q"), product((0, 1), repeat=3)))
    if set(paths) != allkeys:
        raise ValueError("V4_32_PATHS_NOT_COMPLETE")
    results, legacy_results, quality = {}, {}, {}
    for mode, order in product(("forward", "adjoint"), ("p", "q")):
        ordered_aux_states = [
            paths[mode, order, aux] for aux in product((0, 1), repeat=3)
        ]
        cov = sum_states_exact(ordered_aux_states)
        legacy_cov = sum_states_legacy_pruned(ordered_aux_states, LP.TOL)
        # Match the independent V2 scalar-projection tolerance, not the
        # looser display/worker threshold in pl_lorentzian_triple_worker.
        gauss, accepted2, rejected2 = LGA.project_scalar_gauss(cov, source, W.TOL)
        legacy_gauss, legacy_accepted2, legacy_rejected2 = LGA.project_scalar_gauss(
            legacy_cov, source, W.TOL)
        numeric_checks = (accepted2, rejected2, legacy_accepted2, legacy_rejected2)
        if not all(math.isfinite(float(value)) for value in numeric_checks):
            raise ValueError("NONFINITE_GAUSS_PROJECTION_NORM")
        rejected_norm = projection_rejected_norm(cov, source, LGA)
        legacy_rejected_norm = projection_rejected_norm(legacy_cov, source, LGA)
        if rejected_norm >= 1e-8 or legacy_rejected_norm >= 1e-8:
            raise ValueError("NONSCALAR_REJECTION_TOO_LARGE")
        if W.max_spin(gauss) > 3.5 + 1e-12:
            raise ValueError("PHYSICAL_SPIN_WALL")
        name = {("forward", "p"): "fp", ("forward", "q"): "fq",
                ("adjoint", "p"): "ap", ("adjoint", "q"): "aq"}[mode, order]
        results[name] = gauss
        legacy_results[name] = legacy_gauss
        quality[name] = {
            "unpruned_gauss_support": len(gauss),
            "unpruned_accepted2": float(accepted2),
            "unpruned_rejected2": float(rejected2),
            "unpruned_rejected_norm": rejected_norm,
            "legacy_pruned_covariant_support": len(legacy_cov),
            "legacy_gauss_support": len(legacy_gauss),
            "legacy_accepted2": float(legacy_accepted2),
            "legacy_rejected2": float(legacy_rejected2),
            "legacy_rejected_norm": legacy_rejected_norm,
            "legacy_auxiliary_prune_tolerance": float(LP.TOL),
            "gauss_projection_tolerance": float(W.TOL),
        }
    ref_specs = {"fp": ("forward", ip), "fq": ("forward", iq),
                 "ap": ("adjoint", ip), "aq": ("adjoint", iq)}
    if len({manifest["independent_v2"][k]["npz"] for k in ref_specs}) != 4:
        raise ValueError("V2_REFERENCES_NOT_DISTINCT_FILES")
    v2 = {k: read_v2(
              manifest["independent_v2"][k], expected_source=source,
              expected_mode=ref_specs[k][0], expected_index=ref_specs[k][1],
              expected_source_commit=source_commit,
              expected_backend_sha=shared["volume_backend_sha256"])
          for k in ref_specs}
    comparisons = {k: compare_states(results[k], v2[k]) for k in v2}
    legacy_comparisons = {k: compare_states(legacy_results[k], v2[k]) for k in v2}
    if any(x["relative_delta"] >= 1e-10 for x in legacy_comparisons.values()):
        raise ValueError("LEGACY_THRESHOLD_ALIGNED_V2_V4_COMPLEX_MISMATCH")
    cp = int(W.ordered_spec(D, source, ip)[4])
    cq = int(W.ordered_spec(D, source, iq)[4])
    if cp != -cq:
        raise ValueError("WRONG_PAIR_ORIENTATION")
    S, pair_norm, pair_status = assemble_pair_exact(
        results, {"fp": -.5j * cp, "fq": -.5j * cq,
                  "ap": +.5j * cp, "aq": +.5j * cq})
    if not math.isfinite(pair_norm):
        raise ValueError("NONFINITE_PHYSICAL_PAIR_NORM")
    return {"status": "VERIFIED_ORIGINAL_V4_ONE_PHYSICAL_PAIR", "passed": True,
            "paths_completed": 32, "v2_primitives_compared": 4,
            "unpruned_v2_comparisons": comparisons,
            "legacy_threshold_aligned_v2_comparisons": legacy_comparisons,
            "quality": quality,
            "gauss_projection_tolerance": float(W.TOL),
            "legacy_auxiliary_prune_tolerance": float(LP.TOL),
            "auxiliary_sum_policy": "report both unpruned stable fsum and legacy threshold-aligned sum",
            "pair_status": pair_status, "pair_norm": pair_norm,
            "pair_zero_semantics": "exact floating-point sum; not symbolic identity",
            "pair_support": len(S), "HDA_proved": False, "HL_complete": False}, S


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    try:
        result, state = collect(json.loads(args.manifest.read_text("utf-8")))
        import pl_lorentzian_triple_worker as W
        W.save_state(args.output_prefix.with_suffix(".npz"), state, 32, 16)
        result["npz_sha256"] = hashlib.sha256(args.output_prefix.with_suffix(".npz").read_bytes()).hexdigest()
        code = 0
    except Exception as exc:
        result = {"status": "V21_C1_NOT_COMPLETE", "passed": False, "HDA_proved": False,
                  "error_type": type(exc).__name__, "error": str(exc)}
        code = 2
    args.output_prefix.with_suffix(".json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
