#!/usr/bin/env python3
"""Collect genuine complete 32-path V4 and compare with four independent V2."""
from __future__ import annotations
import argparse, hashlib, json, math
from itertools import product
from pathlib import Path
import numpy as np
from v21_stage_logic import operations, predecessor_ok
from v21_checkpoint_guard import validate_physical_stage
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
        meta["diagnostics"], state, meta.get("stage"), signature,
        meta.get("guard", {}).get("threshold", 1e-8))
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


def read_v2(spec, *, expected_source=None, expected_mode=None, expected_index=None):
    file = Path(spec["npz"])
    if hashlib.sha256(file.read_bytes()).hexdigest() != spec["sha256"]:
        raise ValueError("INDEPENDENT_V2_SHA256_MISMATCH")
    if expected_source is not None:
        meta_path = Path(spec["metadata_json"])
        if hashlib.sha256(meta_path.read_bytes()).hexdigest() != spec["metadata_sha256"]:
            raise ValueError("V2_METADATA_SHA256_MISMATCH")
        meta = json.loads(meta_path.read_text("utf-8"))
        if meta.get("passed") is not True or meta.get("mode") != expected_mode:
            raise ValueError("V2_PRIMITIVE_FAILED_OR_WRONG_MODE")
        if meta.get("index") != expected_index or meta.get("source_node") != expected_source:
            raise ValueError("V2_PRIMITIVE_WRONG_INPUT_OR_ORIENTATION")
        if float(meta.get("Jmax", -1)) != 3.5:
            raise ValueError("V2_PRIMITIVE_WRONG_CUTOFF")
        if not meta.get("checks") or not all(meta["checks"].values()):
            raise ValueError("V2_PRIMITIVE_FAILED_ACCEPTANCE")
    with np.load(file, allow_pickle=False) as data:
        if set(data.files) != {"spins", "Ks", "amp"}:
            raise ValueError("INVALID_V2_FORMAT")
        spins, Ks, amp = data["spins"], data["Ks"], data["amp"]
        if (spins.shape != (len(amp), 32) or Ks.shape != (len(amp), 16)
                or amp.ndim != 1 or amp.dtype.kind != "c" or not np.isfinite(amp).all()):
            raise ValueError("V2_BAD_SHAPE_OR_NONFINITE")
        result = {}
        for sp, ks, z in zip(spins, Ks, amp):
            key = tuple(map(int, sp)), tuple(map(int, ks))
            if key in result:
                raise ValueError("DUPLICATE_V2_KEY")
            result[key] = complex(z)
    return result


def compare_states(a, b):
    keys = set(a) | set(b)
    delta = math.sqrt(math.fsum(abs(a.get(k, 0j) - b.get(k, 0j)) ** 2 for k in keys))
    norm_a = math.sqrt(math.fsum(abs(z) ** 2 for z in a.values()))
    norm_b = math.sqrt(math.fsum(abs(z) ** 2 for z in b.values()))
    return {"relative_delta": delta / max(norm_a, norm_b, 1e-15),
            "absolute_delta": delta, "support_v4": len(a), "support_v2": len(b)}


def sum_states_exact(states):
    """Sum complex state maps without a numerical support threshold."""
    total = {}
    for state in states:
        for key, amplitude in state.items():
            total[key] = total.get(key, 0j) + amplitude
    return {key: value for key, value in total.items() if value != 0j}


def assemble_pair_exact(results, coefficients):
    """Assemble the signed pair without tolerance pruning or zero relabeling."""
    total = {}
    for name, coefficient in coefficients.items():
        for key, amplitude in results[name].items():
            total[key] = total.get(key, 0j) + coefficient * amplitude
    total = {key: value for key, value in total.items() if value != 0j}
    norm = math.sqrt(math.fsum(abs(z) ** 2 for z in total.values()))
    return total, norm, ("VERIFIED_EXACT_ZERO_PAIR" if norm == 0.0
                         else "VERIFIED_NONZERO_PAIR")


def collect(manifest):
    from pl_dual_complex import DualComplex, seed_16cell_boundary
    from pl_lorentzian_physical_pair_worker_v4 import pair_indices
    import pl_lorentzian_triple_worker as W
    import peter_weyl_lorentzian_gauss_action_gate as LGA

    D = DualComplex(seed_16cell_boundary())
    source, omit, cycle = (int(manifest[k]) for k in ("source_node", "omit", "cycle"))
    p, q, ip, iq = pair_indices(D, source, omit, cycle)
    targets_p = tuple(D.neighbor[(source, t)] for t in p)
    targets_q = tuple(D.neighbor[(source, t)] for t in q)
    shared = {k: manifest[k] for k in ("source_commit", "source_node", "omit", "cycle",
              "Jmax2", "volume_backend", "volume_backend_sha256", "seed_state_sha256")}
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
    results, quality = {}, {}
    for mode, order in product(("forward", "adjoint"), ("p", "q")):
        cov = sum_states_exact(
            paths[mode, order, aux] for aux in product((0, 1), repeat=3)
        )
        gauss, accepted2, rejected2 = LGA.project_scalar_gauss(cov, source, W.TOL)
        if math.sqrt(max(float(rejected2), 0)) >= 1e-8:
            raise ValueError("NONSCALAR_REJECTION_TOO_LARGE")
        if W.max_spin(gauss) > 3.5 + 1e-12:
            raise ValueError("PHYSICAL_SPIN_WALL")
        name = {("forward", "p"): "fp", ("forward", "q"): "fq",
                ("adjoint", "p"): "ap", ("adjoint", "q"): "aq"}[mode, order]
        results[name] = gauss
        quality[name] = {"gauss_support": len(gauss), "accepted2": float(accepted2),
                         "rejected2": float(rejected2)}
    ref_specs = {"fp": ("forward", ip), "fq": ("forward", iq),
                 "ap": ("adjoint", ip), "aq": ("adjoint", iq)}
    if len({manifest["independent_v2"][k]["npz"] for k in ref_specs}) != 4:
        raise ValueError("V2_REFERENCES_NOT_DISTINCT_FILES")
    v2 = {k: read_v2(manifest["independent_v2"][k], expected_source=source,
                     expected_mode=ref_specs[k][0], expected_index=ref_specs[k][1])
          for k in ref_specs}
    comparisons = {k: compare_states(results[k], v2[k]) for k in v2}
    if any(x["relative_delta"] >= 1e-10 for x in comparisons.values()):
        raise ValueError("ORIGINAL_V2_V4_COMPLEX_MISMATCH")
    cp = int(W.ordered_spec(D, source, ip)[4])
    cq = int(W.ordered_spec(D, source, iq)[4])
    if cp != -cq:
        raise ValueError("WRONG_PAIR_ORIENTATION")
    S, pair_norm, pair_status = assemble_pair_exact(
        results, {"fp": -.5j * cp, "fq": -.5j * cq,
                  "ap": +.5j * cp, "aq": +.5j * cq})
    return {"status": "VERIFIED_ORIGINAL_V4_ONE_PHYSICAL_PAIR", "passed": True,
            "paths_completed": 32, "v2_primitives_compared": 4,
            "comparisons": comparisons, "quality": quality,
            "pair_status": pair_status, "pair_norm": pair_norm,
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
