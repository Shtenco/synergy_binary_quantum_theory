#!/usr/bin/env python3
"""Paired V09/V10 W5 regression on identical seeds, preserving raw amplitudes.

This is an audit harness, not a physics PASS generator. It runs the frozen V09
column adapter and V10 finite closure on the same node, five K seeds and J wall.
It records pre-serialization V09 maps, the legacy serializer loss, epsilon-orbit
identity, V10 direct-adjoint residuals, and an old/new sine-composition probe.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import math
from pathlib import Path


W5_K = (0, 2, 4, 6, 8)
SERIALIZER_CUTOFF = 1e-13


def _norm(state):
    value = 0.0
    for amplitude in state.values():
        value = math.hypot(value, abs(complex(amplitude)))
    return value


def _relative_difference(a, b):
    keys = set(a) | set(b)
    err = 0.0
    for key in keys:
        err = math.hypot(err, abs(complex(a.get(key, 0j)) - complex(b.get(key, 0j))))
    scale = _norm(b)
    return {
        "absolute_l2": err,
        "reference_l2": scale,
        "relative_l2": err / max(scale, 1e-300),
        "union_support": len(keys),
        "left_support": len(a),
        "right_support": len(b),
    }


def _decode_v10_column(result, seed_key):
    basis = result["basis_keys_repr"]
    col_index = basis.index(repr(seed_key))
    encoded = result["forward_B_columns_sparse"][col_index]
    out = {}
    for row_index, real, imag in encoded:
        key_repr = basis[int(row_index)]
        out[key_repr] = complex(float(real), float(imag))
    return out


def _repr_map(state):
    return {repr(key): complex(value) for key, value in state.items() if complex(value) != 0j}


def _sparse_digest(state):
    digest = hashlib.sha256()
    for key, amplitude in sorted(state.items(), key=lambda item: repr(item[0])):
        z = complex(amplitude)
        digest.update(repr(key).encode("utf-8"))
        digest.update(b"\0")
        digest.update(z.real.hex().encode("ascii"))
        digest.update(b",")
        digest.update(z.imag.hex().encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _matrix_digest(basis, columns):
    payload = json.dumps(
        {"basis": basis, "columns": columns},
        sort_keys=True, separators=(",", ":"), allow_nan=False,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _serializer_loss(raw, serialized_rows):
    kept = {row[0] for row in serialized_rows}
    dropped = {key: value for key, value in raw.items() if repr(key) not in kept}
    return {
        "input_support": len(raw),
        "serialized_support": len(serialized_rows),
        "dropped_support": len(dropped),
        "dropped_l2": _norm(dropped),
        "max_dropped_abs": max((abs(z) for z in dropped.values()), default=0.0),
        "cutoff": SERIALIZER_CUTOFF,
    }


def _sine_composition_probe(node, seed_k, jmax2, mode):
    import peter_weyl_euclidean_sine_ordering_gate as OLD_SINE
    import peter_weyl_covariant_K_sine_composition_gate as NEW_SINE
    import peter_weyl_covariant_volume_leg_gate as CV
    import peter_weyl_covariant_composition_gate as COMP
    import bqg_pruning_audit as PRUNE
    import sys
    from pathlib import Path
    archive_dir = Path(__file__).resolve().parent / "archive"
    if str(archive_dir) not in sys.path:
        sys.path.insert(0, str(archive_dir))
    import bqg_v09_w5_lorentzian_raw_column_v09_frozen as V09

    seed = V09.w5_input_key(seed_k)
    context = PRUNE.bypass() if mode == "shadow_no_prune" else contextlib.nullcontext()
    with context:
        old, old_volume_leak = OLD_SINE.allJ_H_sine(seed, node, jmax2)
        covariant = CV.gauss_to_covariant({seed: 1.0 + 0j}, node)
        complete = COMP.covariant_to_complete(covariant, node)
        complete_key = next(iter(complete))
        diagnostics = {}
        new, new_volume_leak, new_basis_leak = NEW_SINE.complete_HE_sine(
            complete_key, node, jmax2, charged_nodes=(node,),
            projection_diagnostics=diagnostics,
        )

    new_scalar = {}
    new_nonscalar = {}
    for (spins, labels), amplitude in new.items():
        Klabels = []
        scalar = True
        for label in labels:
            if label[0] == "G":
                Klabels.append(int(label[1]))
            elif (label[0] == "C" and int(label[1]) == 0
                  and int(label[2]) == 0 and int(label[3]) == int(label[4])):
                Klabels.append(int(label[3]))
            else:
                scalar = False
                break
        target = new_scalar if scalar else new_nonscalar
        if scalar:
            key = (spins, tuple(Klabels))
            target[key] = target.get(key, 0j) + complex(amplitude)
        else:
            target[(spins, labels)] = complex(amplitude)

    old_repr = _repr_map(old)
    new_repr = _repr_map(new_scalar)
    return {
        "node": node,
        "seed_K": seed_k,
        "Jmax2": jmax2,
        "legacy_allJ_projection_support": len(old),
        "new_term_local_complete_support": len(new),
        "new_scalar_reconstructible_support": len(new_scalar),
        "new_non_scalar_or_unclosed_support": len(new_nonscalar),
        "legacy_vs_new_scalar_map": _relative_difference(old_repr, new_repr),
        "legacy_norm": _norm(old),
        "new_scalar_norm": _norm(new_scalar),
        "new_non_scalar_norm": _norm(new_nonscalar),
        "old_volume_leakage": float(old_volume_leak),
        "new_volume_leakage": float(new_volume_leak),
        "new_complete_basis_leakage": float(new_basis_leak),
        "new_projection_diagnostics": diagnostics,
        "interpretation": (
            "A mismatch is diagnostic: it must be localized to changed charge support, "
            "projection boundary, or explicit non-scalar output; it is not silently "
            "reclassified as an equivalent sine operator."
        ),
    }


def run(node, jmax2, mode, output):
    if node not in range(5) or jmax2 < 10:
        raise ValueError("node must be 0..4 and Jmax2 must be >=10")
    if mode not in ("production", "shadow_no_prune"):
        raise ValueError("invalid mode")

    import bqg_pruning_audit as PRUNE
    import k5_peter_weyl_safe_hda_column as PW
    import sys
    from pathlib import Path
    archive_dir = Path(__file__).resolve().parent / "archive"
    if str(archive_dir) not in sys.path:
        sys.path.insert(0, str(archive_dir))
    import bqg_v09_w5_lorentzian_raw_column_v09_frozen as V09
    import peter_weyl_lorentzian_sine_ordered_triple_gate as LSINE
    import bqg_v10_w5_lorentzian_adjoint_closure as V10

    v09_epsilon = V09.ordered_epsilon_triples(PW.NEIG[node])
    v10_epsilon = V10.epsilon_triples(PW.NEIG[node])
    epsilon_equal = v09_epsilon == v10_epsilon
    epsilon_audit = {
        "V09_count": len(v09_epsilon),
        "V10_count": len(v10_epsilon),
        "V09_unique_ordered_triples": len({p for _, p in v09_epsilon}),
        "V10_unique_ordered_triples": len({p for _, p in v10_epsilon}),
        "ordered_sign_and_permutation_lists_equal": epsilon_equal,
        "V09_epsilon_sign_sum": sum(sign for sign, _ in v09_epsilon),
        "V10_epsilon_sign_sum": sum(sign for sign, _ in v10_epsilon),
    }

    original_serializer = V09.serialize_sparse
    original_lsine_run = LSINE.run
    captured = []
    sine_capture = {}

    def capture_serializer(state):
        captured.append(dict(state))
        return original_serializer(state)

    def capture_lsine_run(*args, **kwargs):
        result = original_lsine_run(*args, **kwargs)
        sine_capture.clear()
        sine_capture.update(result)
        return result

    V09.serialize_sparse = capture_serializer
    LSINE.run = capture_lsine_run
    v09_rows = []
    try:
        context = PRUNE.bypass() if mode == "shadow_no_prune" else contextlib.nullcontext()
        with context:
            for seed_k in W5_K:
                captured.clear()
                sine_capture.clear()
                result = V09.run_column(node, seed_k, jmax2, 0, 24)
                if len(captured) != 3:
                    raise AssertionError(f"expected three V09 serialization calls, got {len(captured)}")
                raw_A, raw_B_gauss, raw_B_covariant = [dict(x) for x in captured]
                v09_rows.append({
                    "seed_K": seed_k,
                    "status": result["status"],
                    "epsilon_terms_completed": result["epsilon_terms_completed"],
                    "term_indices": [int(row["term_index"]) for row in result["term_diagnostics"]],
                    "term_epsilon_order_matches": [
                        (int(row["epsilon_sign"]), tuple(row["ordered_edges"]))
                        for row in result["term_diagnostics"]
                    ] == v09_epsilon,
                    "sine_wrapper_gate_passed": sine_capture.get("passed") is True,
                    "sine_wrapper_diagnostics": sine_capture.get("max_diagnostics", {}),
                    "A_serializer_loss": _serializer_loss(raw_A, result["A_Euclidean_sine_Gauss"]),
                    "B_gauss_serializer_loss": _serializer_loss(raw_B_gauss, result["B_Lorentzian_raw_KKV_Gauss"]),
                    "B_covariant_serializer_loss": _serializer_loss(raw_B_covariant, result["B_Lorentzian_raw_KKV_covariant"]),
                    "raw_B_gauss": _repr_map(raw_B_gauss),
                    "serialized_B_gauss": {row[0]: complex(float(row[1]), float(row[2])) for row in result["B_Lorentzian_raw_KKV_Gauss"]},
                })
    finally:
        V09.serialize_sparse = original_serializer
        LSINE.run = original_lsine_run

    amp_tol = 0.0 if mode == "shadow_no_prune" else 1e-12
    v10 = V10.run_closure(
        source_nodes=(node,),
        Jmax2=jmax2,
        max_states=500,
        max_depth=8,
        amp_tol=amp_tol,
        pruning_mode=mode,
        spectral_zero_factor=1000.0,
        cache_mode="off",
        seed_Ks=W5_K,
    )
    closure_complete = "forward_B_columns_sparse" in v10 and "basis_keys_repr" in v10
    comparisons = []
    if closure_complete:
        for row in v09_rows:
            seed_key = V10.w5_input_key(row["seed_K"])
            v10_column = _decode_v10_column(v10, seed_key)
            comparisons.append({
                "seed_K": row["seed_K"],
                "V09_raw_pre_serialization_vs_V10": _relative_difference(row["raw_B_gauss"], v10_column),
                "V09_serialized_vs_V10": _relative_difference(row["serialized_B_gauss"], v10_column),
                "V09_raw_support": len(row["raw_B_gauss"]),
                "V09_raw_norm": _norm(row["raw_B_gauss"]),
                "V09_raw_digest": _sparse_digest(row["raw_B_gauss"]),
                "V09_serialized_digest": _sparse_digest(row["serialized_B_gauss"]),
                "V10_seed_column_support": len(v10_column),
                "V10_seed_column_norm": _norm(v10_column),
                "V10_seed_column_digest": _sparse_digest(v10_column),
            })

    sine_rows = [
        _sine_composition_probe(node, seed_k, jmax2, mode)
        for seed_k in W5_K
    ]

    report = {
        "schema": "bqg-w5-v09-v10-paired-regression-v1",
        "status": "COMPUTED_DIAGNOSTICS_NOT_PHYSICAL_PASS",
        "source_node": node,
        "seed_Ks": list(W5_K),
        "Jmax2": jmax2,
        "mode": mode,
        "amp_tol": amp_tol,
        "cache_mode": "off",
        "epsilon_audit": epsilon_audit,
        "V09_columns": [
            {k: v for k, v in row.items() if k not in ("raw_B_gauss", "serialized_B_gauss")}
            for row in v09_rows
        ],
        "V10": {
            "status": v10.get("status"),
            "algebraic_status": v10.get("algebraic_status"),
            "closure_complete": closure_complete,
            "closure_size": v10.get("closure_size", v10.get("closure_size_discovered")),
            "basis_id_sha256": v10.get("basis_id_sha256"),
            "forward_matrix_digest": (
                _matrix_digest(v10.get("basis_keys_repr", []), v10.get("forward_B_columns_sparse", []))
                if closure_complete else None
            ),
            "direct_adjoint_matrix_digest": (
                _matrix_digest(v10.get("basis_keys_repr", []), v10.get("direct_Bdag_columns_sparse", []))
                if closure_complete else None
            ),
            "diagnostics": v10.get("diagnostics", v10.get("partial_diagnostics")),
            "coefficient_loss_ledger": v10.get("coefficient_loss_ledger"),
            "lower_level_pruning_ledger": v10.get("lower_level_pruning_ledger"),
        },
        "V09_vs_V10_forward_columns": comparisons,
        "old_vs_new_sine_composition": sine_rows,
        "acceptance": {
            "epsilon_orbits_identical": epsilon_equal and len(v09_epsilon) == 24,
            "all_five_V09_columns_complete": all(r["epsilon_terms_completed"] == 24 for r in v09_rows),
            "all_five_V09_term_orders_match": all(r["term_epsilon_order_matches"] for r in v09_rows),
            "all_five_V09_sine_wrapper_gates_pass": all(r["sine_wrapper_gate_passed"] for r in v09_rows),
            "V10_common_basis_complete": closure_complete,
            "V10_algebraic_closure_passed": v10.get("algebraic_status") == "PASS_ADJOINT_CONSISTENCY",
            "direct_adjoint_residual_below_1e-8": (
                isinstance(v10.get("diagnostics"), dict)
                and isinstance(v10["diagnostics"].get("direct_adjoint_relative_residual"), (int, float))
                and v10["diagnostics"]["direct_adjoint_relative_residual"] < 1e-8
            ),
            "V09_V10_forward_columns_match_below_1e-8": (
                len(comparisons) == 5
                and all(x["V09_raw_pre_serialization_vs_V10"]["relative_l2"] < 1e-8 for x in comparisons)
            ),
            "all_sine_composition_probes_match_below_1e-8": all(
                x["legacy_norm"] > 1e-10
                and x["new_scalar_norm"] > 1e-10
                and x["legacy_vs_new_scalar_map"]["relative_l2"] < 1e-8
                and x["new_non_scalar_norm"] < 1e-9
                for x in sine_rows
            ),
        },
        "scope": (
            "Finite raw W5 KKV structural comparison only. No physical Hermitian H_L, "
            "beta/kappa/hbar prefactor, physical joint-kernel, or continuum claim."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, default=lambda x: [x.real, x.imag] if isinstance(x, complex) else str(x)) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "node": node,
        "Jmax2": jmax2,
        "mode": mode,
        "epsilon_equal": epsilon_equal,
        "V10_status": v10.get("status"),
        "closure_complete": closure_complete,
        "comparison_count": len(comparisons),
        "acceptance": report["acceptance"],
        "output": str(output),
    }, indent=2))
    return 0 if all(report["acceptance"].values()) else 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--node", type=int, default=0)
    parser.add_argument("--Jmax2", type=int, default=10)
    parser.add_argument("--mode", choices=("production", "shadow_no_prune"), default="shadow_no_prune")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    return run(args.node, args.Jmax2, args.mode, args.output)


if __name__ == "__main__":
    raise SystemExit(main())
