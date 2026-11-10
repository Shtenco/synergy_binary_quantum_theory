#!/usr/bin/env python3
"""Capture the exact two-endpoint state at the first W5 scalar-covariance failure.

This probe reuses the real W5 KKV path, but restricts the input to one seed and
one source node. It does not alter amplitudes, thresholds, projections, or the
scalar-preservation guard. Its only purpose is to save the exact input key,
operator stage, and full output J,M distribution when the guard fires.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))


class ProbeCaptured(Exception):
    pass


def json_default(value):
    if hasattr(value, "item"):
        return value.item()
    if isinstance(value, (set, tuple)):
        return list(value)
    return repr(value)


def output_sector_snapshot(state, charged_nodes, original):
    by_node = {}
    for (_, labels), amp in state.items():
        for node in charged_nodes:
            lab = labels[node]
            if lab[0] != "C":
                label = repr(lab)
            else:
                _, j2, m2, _, _ = lab
                label = f"J2={int(j2)},M2={int(m2)}"
            by_node.setdefault(str(node), {})
            by_node[str(node)][label] = (
                by_node[str(node)].get(label, 0.0) + float(abs(amp) ** 2)
            )

    wrong = {}
    for node in charged_nodes:
        expected_j2, expected_m2 = original[node]
        rows = by_node.get(str(node), {})
        wrong[str(node)] = {
            label: weight
            for label, weight in rows.items()
            if label != f"J2={expected_j2},M2={expected_m2}"
        }
    return by_node, wrong


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--K", type=int, choices=(0, 2, 4, 6, 8), default=0)
    ap.add_argument("--Jmax2", type=int, default=10)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    import bqg_v10_w5_lorentzian_adjoint_closure as PILOT
    import peter_weyl_covariant_K_composition_cached_gate as CACHE
    import peter_weyl_covariant_K_composition_gate as KC
    import peter_weyl_covariant_K_sine_composition_gate as SINEK

    context = {}
    original_canonicalize = CACHE.canonicalize_scalar_charge_M
    original_restore = CACHE.restore_scalar_charge_M
    original_complete = SINEK.complete_HE_sine
    original_k = KC.apply_K_complete_custom
    original_he_state = KC.apply_HE_complete_state_custom
    original_ck = KC.C_K_component
    original_w5 = PILOT.W5_K

    def traced_canonicalize(key, charged_nodes):
        canonical, original = original_canonicalize(key, charged_nodes)
        context.update({
            "input_key_repr": repr(key),
            "canonical_key_repr": repr(canonical),
            "charged_nodes": list(charged_nodes),
            "original_JM": {str(k): list(v) for k, v in original.items()},
        })
        return canonical, original

    def traced_restore(state, charged_nodes, original):
        try:
            return original_restore(state, charged_nodes, original)
        except RuntimeError as exc:
            by_node, wrong = output_sector_snapshot(state, charged_nodes, original)
            ranked = sorted(state.items(), key=lambda item: abs(item[1]), reverse=True)
            detail = {
                "status": "CAPTURED_SCALAR_COVARIANCE_FAILURE",
                "K_seed": args.K,
                "Jmax2": args.Jmax2,
                **context,
                "source_v": context.get("source_v"),
                "he_stage": context.get("he_stage"),
                "K_call_kind": context.get("K_call_kind"),
                "K_call_index": context.get("K_call_index"),
                "C_K_arguments": context.get("C_K_arguments"),
                "original_exception": str(exc),
                "output_weight_by_node_JM": by_node,
                "wrong_sector_weight_by_node_JM": wrong,
                "output_state_support": len(state),
                "largest_output_components": [
                    {
                        "key_repr": repr(key),
                        "real": float(amp.real),
                        "imag": float(amp.imag),
                        "abs": float(abs(amp)),
                    }
                    for key, amp in ranked[:25]
                ],
                "operator_changed": False,
                "projection_changed": False,
                "guard_disabled": False,
            }
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(detail, indent=2, default=json_default) + "\n",
                encoding="utf-8",
            )
            print(json.dumps(detail, indent=2, default=json_default), flush=True)
            raise ProbeCaptured("Expected scalar-covariance failure captured") from exc

    def traced_complete(key, source_v, Jmax2, charged_nodes=(0, 1)):
        # Reproduce complete_HE_sine's exact loop and accumulation order while
        # observing each oriented specification, adjoint choice and T sequence.
        # The physical output accumulator receives precisely the same additions
        # in the same order as the original function.
        import k5_peter_weyl_safe_hda_column as PW

        context.update({
            "source_v": int(source_v),
            "Jmax2": int(Jmax2),
            "charged_nodes": list(charged_nodes),
            "canonical_key_repr": repr(key),
        })

        def jm_weights(state):
            result = {}
            for (_, labels), amp in state.items():
                for node in charged_nodes:
                    lab = labels[node]
                    if lab[0] == "C":
                        _, j2, m2, _, _ = lab
                        slot = result.setdefault(str(node), {})
                        label = f"J2={int(j2)},M2={int(m2)}"
                        slot[label] = slot.get(label, 0.0) + float(abs(amp) ** 2)
            return result

        def wrong_weight(weights):
            total = 0.0
            wrong = 0.0
            for node in charged_nodes:
                expected_j2, expected_m2 = context["original_JM"][str(node)]
                row = weights.get(str(node), {})
                for label, value in row.items():
                    total += value
                    if label != f"J2={expected_j2},M2={expected_m2}":
                        wrong += value
            return {"total_endpoint_weight": total, "wrong_endpoint_weight": wrong}

        base = KC.CK.branch_from_key(key)
        out = {}
        max_v = 0.0
        max_b = 0.0
        rows = []
        for spec_index, (sign, spec) in enumerate(PW.oriented_specs(source_v)):
            spec_acc = {}
            for adj in (False, True):
                adj_acc = {}
                for sequence_index, (coef, seq0) in enumerate(PW.T_sequences(*spec)):
                    seq = PW.adjoint_sequence(seq0) if adj else seq0
                    branches, vleak = KC.CK.apply_sequence_to_branch(
                        base, seq, source_v, Jmax2
                    )
                    max_v = max(max_v, float(vleak))
                    sequence_acc = {}
                    for branch_index, branch in enumerate(branches):
                        projected, bleak = KC.CK.project_branch_complete_charges(
                            branch, charged_nodes
                        )
                        max_b = max(max_b, float(bleak))
                        scale = (+0.5j if adj else -0.5j) * sign * coef
                        KC.CK.add_state(out, projected, scale)
                        KC.CK.add_state(sequence_acc, projected, scale)
                    weights = jm_weights(sequence_acc)
                    rows.append({
                        "spec_index": spec_index,
                        "oriented_spec": list(spec),
                        "epsilon_sign": int(sign),
                        "adjoint": bool(adj),
                        "sequence_index": sequence_index,
                        "sequence": repr(seq),
                        "sequence_coefficient": [float(complex(coef).real), float(complex(coef).imag)],
                        "branch_count": len(branches),
                        "sequence_output_support": len(sequence_acc),
                        "sequence_JM_weights": weights,
                        **wrong_weight(weights),
                    })
                    KC.CK.add_state(adj_acc, sequence_acc)
                KC.CK.add_state(spec_acc, adj_acc)
            spec_weights = jm_weights(spec_acc)
            rows.append({
                "spec_index": spec_index,
                "oriented_spec": list(spec),
                "epsilon_sign": int(sign),
                "stage": "after_complete_oriented_spec",
                "spec_output_support": len(spec_acc),
                "spec_JM_weights": spec_weights,
                **wrong_weight(spec_weights),
            })

        total_weights = jm_weights(out)
        context["operator_stage_trace"] = {
            "construction": "complete_HE_sine exact loop with passive per-sequence observers",
            "sequence_rows": rows,
            "total_output_support": len(out),
            "total_output_JM_weights": total_weights,
            "total_output_wrong_weight": wrong_weight(total_weights),
            "max_volume_leakage": max_v,
            "max_basis_leakage": max_b,
        }
        return out, max_v, max_b

    def traced_he_state(state, source_v, Jmax2, charged_nodes):
        idx = int(context.get("he_call_index", 0)) + 1
        context["he_call_index"] = idx
        context["he_stage"] = (
            "first_H_E_before_volume" if idx == 1
            else "second_H_E_after_volume" if idx == 2
            else f"H_E_call_{idx}"
        )
        return original_he_state(state, source_v, Jmax2, charged_nodes)

    def traced_k(state, source_v, Jmax2, charged_nodes):
        context["K_call_index"] = int(context.get("K_call_index", 0)) + 1
        context["K_call_kind"] = (
            "direct_K_single_charged_endpoint"
            if tuple(charged_nodes) == (source_v,)
            else "inverse_holonomy_two_charged_endpoints"
        )
        context["he_call_index"] = 0
        return original_k(state, source_v, Jmax2, charged_nodes)

    def traced_ck(state_cov, source_v, target_v, i, j, Jmax2):
        context["C_K_arguments"] = {
            "source_v": int(source_v), "target_v": int(target_v),
            "i": int(i), "j": int(j), "Jmax2": int(Jmax2),
        }
        return original_ck(state_cov, source_v, target_v, i, j, Jmax2)

    CACHE.canonicalize_scalar_charge_M = traced_canonicalize
    CACHE.restore_scalar_charge_M = traced_restore
    SINEK.complete_HE_sine = traced_complete
    KC.apply_K_complete_custom = traced_k
    KC.apply_HE_complete_state_custom = traced_he_state
    KC.C_K_component = traced_ck
    PILOT.W5_K = (args.K,)

    try:
        result = PILOT.run_closure(
            source_nodes=(0,),
            Jmax2=args.Jmax2,
            max_states=1,
            max_depth=0,
            amp_tol=1e-12,
        )
        detail = {
            "status": "NO_SCALAR_COVARIANCE_FAILURE_REPRODUCED",
            "K_seed": args.K,
            "Jmax2": args.Jmax2,
            "pilot_status": result.get("status"),
            "operator_changed": False,
            "projection_changed": False,
            "guard_disabled": False,
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(detail, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(detail, indent=2), flush=True)
    except ProbeCaptured:
        return 0
    finally:
        CACHE.canonicalize_scalar_charge_M = original_canonicalize
        CACHE.restore_scalar_charge_M = original_restore
        SINEK.complete_HE_sine = original_complete
        KC.apply_K_complete_custom = original_k
        KC.apply_HE_complete_state_custom = original_he_state
        KC.C_K_component = original_ck
        PILOT.W5_K = original_w5

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
