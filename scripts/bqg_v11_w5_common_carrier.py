#!/usr/bin/env python3
"""Assemble a common W5 Euclidean/Lorentzian carrier from existing V09/V10 engines.

This is a structural integration builder. It deliberately does not choose the
physical Lorentzian coefficient, the unique Thiemann ordering, or a physical
joint-kernel verdict.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import importlib
import json
import math
import sys
import traceback
from pathlib import Path

W5_K = (0, 2, 4, 6, 8)
NODES = (0, 1, 2, 3, 4)


def _key_repr(key):
    return repr(key)


def _decode_v09_sparse(rows, label):
    out = {}
    for row in rows:
        if len(row) != 3:
            raise ValueError(f"{label}: sparse row must be [key_repr, real, imag]")
        key = ast.literal_eval(row[0])
        if not (isinstance(key, tuple) and len(key) == 2):
            raise ValueError(f"{label}: malformed Gauss key {row[0]!r}")
        z = complex(float(row[1]), float(row[2]))
        if not (math.isfinite(z.real) and math.isfinite(z.imag)):
            raise ValueError(f"{label}: non-finite amplitude")
        out[key] = out.get(key, 0j) + z
    return {k: z for k, z in out.items() if z != 0j}


def _decode_v10_column(encoded_column, basis_keys, label):
    out = {}
    for item in encoded_column:
        if len(item) != 3:
            raise ValueError(f"{label}: V10 sparse entry must be [row, real, imag]")
        row = int(item[0])
        if row < 0 or row >= len(basis_keys):
            raise ValueError(f"{label}: row index outside V10 basis")
        z = complex(float(item[1]), float(item[2]))
        if not (math.isfinite(z.real) and math.isfinite(z.imag)):
            raise ValueError(f"{label}: non-finite amplitude")
        key = ast.literal_eval(basis_keys[row])
        out[key] = out.get(key, 0j) + z
    return {k: z for k, z in out.items() if z != 0j}


def _add(dst, src):
    for key, z in src.items():
        dst[key] = dst.get(key, 0j) + z


def _norm(state):
    return math.sqrt(sum(abs(z) ** 2 for z in state.values()))


def _relative_residual(a, b):
    keys = set(a) | set(b)
    delta = math.sqrt(sum(abs(a.get(k, 0j) - b.get(k, 0j)) ** 2 for k in keys))
    scale = max(_norm(a), _norm(b), 1e-30)
    return delta / scale, delta, scale


def _encode_map(state, row_index):
    return [
        [row_index[key], float(z.real), float(z.imag)]
        for key, z in sorted(state.items(), key=lambda kv: row_index[kv[0]])
        if z != 0j
    ]


def build_carrier(
    Jmax2=10,
    source_nodes=NODES,
    max_states=500,
    max_depth=8,
    amp_tol=1e-12,
    crosscheck_tol=1e-8,
    leakage_tol=1e-9,
    v09_module=None,
    v10_module=None,
):
    nodes = tuple(sorted(set(int(v) for v in source_nodes)))
    if nodes != tuple(source_nodes) or not nodes or any(v not in NODES for v in nodes):
        raise ValueError("source_nodes must be a sorted, nonempty subset of 0..4")
    if Jmax2 < 10 or max_states < len(W5_K) or max_depth < 0:
        raise ValueError("invalid W5 cutoff or closure limits")
    if amp_tol <= 0 or crosscheck_tol <= 0 or leakage_tol <= 0:
        raise ValueError("amp_tol, crosscheck_tol and leakage_tol must be positive")

    here = Path(__file__).resolve().parent
    if str(here) not in sys.path:
        sys.path.insert(0, str(here))
    v09 = v09_module or importlib.import_module("bqg_v09_w5_lorentzian_raw_column")
    v10 = v10_module or importlib.import_module("bqg_v10_w5_lorentzian_adjoint_closure")

    # First obtain the full, two-sided closed domain and its direct-adjoint
    # matrices. This is the authoritative common input/output key universe
    # for the raw node-summed KKV map.
    closure = v10.run_closure(
        source_nodes=nodes,
        Jmax2=Jmax2,
        max_states=max_states,
        max_depth=max_depth,
        amp_tol=amp_tol,
    )
    if closure.get("status") != "PASS_ADJOINT_CONSISTENCY":
        return {
            "status": "BLOCKED_V10_CLOSURE",
            "Jmax2": Jmax2,
            "source_nodes": list(nodes),
            "closure_status": closure.get("status"),
            "closure_diagnostics": closure.get("diagnostics", closure.get("partial_diagnostics", {})),
            "claim_scope": "no common carrier emitted because V10 did not certify a complete adjoint-closed basis",
        }

    basis_repr = closure.get("basis_keys_repr")
    if not isinstance(basis_repr, list) or not basis_repr:
        raise ValueError("V10 closure is missing basis_keys_repr")
    basis_keys = [ast.literal_eval(x) for x in basis_repr]
    if len(set(basis_keys)) != len(basis_keys):
        raise ValueError("V10 basis contains duplicate keys")
    row_index = {key: i for i, key in enumerate(basis_keys)}
    seed_keys = [v10.w5_input_key(K) for K in W5_K]
    if any(k not in row_index for k in seed_keys):
        raise ValueError("V10 closed basis does not contain all five W5 seeds")

    forward_cols = closure["forward_B_columns_sparse"]
    sym_cols = closure["hermitian_ordering_candidate"]["B_comm_symmetric_columns"]
    anti_cols = closure["hermitian_ordering_candidate"]["B_comm_antihermitian_columns"]
    if (len(forward_cols) != len(basis_keys)
            or len(sym_cols) != len(basis_keys)
            or len(anti_cols) != len(basis_keys)):
        raise ValueError("V10 matrix column count does not match its basis")

    A_by_K = {K: {} for K in W5_K}
    B_v09_by_K = {K: {} for K in W5_K}
    v09_diagnostics = []
    for K in W5_K:
        for node in nodes:
            result = v09.run_column(node, K, Jmax2=Jmax2, start=0, count=24)
            if result.get("status") != "COMPUTED_RAW_W5_KKV_COLUMN":
                raise RuntimeError(f"V09 incomplete column: node={node}, K={K}, status={result.get('status')}")
            if result.get("epsilon_terms_completed") != 24:
                raise RuntimeError(f"V09 did not complete all 24 epsilon terms: node={node}, K={K}")
            if result.get("B_scalar_projection_rejections"):
                raise RuntimeError(f"V09 scalar projection rejected terms: node={node}, K={K}")
            max_term_leakage = max(
                (float(x["max_basis_leakage"]) for x in result["term_diagnostics"]),
                default=0.0,
            )
            non_scalar_ratio = float(result["B_non_scalar_covariant_norm"]) / max(
                float(result["B_raw_norm"]), 1e-30
            )
            if max_term_leakage > leakage_tol:
                raise RuntimeError(
                    f"V09 component leakage {max_term_leakage} exceeds {leakage_tol}: node={node}, K={K}"
                )
            max_spin_reached = max(
                (max(ast.literal_eval(row[0])[0]) for row in result["B_Lorentzian_raw_KKV_covariant"]),
                default=0,
            ) / 2
            if non_scalar_ratio > leakage_tol:
                raise RuntimeError(
                    f"V09 non-scalar covariant ratio {non_scalar_ratio} exceeds {leakage_tol}: node={node}, K={K}"
                )
            if max_spin_reached > 3.5 + 1e-12:
                raise RuntimeError(
                    f"V09 max spin {max_spin_reached} exceeds the frozen 3.5 wall: node={node}, K={K}"
                )
            a = _decode_v09_sparse(result["A_Euclidean_sine_Gauss"], "V09 A")
            b = _decode_v09_sparse(result["B_Lorentzian_raw_KKV_Gauss"], "V09 B")
            # V10 closes/stores support after the declared amp_tol. Apply the
            # same post-node-sum threshold before the independent comparison.
            b = {key: z for key, z in b.items() if abs(z) > amp_tol}
            _add(A_by_K[K], a)
            _add(B_v09_by_K[K], b)
            v09_diagnostics.append({
                "node": node,
                "input_K": K,
                "epsilon_terms_completed": result["epsilon_terms_completed"],
                "A_support": result["A_support"],
                "B_support": result["B_support"],
                "A_norm": result["A_norm"],
                "B_raw_norm": result["B_raw_norm"],
                "B_non_scalar_covariant_norm": result["B_non_scalar_covariant_norm"],
                "max_term_leakage": max_term_leakage,
                "non_scalar_covariant_ratio": non_scalar_ratio,
                "max_spin_reached": max_spin_reached,
            })

    B_forward_by_K = {}
    B_sym_by_K = {}
    B_anti_by_K = {}
    crosschecks = []
    for K, seed in zip(W5_K, seed_keys):
        col = row_index[seed]
        b_forward = _decode_v10_column(forward_cols[col], basis_repr, "V10 forward")
        b_sym = _decode_v10_column(sym_cols[col], basis_repr, "V10 symmetric")
        b_anti = _decode_v10_column(anti_cols[col], basis_repr, "V10 anti-Hermitian")
        B_forward_by_K[K] = b_forward
        B_sym_by_K[K] = b_sym
        B_anti_by_K[K] = b_anti
        residual, absolute, scale = _relative_residual(B_v09_by_K[K], b_forward)
        crosschecks.append({
            "input_K": K,
            "v09_node_sum_vs_v10_forward_relative_residual": residual,
            "absolute_residual": absolute,
            "reference_scale": scale,
            "passed": residual <= crosscheck_tol,
        })

    decomposition_checks = []
    for K in W5_K:
        sym_plus_anti = dict(B_sym_by_K[K])
        _add(sym_plus_anti, B_anti_by_K[K])
        residual, absolute, scale = _relative_residual(B_forward_by_K[K], sym_plus_anti)
        decomposition_checks.append({
            "input_K": K,
            "B_forward_vs_symmetric_plus_anti_relative_residual": residual,
            "absolute_residual": absolute,
            "reference_scale": scale,
            "passed": residual <= crosscheck_tol,
        })

    failed = [x for x in crosschecks if not x["passed"]]
    failed_decomposition = [x for x in decomposition_checks if not x["passed"]]
    if failed or failed_decomposition:
        return {
            "status": "FAIL_V09_V10_RAW_AMPLITUDE_CROSSCHECK",
            "Jmax2": Jmax2,
            "source_nodes": list(nodes),
            "crosscheck_tolerance": crosscheck_tol,
            "leakage_tolerance": leakage_tol,
            "crosschecks": crosschecks,
            "decomposition_checks": decomposition_checks,
            "v10_closure_basis_id_sha256": closure["basis_id_sha256"],
            "claim_scope": "raw engine mismatch; no common carrier or joint-Gram input emitted",
        }

    # The common output basis is the exact union of all Euclidean A supports
    # and V10's symmetrized raw Lorentzian supports. It is sorted by repr for
    # deterministic serialization; exact tuple keys and complex phases survive.
    common_keys = set()
    for K in W5_K:
        common_keys.update(A_by_K[K])
        common_keys.update(B_sym_by_K[K])
        common_keys.update(B_forward_by_K[K])
        common_keys.update(B_anti_by_K[K])
    common_keys = sorted(common_keys, key=repr)
    if not common_keys:
        raise ValueError("common W5 output carrier is empty")
    common_id = hashlib.sha256(
        "\n".join(repr(k) for k in common_keys).encode("utf-8")
    ).hexdigest()
    common_index = {key: i for i, key in enumerate(common_keys)}

    A_matrix = []
    Bsym_matrix = []
    Braw_matrix = []
    Banti_matrix = []
    for key in common_keys:
        A_matrix.append([
            [float(A_by_K[K].get(key, 0j).real), float(A_by_K[K].get(key, 0j).imag)]
            for K in W5_K
        ])
        Bsym_matrix.append([
            [float(B_sym_by_K[K].get(key, 0j).real), float(B_sym_by_K[K].get(key, 0j).imag)]
            for K in W5_K
        ])
        Braw_matrix.append([
            [float(B_forward_by_K[K].get(key, 0j).real), float(B_forward_by_K[K].get(key, 0j).imag)]
            for K in W5_K
        ])
        Banti_matrix.append([
            [float(B_anti_by_K[K].get(key, 0j).real), float(B_anti_by_K[K].get(key, 0j).imag)]
            for K in W5_K
        ])

    seed_payload = "\n".join(repr(v10.w5_input_key(K)) for K in W5_K).encode("utf-8")
    source_names = (
        "bqg_v09_w5_lorentzian_raw_column.py",
        "bqg_v10_w5_lorentzian_adjoint_closure.py",
        "peter_weyl_lorentzian_sine_ordered_triple_gate.py",
        "peter_weyl_lorentzian_ordered_triple_gate.py",
        "peter_weyl_euclidean_sine_ordering_gate.py",
        "peter_weyl_covariant_volume_leg_gate.py",
        "peter_weyl_zeroaware_volume_migration_experiment.py",
        "k5_peter_weyl_safe_hda_column.py",
    )
    source_hashes = {
        name: hashlib.sha256((here / name).read_bytes()).hexdigest()
        for name in source_names
    }

    result = {
        "status": "PASS_STRUCTURAL_COMMON_CARRIER",
        "claim_scope": "finite W5 structural carrier; raw Lorentzian KKV symmetric completion only; not a physical Hamiltonian or physical kernel",
        "input_basis_K": list(W5_K),
        "source_nodes": list(nodes),
        "Jmax2": Jmax2,
        "Jmax": Jmax2 / 2,
        "amp_tol": amp_tol,
        "crosscheck_tolerance": crosscheck_tol,
        "leakage_tolerance": leakage_tol,
        "node_weights": {str(node): 1.0 for node in nodes},
        "node_sum_convention": "unweighted coherent sum; unit lapse at every selected K5 source node",
        "input_seed_id_sha256": hashlib.sha256(seed_payload).hexdigest(),
        "source_modules_sha256": source_hashes,
        "v10_closure": {
            "basis_id_sha256": closure["basis_id_sha256"],
            "closure_size": closure["closure_size"],
            "closure_depth_max": closure["closure_depth_max"],
            "direct_adjoint_relative_residual": closure["diagnostics"]["direct_adjoint_relative_residual"],
            "max_component_leakage": closure["diagnostics"]["max_component_leakage"],
            "symmetric_completion_definition": "0.5*(B_comm+B_comm_dagger)",
            "physical_ordering_selected": False,
        },
        "common_output_basis": {
            "id_sha256": common_id,
            "size": len(common_keys),
            "keys_repr": [repr(k) for k in common_keys],
        },
        "A": A_matrix,
        "B": Bsym_matrix,
        "B_raw_forward": Braw_matrix,
        "B_antihermitian": Banti_matrix,
        "normalization": {
            "A": "V09 sum over declared source nodes; Euclidean sine-Hermitian structural map",
            "B": "V10 direct-adjoint symmetric completion of raw epsilon KKV map",
            "c_E": None,
            "c_L": None,
            "normalization_source": None,
            "physical_absolute_prefactor_fixed": False,
        },
        "crosschecks": {
            "V09_per_node_then_node_sum_vs_V10_all_node_forward": crosschecks,
            "V10_forward_equals_symmetric_plus_anti": decomposition_checks,
            "all_passed": True,
        },
        "v09_column_diagnostics": v09_diagnostics,
        "physical_status": {
            "unique_hermitian_ordering_selected": False,
            "absolute_beta_kappa_hbar_regulator_prefactor_fixed": False,
            "full_physical_HL": False,
            "physical_joint_kernel_verdict": False,
            "physical_causal_cone": False,
        },
    }
    return result


def joint_gram_input(carrier, c_E, c_L, normalization_source):
    if carrier.get("status") != "PASS_STRUCTURAL_COMMON_CARRIER":
        raise ValueError("only a passed structural common carrier can be adapted")
    if not isinstance(normalization_source, str) or not normalization_source.strip():
        raise ValueError("normalization_source must be explicit")
    return {
        "input_basis_K": carrier["input_basis_K"],
        "output_basis_id": carrier["common_output_basis"]["id_sha256"],
        "A_output_basis_id": carrier["common_output_basis"]["id_sha256"],
        "B_output_basis_id": carrier["common_output_basis"]["id_sha256"],
        "A": carrier["A"],
        "B": carrier["B"],
        "normalization": {"c_E": c_E, "c_L": c_L},
        "normalization_source": normalization_source,
        "claim_scope": "caller-declared finite W5 combination; physical only if normalization_source and coefficients are independently derived from the frozen canonical operator convention",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--Jmax2", type=int, default=10)
    ap.add_argument("--nodes", type=int, nargs="+", default=list(NODES))
    ap.add_argument("--max-states", type=int, default=500)
    ap.add_argument("--max-depth", type=int, default=8)
    ap.add_argument("--amp-tol", type=float, default=1e-12)
    ap.add_argument("--crosscheck-tol", type=float, default=1e-8)
    ap.add_argument("--leakage-tol", type=float, default=1e-9)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--joint-output", type=Path)
    ap.add_argument("--c-e", type=float)
    ap.add_argument("--c-l", type=float)
    ap.add_argument("--normalization-source")
    args = ap.parse_args()

    try:
        carrier = build_carrier(
            Jmax2=args.Jmax2,
            source_nodes=args.nodes,
            max_states=args.max_states,
            max_depth=args.max_depth,
            amp_tol=args.amp_tol,
            crosscheck_tol=args.crosscheck_tol,
            leakage_tol=args.leakage_tol,
        )
    except Exception as exc:
        carrier = {
            "status": "V11_BUILD_ERROR",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
            "Jmax2": args.Jmax2,
            "source_nodes": args.nodes,
            "physical_status": {
                "unique_hermitian_ordering_selected": False,
                "absolute_beta_kappa_hbar_regulator_prefactor_fixed": False,
                "physical_joint_kernel_verdict": False,
                "physical_causal_cone": False,
            },
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(carrier, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if args.joint_output:
        if carrier.get("status") != "PASS_STRUCTURAL_COMMON_CARRIER":
            raise ValueError("cannot create a joint-Gram input from a failed V11 carrier")
        if args.c_e is None or args.c_l is None or not args.normalization_source:
            raise ValueError("--joint-output requires --c-e, --c-l and --normalization-source")
        joint = joint_gram_input(carrier, args.c_e, args.c_l, args.normalization_source)
        args.joint_output.write_text(json.dumps(joint, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": carrier["status"],
        "Jmax2": args.Jmax2,
        "common_output_rows": carrier.get("common_output_basis", {}).get("size"),
        "crosschecks_passed": carrier.get("crosschecks", {}).get("all_passed"),
        "output": str(args.output),
        "physical_joint_kernel_verdict": False,
    }, indent=2))
    return 0 if carrier["status"] == "PASS_STRUCTURAL_COMMON_CARRIER" else 2


if __name__ == "__main__":
    raise SystemExit(main())
