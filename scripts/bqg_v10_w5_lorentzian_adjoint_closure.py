#!/usr/bin/env python3
"""Build a W5-seeded, adjoint-closed finite basis for the raw Lorentzian KKV map.

The script computes BOTH the forward epsilon K-K-V map and its independently
ordered adjoint. It closes the reachable Gauss basis under both maps before
comparing the direct adjoint columns with the conjugate transpose of the
assembled forward matrix.

This is a structural adjoint-consistency gate only. It does not choose the
physical Hermitian ordering, restore beta/kappa/hbar factors, or claim a
physical Hamiltonian/kernel result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import deque
from pathlib import Path
import sys


W5_K = (0, 2, 4, 6, 8)
NODES = (0, 1, 2, 3, 4)
TOL = 1e-13


def epsilon_triples(neighbors):
    import itertools
    out = []
    for r, omitted in enumerate(neighbors):
        tri = tuple(x for x in neighbors if x != omitted)
        for p in itertools.permutations(tri):
            positions = [tri.index(x) for x in p]
            inversions = sum(
                positions[i] > positions[j]
                for i in range(3) for j in range(i + 1, 3)
            )
            sign = (-1) ** r * (-1 if inversions % 2 else 1)
            out.append((sign, p))
    if len(out) != 24 or len({p for _, p in out}) != 24:
        raise AssertionError("the K5 node epsilon orbit must contain 24 unique triples")
    return out


def w5_input_key(K):
    if K not in W5_K:
        raise ValueError(f"W5 K must be one of {W5_K}")
    return ((4,) * 10, (K, 0, 0, 0, 0))


def add_sparse(dst, src, scale=1.0, tol=TOL):
    for key, value in src.items():
        z = dst.get(key, 0j) + scale * value
        if abs(z) > tol:
            dst[key] = z
        elif key in dst:
            del dst[key]


def sparse_norm(state):
    return math.sqrt(sum(abs(z) ** 2 for z in state.values()))


def prune_scalar_channel(state):
    # The final rank-(0+1) covariant leg can reach total source J=0 only
    # from J=0 or J=1. Apply this exact selection rule after the middle leg.
    return {key: amp for key, amp in state.items() if key[2] in (0, 2)}


def triple_forward(psi_cov, source_v, a, b, c, Jmax2, cv_component, ck_component):
    """Tr_aux[C_a(K) C_b(K) C_c(V)] acting right-to-left."""
    out = {}
    max_leak = 0.0
    max_diag = {}
    for i in range(2):
        for j in range(2):
            for k in range(2):
                q, leak = cv_component(psi_cov, source_v, c, k, i, Jmax2)
                max_leak = max(max_leak, float(leak))
                if not q:
                    continue
                q, d1 = ck_component(q, source_v, b, j, k, Jmax2)
                _merge_diag(max_diag, d1)
                if not q:
                    continue
                q = prune_scalar_channel(q)
                if not q:
                    continue
                q, d2 = ck_component(q, source_v, a, i, j, Jmax2)
                _merge_diag(max_diag, d2)
                add_sparse(out, q)
    return out, max_leak, max_diag


def triple_adjoint(psi_cov, source_v, a, b, c, Jmax2, cv_component, ck_component):
    """Directly apply T_abc^dagger = Tr_aux[C_c(V) C_b(K) C_a(K)].

    K_sine is anti-Hermitian and V is Hermitian. The two K-adjoint minus
    signs cancel. The auxiliary indices are transposed explicitly:
      K_a[j,i] -> K_b[k,j] -> V_c[i,k]
    in right-to-left action order.
    """
    out = {}
    max_leak = 0.0
    max_diag = {}
    for i in range(2):
        for j in range(2):
            for k in range(2):
                q, d1 = ck_component(psi_cov, source_v, a, j, i, Jmax2)
                _merge_diag(max_diag, d1)
                if not q:
                    continue
                q, d2 = ck_component(q, source_v, b, k, j, Jmax2)
                _merge_diag(max_diag, d2)
                if not q:
                    continue
                q = prune_scalar_channel(q)
                if not q:
                    continue
                q, leak = cv_component(q, source_v, c, i, k, Jmax2)
                max_leak = max(max_leak, float(leak))
                add_sparse(out, q)
    return out, max_leak, max_diag


def node_kkv_maps(gauss_state, source_v, Jmax2, PW, CV, RAW):
    """Compute forward B_v and direct reversed-order B_v^dagger on one input."""
    psi_cov = CV.gauss_to_covariant(gauss_state, source_v)
    triples = epsilon_triples(PW.NEIG[source_v])
    B, Bdag = {}, {}
    diagnostics = {
        "CV_complete_basis_leakage": 0.0,
        "CK_outer_complete_basis_leakage": 0.0,
        "CK_internal_volume_sector_leakage": 0.0,
        "CK_complete_charge_basis_leakage": 0.0,
    }
    for epsilon, (a, b, c) in triples:
        forward, leak_f, diag_f = triple_forward(
            psi_cov, source_v, a, b, c, Jmax2,
            RAW.COMP.C_volume_component, RAW.KCOMP.C_K_component
        )
        adjoint, leak_a, diag_a = triple_adjoint(
            psi_cov, source_v, a, b, c, Jmax2,
            RAW.COMP.C_volume_component, RAW.KCOMP.C_K_component
        )
        add_sparse(B, forward, epsilon)
        add_sparse(Bdag, adjoint, epsilon)
        diagnostics["CV_complete_basis_leakage"] = max(
            diagnostics["CV_complete_basis_leakage"], leak_f, leak_a
        )
        _merge_diag(diagnostics, diag_f)
        _merge_diag(diagnostics, diag_a)
    diagnostics["CK_outer_complete_basis_leakage"] = max(
        float(diagnostics.get("outer_complete_basis_leakage", 0.0)),
        float(diagnostics.get("CK_outer_complete_basis_leakage", 0.0)),
    )
    diagnostics["CK_internal_volume_sector_leakage"] = max(
        float(diagnostics.get("internal_volume_sector_leakage", 0.0)),
        float(diagnostics.get("CK_internal_volume_sector_leakage", 0.0)),
    )
    diagnostics["CK_complete_charge_basis_leakage"] = max(
        float(diagnostics.get("complete_charge_basis_leakage", 0.0)),
        float(diagnostics.get("CK_complete_charge_basis_leakage", 0.0)),
    )
    return B, Bdag, diagnostics


def _merge_diag(dst, src):
    for name, value in src.items():
        if isinstance(value, (int, float)):
            dst[name] = max(dst.get(name, 0.0), float(value))


def covariant_scalar_to_gauss(key, source_v):
    spins, Kother, J2, M2, K12, K34 = key
    if J2 != 0 or M2 != 0 or K12 != K34:
        raise ValueError("epsilon-traced output is not a Gauss scalar")
    labels = list(Kother)
    if len(labels) != 5 or labels[source_v] != -1:
        raise ValueError("malformed covariant source sentinel")
    labels[source_v] = K12
    if any(x < 0 for x in labels):
        raise ValueError("unclosed non-source Gauss label")
    return (tuple(spins), tuple(labels))


def covariant_to_gauss_map(state, source_v):
    out, rejected, rejected_norm2 = {}, 0, 0.0
    for key, amp in state.items():
        try:
            gkey = covariant_scalar_to_gauss(key, source_v)
        except ValueError:
            rejected += 1
            rejected_norm2 += abs(amp) ** 2
            continue
        out[gkey] = out.get(gkey, 0j) + amp
    return out, rejected, math.sqrt(rejected_norm2)


def _matrix_adjoint_expected_columns(forward_columns, basis_keys, key_index):
    n = len(basis_keys)
    expected = [dict() for _ in range(n)]
    for col, state in enumerate(forward_columns):
        for row_key, amp in state.items():
            row = key_index.get(row_key)
            if row is None:
                continue
            expected[row][basis_keys[col]] = amp.conjugate()
    return expected


def sparse_matrix_residual(actual_columns, expected_columns):
    err2 = 0.0
    scale2 = 0.0
    for actual, expected in zip(actual_columns, expected_columns):
        for key in set(actual) | set(expected):
            a = actual.get(key, 0j)
            b = expected.get(key, 0j)
            err2 += abs(a - b) ** 2
            scale2 += abs(b) ** 2
    return math.sqrt(err2), math.sqrt(scale2)


def sparse_frobenius(columns):
    return math.sqrt(sum(abs(z) ** 2 for col in columns for z in col.values()))


def sparse_sum_difference_norm(columns_a, columns_b, sign=1.0):
    total2 = 0.0
    for a, b in zip(columns_a, columns_b):
        keys = set(a) | set(b)
        total2 += sum(abs(a.get(k, 0j) + sign * b.get(k, 0j)) ** 2 for k in keys)
    return math.sqrt(total2)


def run_closure(source_v=0, Jmax2=10, max_states=500, max_depth=8, amp_tol=1e-12):
    if source_v not in NODES:
        raise ValueError("source_v must be 0..4")
    if Jmax2 < 10:
        raise ValueError("W5 j=2 inputs require Jmax2>=10 as an initial support wall")
    if max_states < len(W5_K) or max_depth < 0 or amp_tol <= 0:
        raise ValueError("invalid closure limits")

    HERE = Path(__file__).resolve().parent
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    import k5_peter_weyl_safe_hda_column as PW
    import peter_weyl_covariant_volume_leg_gate as CV
    import peter_weyl_lorentzian_sine_ordered_triple_gate as LSINE
    import peter_weyl_lorentzian_ordered_triple_gate as RAW

    seeds = [w5_input_key(K) for K in W5_K]
    closure_keys = list(seeds)
    seen = set(closure_keys)
    depth = {key: 0 for key in seeds}
    queue = deque(seeds)
    forward_columns, adjoint_columns = [], []
    leakage_max = {}
    rejected_count = 0
    rejected_norm2 = 0.0
    truncated = False
    omitted_support = 0

    old_raw_run = RAW.run

    def closure_executor(source_v=0, a=None, b=None, c=None):
        nonlocal rejected_count, rejected_norm2, truncated, omitted_support
        while queue:
            key = queue.popleft()
            f_cov, a_cov, diag = node_kkv_maps(
                {key: 1.0 + 0j}, source_v, Jmax2, PW, CV, RAW
            )
            f, nf, nfnorm = covariant_to_gauss_map(f_cov, source_v)
            ad, na, nanorm = covariant_to_gauss_map(a_cov, source_v)
            rejected_count += nf + na
            rejected_norm2 += nfnorm**2 + nanorm**2
            for name, value in diag.items():
                leakage_max[name] = max(leakage_max.get(name, 0.0), float(value))

            # Store only amplitudes above the declared support tolerance.
            f = {k: z for k, z in f.items() if abs(z) > amp_tol}
            ad = {k: z for k, z in ad.items() if abs(z) > amp_tol}
            forward_columns.append(f)
            adjoint_columns.append(ad)

            for output in set(f) | set(ad):
                if output not in seen:
                    if len(closure_keys) >= max_states or depth[key] >= max_depth:
                        truncated = True
                        omitted_support += 1
                        continue
                    seen.add(output)
                    depth[output] = depth[key] + 1
                    closure_keys.append(output)
                    queue.append(output)

        # The intercepted result is adapted to the return schema expected by
        # LSINE.run. That wrapper supplies the physical sine K=[V,H_E] patch
        # and restores all monkey-patched functions in a finally block.
        norm_f = sparse_frobenius(forward_columns)
        return {
            "passed": bool(not truncated and rejected_count == 0),
            "output_support": sum(len(x) for x in forward_columns),
            "output_norm": norm_f,
            "final_nonscalar_J_weight_fraction": 0.0 if rejected_norm2 == 0 else 1.0,
            "max_spin_reached": max(
                (max(key[0]) for key in closure_keys), default=0
            ) / 2,
            "max_diagnostics": {
                "CV_complete_basis_leakage": leakage_max.get("CV_complete_basis_leakage", 0.0),
                "CK_outer_complete_basis_leakage": leakage_max.get("CK_outer_complete_basis_leakage", 0.0),
                "CK_internal_volume_sector_leakage": leakage_max.get("CK_internal_volume_sector_leakage", 0.0),
                "CK_complete_charge_basis_leakage": leakage_max.get("CK_complete_charge_basis_leakage", 0.0),
            },
        }

    RAW.run = closure_executor
    try:
        # LSINE.run installs the exact sine-Hermitian K=[V,H_E] engine and
        # cache wrappers once for the entire closure walk.
        LSINE.run(source_v)
    finally:
        RAW.run = old_raw_run

    complete = not truncated and not queue and len(forward_columns) == len(closure_keys)
    if not complete:
        return {
            "status": "INCOMPLETE_CLOSURE_LIMIT_REACHED",
            "source_node": source_v,
            "Jmax2": Jmax2,
            "Jmax": Jmax2 / 2,
            "seed_count": len(seeds),
            "closure_size_discovered": len(closure_keys),
            "columns_computed": len(forward_columns),
            "max_states": max_states,
            "max_depth": max_depth,
            "omitted_new_support_keys": omitted_support,
            "adjoint_residual": None,
            "claim_scope": "truncated reachable subspace only; no full adjoint-consistency verdict",
        }

    basis_keys = sorted(closure_keys, key=repr)
    key_index = {key: i for i, key in enumerate(basis_keys)}
    # Reorder computed columns from BFS order to the deterministic basis order.
    bfs_index = {key: i for i, key in enumerate(closure_keys)}
    forward_sorted = [forward_columns[bfs_index[key]] for key in basis_keys]
    adjoint_sorted = [adjoint_columns[bfs_index[key]] for key in basis_keys]

    # Verify every output is represented in the same closed basis.
    outside = [
        key for cols in (forward_sorted, adjoint_sorted)
        for col in cols for key in col if key not in key_index
    ]
    if outside:
        raise AssertionError("closure algorithm terminated with outputs outside basis")

    expected_adj = _matrix_adjoint_expected_columns(forward_sorted, basis_keys, key_index)
    err, scale = sparse_matrix_residual(adjoint_sorted, expected_adj)
    relative = err / max(scale, 1e-30)
    bnorm = sparse_frobenius(forward_sorted)
    bdag_norm = sparse_frobenius(adjoint_sorted)
    sym_norm = sparse_sum_difference_norm(forward_sorted, expected_adj, sign=+1.0)
    anti_norm = sparse_sum_difference_norm(forward_sorted, expected_adj, sign=-1.0)
    basis_payload = "\n".join(repr(k) for k in basis_keys).encode("utf-8")
    basis_id = hashlib.sha256(basis_payload).hexdigest()

    def encode_columns(cols):
        encoded = []
        for col in cols:
            encoded.append([
                [key_index[key], float(z.real), float(z.imag)]
                for key, z in sorted(col.items(), key=lambda kv: key_index[kv[0]])
            ])
        return encoded

    max_leak = max(leakage_max.values(), default=0.0)
    passed = (
        bnorm > 1e-10 and bdag_norm > 1e-10
        and relative < 1e-8 and rejected_count == 0 and max_leak < 1e-9
    )
    return {
        "status": "PASS_ADJOINT_CONSISTENCY" if passed else "FAIL_ADJOINT_OR_LEAKAGE_GATE",
        "source_node": source_v,
        "Jmax2": Jmax2,
        "Jmax": Jmax2 / 2,
        "seed_count": len(seeds),
        "closure_size": len(basis_keys),
        "closure_depth_max": max(depth.values(), default=0),
        "basis_id_sha256": basis_id,
        "basis_keys_repr": [repr(k) for k in basis_keys],
        "forward_B_columns_sparse": encode_columns(forward_sorted),
        "direct_Bdag_columns_sparse": encode_columns(adjoint_sorted),
        "diagnostics": {
            "direct_adjoint_relative_residual": relative,
            "direct_adjoint_absolute_residual": err,
            "forward_frobenius_norm": bnorm,
            "direct_adjoint_frobenius_norm": bdag_norm,
            "hermitian_part_frobenius_norm": 0.5 * sym_norm,
            "antihermitian_part_frobenius_norm": 0.5 * anti_norm,
            "scalar_projection_rejected_terms": rejected_count,
            "scalar_projection_rejected_norm": math.sqrt(rejected_norm2),
            "max_component_leakage": max_leak,
            "component_leakage_by_type": leakage_max,
        },
        "physical_status": {
            "physical_hermitian_ordering_selected": False,
            "beta_kappa_hbar_prefactor_fixed": False,
            "full_physical_HL": False,
            "physical_joint_kernel_verdict": False,
        },
        "claim_scope": "finite W5-seeded reachable component at declared Jmax2; raw structural KKV and independently ordered adjoint only",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--node", type=int, default=0)
    ap.add_argument("--Jmax2", type=int, default=10)
    ap.add_argument("--max-states", type=int, default=500)
    ap.add_argument("--max-depth", type=int, default=8)
    ap.add_argument("--amp-tol", type=float, default=1e-12)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    result = run_closure(args.node, args.Jmax2, args.max_states, args.max_depth, args.amp_tol)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "source_node": result["source_node"],
        "Jmax2": result["Jmax2"],
        "closure_size": result.get("closure_size", result.get("closure_size_discovered")),
        "direct_adjoint_relative_residual": result.get("diagnostics", {}).get("direct_adjoint_relative_residual"),
        "max_component_leakage": result.get("diagnostics", {}).get("max_component_leakage"),
        "output": str(args.output),
        "claim_scope": result["claim_scope"],
    }, indent=2))
    return 0 if result["status"] == "PASS_ADJOINT_CONSISTENCY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
