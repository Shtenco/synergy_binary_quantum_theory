#!/usr/bin/env python3
"""Common Gauss-key seed-column audit for Euclidean and Lorentzian BQG maps.

This is a bounded compatibility experiment, NOT a joint physical projector.
It applies the existing graph-changing Euclidean H_E and the existing raw
Lorentzian KKV map plus its independently ordered adjoint to the same W5
Gauss-key seed states. It records support overlap without assuming that the
two maps are already one anomaly-free, physical constraint family.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))


def norm2(state: dict) -> float:
    return float(sum(abs(v) ** 2 for v in state.values()))


def overlap(a: dict, b: dict) -> complex:
    return sum(a.get(k, 0j).conjugate() * b.get(k, 0j)
               for k in set(a) | set(b))


def zero_spin_fraction(state: dict) -> float:
    total = norm2(state)
    if total <= 1e-30:
        return 0.0
    changed = sum(abs(a) ** 2 for (spins, _), a in state.items()
                  if any(int(j2) == 0 for j2 in spins))
    return float(changed / total)


def valid_gauss_key(key) -> bool:
    if not isinstance(key, tuple) or len(key) != 2:
        return False
    spins, labels = key
    return (
        isinstance(spins, tuple) and len(spins) == 10
        and isinstance(labels, tuple) and len(labels) == 5
        and all(isinstance(x, int) and x >= 0 for x in spins)
        and all(isinstance(x, int) and x >= 0 for x in labels)
    )


def run(seed_K=(0, 2, 4, 6, 8), source_node=0, Jmax2=10):
    import k5_peter_weyl_safe_hda_column as PW
    import peter_weyl_covariant_volume_leg_gate as CV
    import peter_weyl_lorentzian_ordered_triple_gate as RAW
    import bqg_v10_w5_lorentzian_adjoint_closure as V10
    import peter_weyl_zeroaware_volume_migration_experiment as ZVM

    if Jmax2 < 10:
        raise ValueError("W5 j=2 seeds require Jmax2 >= 10")
    if source_node not in V10.NODES:
        raise ValueError("source_node must be in 0..4")
    if not seed_K or any(k not in V10.W5_K for k in seed_K):
        raise ValueError(f"seed_K must be a nonempty subset of {V10.W5_K}")

    # Match the exact zero-aware convention used by the frozen Lorentzian
    # closure runner before applying either Lorentzian map.
    ZVM.patch_and_clear()

    rows = []
    all_euclidean_keys = set()
    all_lorentzian_keys = set()
    all_adjoint_keys = set()
    rejected_norm2_total = 0.0
    max_leakage = 0.0

    for K in seed_K:
        key = V10.w5_input_key(K)
        seed = {key: 1.0 + 0j}

        # Existing graph-changing Euclidean operator, not a surrogate matrix.
        E = PW.prune_state(PW.apply_H_cached_state(seed, source_node, Jmax2), 1e-11)

        # Existing Lorentzian KKV map and its independently ordered adjoint.
        Lcov, Ladjcov, diag = V10.node_kkv_maps(
            seed, source_node, Jmax2, PW, CV, RAW
        )
        L, rejected_L, rejected_norm_L = V10.covariant_to_gauss_map(
            Lcov, source_node
        )
        Ladj, rejected_A, rejected_norm_A = V10.covariant_to_gauss_map(
            Ladjcov, source_node
        )
        L = {k: v for k, v in L.items() if abs(v) > 1e-11}
        Ladj = {k: v for k, v in Ladj.items() if abs(v) > 1e-11}

        ekeys, lkeys, akeys = set(E), set(L), set(Ladj)
        all_euclidean_keys.update(ekeys)
        all_lorentzian_keys.update(lkeys)
        all_adjoint_keys.update(akeys)
        rejected_norm2_total += rejected_norm_L**2 + rejected_norm_A**2
        max_leakage = max(max_leakage, *[
            float(v) for v in diag.values() if isinstance(v, (int, float))
        ])

        ne, nl, na = math.sqrt(norm2(E)), math.sqrt(norm2(L)), math.sqrt(norm2(Ladj))
        ov = overlap(E, L)
        cosine = abs(ov) / max(ne * nl, 1e-30)
        common = ekeys & lkeys
        rows.append({
            "seed_K": K,
            "seed_key_valid": valid_gauss_key(key),
            "euclidean_support": len(E),
            "lorentzian_forward_support": len(L),
            "lorentzian_direct_adjoint_support": len(Ladj),
            "euclidean_norm": ne,
            "lorentzian_forward_norm": nl,
            "lorentzian_direct_adjoint_norm": na,
            "euclidean_lorentzian_common_keys": len(common),
            "euclidean_lorentzian_support_union": len(ekeys | lkeys),
            "euclidean_lorentzian_normalized_overlap": float(cosine),
            "euclidean_j0_graph_change_fraction": zero_spin_fraction(E),
            "lorentzian_j0_graph_change_fraction": zero_spin_fraction(L),
            "direct_adjoint_j0_graph_change_fraction": zero_spin_fraction(Ladj),
            "lorentzian_scalar_rejections": rejected_L,
            "adjoint_scalar_rejections": rejected_A,
            "lorentzian_scalar_rejected_norm": rejected_norm_L,
            "adjoint_scalar_rejected_norm": rejected_norm_A,
            "max_operator_diagnostic": diag,
            "all_outputs_use_gauss_key_schema": all(
                valid_gauss_key(k) for k in ekeys | lkeys | akeys
            ),
        })

    nonzero = all(
        r["euclidean_norm"] > 1e-12
        and r["lorentzian_forward_norm"] > 1e-12
        and r["lorentzian_direct_adjoint_norm"] > 1e-12
        for r in rows
    )
    schema = all(r["seed_key_valid"] and r["all_outputs_use_gauss_key_schema"]
                 for r in rows)
    scalar_clean = all(
        r["lorentzian_scalar_rejections"] == 0
        and r["adjoint_scalar_rejections"] == 0
        and r["lorentzian_scalar_rejected_norm"] < 1e-10
        and r["adjoint_scalar_rejected_norm"] < 1e-10
        for r in rows
    )
    return {
        "status": "bounded common Gauss-key seed-column compatibility audit",
        "passed": bool(nonzero and schema and scalar_clean),
        "source_node": source_node,
        "Jmax2": Jmax2,
        "seed_K": list(seed_K),
        "seed_rows": rows,
        "union_supports": {
            "euclidean": len(all_euclidean_keys),
            "lorentzian_forward": len(all_lorentzian_keys),
            "lorentzian_direct_adjoint": len(all_adjoint_keys),
            "euclidean_and_lorentzian_forward_intersection": len(
                all_euclidean_keys & all_lorentzian_keys
            ),
            "joint_union": len(all_euclidean_keys | all_lorentzian_keys | all_adjoint_keys),
        },
        "max_lorentzian_component_leakage_diagnostic": max_leakage,
        "total_rejected_scalar_norm_squared": rejected_norm2_total,
        "checks": {
            "all_seed_and_output_keys_share_gauss_schema": bool(schema),
            "all_three_actual_maps_nonzero_on_each_seed": bool(nonzero),
            "lorentzian_forward_and_direct_adjoint_scalar_reductions_clean": bool(scalar_clean),
        },
        "claim_boundary": (
            "PASS means only that these existing finite maps act on a common "
            "Gauss-key seed domain and their scalar outputs can be represented "
            "in the same kinematical key schema. It does not establish a common "
            "self-adjoint constraint family, closure of the combined operator "
            "domain, a joint kernel, a physical inner product, or a rigging-map limit."
        ),
        "physical_gates_still_open": [
            "common graph-changing constraint family and anomaly-controlled domains",
            "joint Euclidean/Lorentzian constraint operator and its physical ordering",
            "closed/refinement-stable common operator domain",
            "nontrivial common physical kernel or rigging map",
            "normalized Lorentzian history and relational metric",
        ],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-node", type=int, default=0)
    ap.add_argument("--Jmax2", type=int, default=10)
    ap.add_argument("--seed-K", type=int, nargs="+", default=[0, 2, 4, 6, 8])
    ap.add_argument("--output", type=Path, default=Path(
        "verification_results/BQG_COMMON_GAUSS_HABITAT_AUDIT.json"
    ))
    args = ap.parse_args()
    result = run(tuple(args.seed_K), args.source_node, args.Jmax2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "passed": result["passed"],
        "source_node": result["source_node"],
        "Jmax2": result["Jmax2"],
        "union_supports": result["union_supports"],
        "checks": result["checks"],
        "claim_boundary": result["claim_boundary"],
    }, indent=2))
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
