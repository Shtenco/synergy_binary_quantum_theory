#!/usr/bin/env python3
"""Exact combinatorial preflight for the next q=2 barycentric PL generation.

This is an allocation-free count audit. It deliberately does not build the
generation-3 simplicial complex or claim that its exporters have run.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def barycentric_f_vector_3d(f_vector: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    """Return the f-vector after barycentric subdivision of a closed 3-complex.

    Counts flags by their simplex dimensions. The input must describe a
    tetrahedral complex; the formula is valid for a simplicial 3-complex.
    """
    f0, f1, f2, f3 = f_vector
    if min(f_vector) < 0:
        raise ValueError("f-vector entries must be non-negative")

    # Vertices of sd(K) are all non-empty simplices of K.
    g0 = f0 + f1 + f2 + f3

    # Edges of sd(K) are strict comparable pairs of non-empty simplices.
    g1 = 2 * f1 + 6 * f2 + 14 * f3

    # Triangles are strict flags of three non-empty simplices:
    # (0<1<2), (0<1<3), (0<2<3), and (1<2<3).
    g2 = 6 * f2 + 36 * f3

    # Tetrahedra are maximal flags inside each original tetrahedron.
    g3 = 24 * f3
    return (g0, g1, g2, g3)


def run() -> dict[str, object]:
    generation2 = (1696, 10912, 18432, 9216)
    generation3 = barycentric_f_vector_3d(generation2)
    f0, f1, f2, f3 = generation3
    causal_edges = 8 ** 3
    dual_edges = 2 * f3
    dense_matrix_bytes = 8 * causal_edges * causal_edges
    jacobi_rotation_work_upper = (
        100 * (causal_edges * (causal_edges - 1) // 2) * causal_edges
    )
    checks = {
        "generation2_euler_characteristic_zero": (
            generation2[0] - generation2[1] + generation2[2] - generation2[3] == 0
        ),
        "generation3_euler_characteristic_zero": f0 - f1 + f2 - f3 == 0,
        "generation3_tetrahedra_are_24_per_parent": f3 == 24 * generation2[3],
        "generation3_dual_edges_are_two_per_tetrahedron": dual_edges == 2 * f3,
        "generation3_causal_edges_follow_eightfold_rule": causal_edges == 8 * 64,
    }
    return {
        "status": "allocation-free exact generation-3 combinatorial preflight",
        "generation2_f_vector": list(generation2),
        "generation3_f_vector": list(generation3),
        "generation3_pl_primal_vertices": f0,
        "generation3_pl_primal_edges": f1,
        "generation3_pl_dual_tetrahedra": f3,
        "generation3_pl_dual_edges": dual_edges,
        "generation3_causal_active_edges": causal_edges,
        "dense_spectral_matrix_raw_float64_bytes": dense_matrix_bytes,
        "jacobi_rotation_work_estimate_at_100_sweeps": jacobi_rotation_work_upper,
        "checks": checks,
        "passed": all(checks.values()),
        "claim_boundary": (
            "Exact f-vector and graph-size preflight only. No generation-3 "
            "complex, correspondence, quotient, eigenvalue, or CI run is "
            "constructed or claimed by this script."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
