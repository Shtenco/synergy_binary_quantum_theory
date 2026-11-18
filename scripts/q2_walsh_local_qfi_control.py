#!/usr/bin/env python3
"""QFI control derived from the canonical q=2 Walsh tetrahedral face qubits.

This is a local information-geometric control, not a physical edge-conductance
assignment. It computes:
  * the 3x3 pure-qubit QFI for SU(2) frame rotations of each Walsh face state;
  * the sum over four faces;
  * the QFI of the Gauss-projected four-face singlet under a global SU(2)
    gauge rotation (which must vanish);
  * a reduced-model relative z-twist probe K=(Z_i-Z_j)/2 on product face
    states, for which F_Q/4=1/3 for every pair.

The result identifies the missing ingredient: a global gauge rotation cannot
supply positive network conductance after Gauss projection. A relational
inter-block generator and its physical reduced state are still required.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import micro_walsh_qgeom_gate as MW
import spatial_qubit_geometry_gate as SQ


def pure_state_qfi(generators, psi):
    means = np.array([np.vdot(psi, G @ psi) for G in generators], complex)
    n = len(generators)
    out = np.zeros((n, n), dtype=float)
    for a, A in enumerate(generators):
        for b, B in enumerate(generators):
            sym = 0.5 * np.vdot(psi, (A @ B + B @ A) @ psi)
            out[a, b] = 4.0 * float((sym - means[a] * means[b]).real)
    return 0.5 * (out + out.T)


def run():
    labels, normals = MW.route_flux_frame()
    I3 = np.eye(3)
    single_qfi = [I3 - np.outer(n, n) for n in normals]
    qfi_sum = np.sum(single_qfi, axis=0)
    qfi_mean = qfi_sum / len(normals)

    spinors = [MW.spinor_from_bloch(n) for n in normals]
    product = MW.kron_all(spinors)
    logical = SQ.logical_basis()
    p_singlet = logical @ logical.conj().T
    projected = p_singlet @ product
    singlet_weight = float(np.vdot(projected, projected).real)
    if singlet_weight <= 1e-14:
        raise RuntimeError("Walsh product has no Gauss-singlet support")
    singlet = projected / np.sqrt(singlet_weight)

    total_generators = []
    for axis, pauli in enumerate(SQ.PAULI):
        total_generators.append(sum(
            (0.5 * SQ.one_site(pauli, site) for site in range(4)),
            np.zeros((16, 16), dtype=complex),
        ))
    global_singlet_qfi = pure_state_qfi(total_generators, singlet)
    global_generator_residual = max(
        float(np.linalg.norm(G @ singlet)) for G in total_generators
    )

    relative_probe_rows = []
    z_generators = [0.5 * SQ.one_site(SQ.Z, site) for site in range(4)]
    for i in range(4):
        for j in range(i + 1, 4):
            K = z_generators[i] - z_generators[j]
            mean = np.vdot(product, K @ product)
            second = np.vdot(product, K @ K @ product)
            fq = 4.0 * float((second - mean * mean).real)
            relative_probe_rows.append({
                "face_pair": [i, j],
                "F_Q": fq,
                "reduced_probe_conductance_FQ_over_4": fq / 4.0,
            })

    checks = {
        "four_walsh_face_qubits": len(labels) == 4 and normals.shape == (4, 3),
        "single_face_qfi_is_positive_semidefinite": all(
            float(np.linalg.eigvalsh(F).min()) > -1e-12 for F in single_qfi
        ),
        "single_face_qfi_matches_tetrahedral_formula": all(
            np.linalg.norm(F - (I3 - np.outer(normals[i], normals[i]))) < 1e-12
            for i, F in enumerate(single_qfi)
        ),
        "four_face_qfi_sum_isotropic_8_over_3": np.linalg.norm(
            qfi_sum - (8.0 / 3.0) * I3
        ) < 1e-12,
        "gauss_singlet_weight_is_2_over_9": abs(singlet_weight - 2.0 / 9.0) < 1e-12,
        "global_su2_generator_annihilates_singlet": global_generator_residual < 1e-12,
        "global_gauge_qfi_vanishes_after_projection": np.linalg.norm(global_singlet_qfi) < 1e-12,
        "all_relative_product_probe_conductances_equal_one_third": all(
            abs(row["reduced_probe_conductance_FQ_over_4"] - 1.0 / 3.0) < 1e-12
            for row in relative_probe_rows
        ),
    }
    return {
        "status": "q2 Walsh face-qubit local QFI control",
        "passed": bool(all(checks.values())),
        "route_labels": [list(x) for x in labels],
        "tetrahedral_normals": normals.tolist(),
        "single_face_QFI_matrices": [F.tolist() for F in single_qfi],
        "four_face_QFI_sum": qfi_sum.tolist(),
        "four_face_QFI_mean": qfi_mean.tolist(),
        "gauss_singlet_weight": singlet_weight,
        "global_SU2_singlet_generator_residual": global_generator_residual,
        "global_SU2_QFI_after_Gauss_projection": global_singlet_qfi.tolist(),
        "relative_product_probe": relative_probe_rows,
        "checks": checks,
        "interpretation": (
            "The unprojected tetrahedral face frame has an isotropic summed local "
            "QFI tensor 8/3 I. The Gauss-projected state is invariant under global "
            "SU(2), so its global gauge-orbit QFI is exactly zero. The reduced-model "
            "relative z-twist gives FQ/4=1/3 for each product-state face pair, a "
            "uniform scalar that only rescales a chosen graph Laplacian."
        ),
        "claim_boundary": (
            "No physical inter-block edge state or relational generator is derived. "
            "The value 1/3 is a local product-state control, not a BQG QFI conductance. "
            "It cannot be inserted as the physical weight of the PL dual graph."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
