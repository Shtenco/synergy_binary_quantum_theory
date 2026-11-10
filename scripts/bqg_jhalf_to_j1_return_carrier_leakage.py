#!/usr/bin/env python3
"""Measure carrier leakage for the finite j=1/2 -> j=1 positive return operator.

This is deliberately narrower than the full graph-changing master. It takes
A = H_E,0^sine + H_E,1^sine and M=A^dagger A on the 243-dimensional product
of the five coarse local j=1 singlet multiplicity spaces (3^5). The exact
32-dimensional fine logical carrier is embedded by W_global.

The script computes M_c W_global by applying A twice, but accepts that route
only if its pullback agrees with the independently assembled Gram matrix
(AW)^dagger(AW) on the embedded carrier. This is an explicit finite
Hermiticity/cutoff consistency check, not an assumption hidden in the result.

The resulting L and D are therefore carrier-leakage/internal-mismatch blocks
for the compressed 243D return operator. They are not the full graph-changing
master residual and do not include leakage into spin-changed habitats.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import bqg_jhalf_to_j1_return_rg_gate as RG


def run():
    RG.ZVM.patch_and_clear()
    bits_all = RG.LOGICAL_BITS
    coarse_labels = tuple(itertools.product((0, 2, 4), repeat=5))
    coarse_index = {label: i for i, label in enumerate(coarse_labels)}
    Wg = np.zeros((len(coarse_labels), len(bits_all)), dtype=complex)

    fine_images = []
    coarse_images = []
    coarse_second_images = []
    rows = []

    coarse_spins = (2,) * len(RG.PW.EDGES)

    for q, bits in enumerate(bits_all):
        sf = RG.fine_logical_state(bits)
        sc = RG.coarse_logical_state(bits)
        for (spins, Ks), amp in sc.items():
            if spins != coarse_spins or Ks not in coarse_index:
                raise RuntimeError(f"embedded coarse vector left declared 3^5 carrier: {spins}, {Ks}")
            Wg[coarse_index[Ks], q] = amp

        af = RG.apply_H01(sf, RG.JMAX2_FINE)
        ac = RG.apply_H01(sc, RG.JMAX2_COARSE)
        a2c = RG.apply_H01(ac, RG.JMAX2_COARSE)
        fine_images.append(af)
        coarse_images.append(ac)
        coarse_second_images.append(a2c)
        rows.append({
            "logical_bits": list(bits),
            "fine_image_support": len(af),
            "coarse_image_support": len(ac),
            "coarse_second_image_support": len(a2c),
            "coarse_base_sector_support": sum(
                1 for (spins, Ks) in a2c
                if spins == coarse_spins and Ks in coarse_index
            ),
        })

    isometry_error = float(np.linalg.norm(Wg.conj().T @ Wg - np.eye(len(bits_all))))
    if isometry_error > 1e-11:
        raise RuntimeError(f"global refinement map is not isometric: {isometry_error}")

    Kf = RG.gram(fine_images)
    Kc_pullback_gram = RG.gram(coarse_images)

    # Coordinates of P_coarse A^2 W in the orthonormal 3^5 coarse singlet basis.
    # A^2 is accepted as A^dagger A only to the extent independently validated
    # by W^dagger A^2 W == (AW)^dagger(AW).
    McW = np.zeros((len(coarse_labels), len(bits_all)), dtype=complex)
    for q, image in enumerate(coarse_second_images):
        for (spins, Ks), amp in image.items():
            if spins == coarse_spins and Ks in coarse_index:
                McW[coarse_index[Ks], q] += amp

    pullback_from_second_action = Wg.conj().T @ McW
    hermiticity_pullback_defect = float(
        np.linalg.norm(pullback_from_second_action - Kc_pullback_gram, ord="fro")
        / max(np.linalg.norm(Kc_pullback_gram, ord="fro"), 1e-30)
    )

    P = Wg @ Wg.conj().T
    L = (np.eye(len(coarse_labels)) - P) @ McW
    D = pullback_from_second_action - Kf
    R = McW - Wg @ Kf

    split_defect = float(np.linalg.norm(
        R.conj().T @ R - (L.conj().T @ L + D.conj().T @ D), ord="fro"
    ))
    orthogonality_defect = float(np.linalg.norm(Wg.conj().T @ L, ord="fro"))

    def norms(x):
        return {
            "spectral": float(np.linalg.norm(x, ord=2)),
            "frobenius": float(np.linalg.norm(x, ord="fro")),
        }

    result = {
        "status": "finite compressed-return carrier leakage diagnostic",
        "fine_logical_dimension": len(bits_all),
        "coarse_logical_dimension": len(coarse_labels),
        "operator": "A^dagger A, with A=H_E,0^sine+H_E,1^sine; A^2 route cross-checked against Gram pullback",
        "fine_Jmax2": RG.JMAX2_FINE,
        "coarse_Jmax2": RG.JMAX2_COARSE,
        "global_intertwiner": "tensor product of exact local S4-[2,2] W",
        "global_W_isometry_error": isometry_error,
        "hermiticity_pullback_relative_defect": hermiticity_pullback_defect,
        "Kf_spectrum": RG.spectrum_summary(Kf),
        "Kc_pullback_spectrum": RG.spectrum_summary(Kc_pullback_gram),
        "residual_decomposition": {
            "L_carrier_leakage": norms(L),
            "D_internal_mismatch": norms(D),
            "R_total_compressed_residual": norms(R),
            "orthogonality_defect_WdaggerL": orthogonality_defect,
            "identity_defect_RdaggerR_minus_LdaggerL_minus_DdaggerD": split_defect,
            "relative_L2_pullback_vs_Gram": hermiticity_pullback_defect,
        },
        "claim_boundary": (
            "L and D refer only to the 243D product of coarse local j=1 singlet "
            "multiplicity spaces and the positive two-node Euclidean return operator. "
            "This does not include spin-changed/off-carrier graph-changing habitats, "
            "the full constraint family, the full master gap, or continuum scaling. "
            "No scale-free physical-projector epsilon is reported without the full "
            "coarse master gap."
        ),
        "columns": rows,
    }
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
