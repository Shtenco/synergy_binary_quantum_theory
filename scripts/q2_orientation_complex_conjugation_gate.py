#!/usr/bin/env python3
"""Exact q=2 orientation <-> complex-conjugation compatibility gate.

The same reflection of the four q=2 labels is tested in two representations:
1. oriented C4 phase carrier: R J R^{-1} = -J;
2. Walsh tetrahedral carrier: the oriented determinant changes sign.

Because the existing geometry-qubit theorem has Q = sqrt(3)/4 Y_L with Q an
oriented triple product, the same odd face-label reflection also sends
Y_L -> -Y_L at the representation-theory level.

The gate proves common covariance under one label reflection. It does NOT claim
that the microscopic Hamiltonian has already identified the phase and geometry
representations dynamically.
"""

import argparse
import json


def mmul(A, B):
    return tuple(
        tuple(sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0])))
        for i in range(len(A))
    )


def transpose(A):
    return tuple(tuple(A[j][i] for j in range(len(A))) for i in range(len(A[0])))


def matvec(A, v):
    return tuple(sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A)))


def det3(M):
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, i = M[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def sub(u, v):
    return tuple(a - b for a, b in zip(u, v))


def walsh(label):
    x, y = label
    # Three nontrivial characters chi_01, chi_10, chi_11.
    return (
        -1 if y else 1,
        -1 if x else 1,
        -1 if (x ^ y) else 1,
    )


def oriented_tetra_det(points):
    p0, p1, p2, p3 = points
    c1 = sub(p1, p0)
    c2 = sub(p2, p0)
    c3 = sub(p3, p0)
    # columns c1,c2,c3
    M = (
        (c1[0], c2[0], c3[0]),
        (c1[1], c2[1], c3[1]),
        (c1[2], c2[2], c3[2]),
    )
    return det3(M)


def permutation_parity(perm):
    inv = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                inv += 1
    return -1 if inv % 2 else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", dest="json_path")
    args = ap.parse_args()

    # Gray-cycle ordering of the q=2 labels.
    labels = [(0, 0), (0, 1), (1, 1), (1, 0)]

    # S|k> = |k+1 mod4>.
    S = (
        (0, 0, 0, 1),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, 0, 1, 0),
    )

    # R:k -> -k mod4. It fixes k=0,2 and swaps k=1,3.
    R = (
        (1, 0, 0, 0),
        (0, 0, 0, 1),
        (0, 0, 1, 0),
        (0, 1, 0, 0),
    )
    Rinv = transpose(R)
    Sinv = transpose(S)

    dihedral_relation = mmul(mmul(R, S), Rinv) == Sinv

    # Real Fourier/phase plane.
    uc = (1, 0, -1, 0)
    us = (0, 1, 0, -1)
    phase_checks = {
        "S_uc_eq_us": matvec(S, uc) == us,
        "S_us_eq_minus_uc": matvec(S, us) == tuple(-x for x in uc),
        "R_uc_eq_uc": matvec(R, uc) == uc,
        "R_us_eq_minus_us": matvec(R, us) == tuple(-x for x in us),
    }
    # In basis (uc,us), J=[[0,-1],[1,0]], R2=diag(1,-1).
    J = ((0, -1), (1, 0))
    R2 = ((1, 0), (0, -1))
    minusJ = tuple(tuple(-x for x in row) for row in J)
    phase_checks["RJR_eq_minus_J"] = mmul(mmul(R2, J), R2) == minusJ
    phase_checks["pass"] = dihedral_relation and all(phase_checks.values())

    walsh_points = [walsh(g) for g in labels]
    det_before = oriented_tetra_det(walsh_points)

    # Same R as a permutation of label positions: [0,3,2,1].
    perm = [0, 3, 2, 1]
    reflected_points = [walsh_points[i] for i in perm]
    det_after = oriented_tetra_det(reflected_points)
    parity = permutation_parity(perm)
    geometry_checks = {
        "det_before": det_before,
        "det_after": det_after,
        "reflection_parity": parity,
        "det_sign_flips": det_after == -det_before and det_before != 0,
        "odd_face_permutation": parity == -1,
    }
    geometry_checks["pass"] = geometry_checks["det_sign_flips"] and geometry_checks["odd_face_permutation"]

    # The logical geometry result Q=(sqrt(3)/4)Y_L is already an exact theorem
    # in SPATIAL_QUBIT_GEOMETRY_BRIDGE.md. For an oriented triple product,
    # an odd permutation flips Q, hence flips Y_L. We record the sign algebra
    # rather than importing floating/SU(2) matrices here.
    common_sign = {
        "phase_complex_structure_sign_under_R": -1 if phase_checks["RJR_eq_minus_J"] else None,
        "walsh_tetra_orientation_sign_under_R": det_after // det_before if det_before and det_after % det_before == 0 else None,
        "logical_Y_pseudoscalar_expected_sign_for_odd_face_permutation": parity,
    }
    common_sign["all_signs_equal_minus_one"] = set(common_sign.values()) == {-1}
    common_sign["pass"] = common_sign["all_signs_equal_minus_one"]

    report = {
        "labels_gray_cycle": ["00", "01", "11", "10"],
        "reflection_index_map": perm,
        "dihedral_R_S_R_equals_S_inverse": dihedral_relation,
        "phase_representation": phase_checks,
        "walsh_tetrahedron": {
            "points": walsh_points,
            **geometry_checks,
        },
        "common_orientation_sign": common_sign,
    }
    report["overall_pass"] = phase_checks["pass"] and geometry_checks["pass"] and common_sign["pass"]

    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.json_path:
        with open(args.json_path, "w", encoding="utf-8") as f:
            f.write(text + "\n")
    if not report["overall_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
