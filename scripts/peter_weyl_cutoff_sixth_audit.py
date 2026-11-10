#!/usr/bin/env python3
"""Sixth Peter-Weyl cutoff audit with an independent exact Haar-moment control.

The CG/sparse-state evolution is checked at Jmax=3 and against an over-cutoff
Jmax=7/2 run.  Independently, the exact control evaluates the SU(2)^3 Haar
integral using the unit-quaternion representation and a rational characteristic-
function recurrence; it does not use the Peter-Weyl CG transition matrices.
"""
from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

import peter_weyl_cutoff_staircase as PW


def _convolve(a, b, limit):
    out = [Fraction(0) for _ in range(limit + 1)]
    for i, x in enumerate(a):
        if not x:
            continue
        for j, y in enumerate(b):
            if i + j > limit:
                break
            if y:
                out[i + j] += x * y
    return out


def exact_quaternion_haar_control(n: int = 6) -> int:
    """Return E[W^(2n)] exactly for W=sum_{i<j} chi(q_i q_j^-1).

    Fix q_4=1 by global SU(2) invariance.  The other three q_i are independent
    uniform unit quaternions in R^4.  Since chi(q_i q_j^-1)=2 q_i.q_j,
    W=2 sum_{i<j} q_i.q_j=|q_1+q_2+q_3+q_4|^2-4.

    The characteristic function of one uniform S^3 vector is 2 J_1(r)/r.
    In four dimensions:
      phi(t)=exp(i t_0) [2 J_1(|t|)/|t|]^3.
    Rational coefficient extraction plus repeated Laplacians gives all
    moments E[|sum q_i|^(2k)] exactly.
    """
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("n must be a positive integer")
    degree = 2 * n
    # b_m are the exact coefficients of 2 J_1(sqrt(u))/sqrt(u).
    b = [
        Fraction((-1) ** m, 4**m * math.factorial(m) * math.factorial(m + 1))
        for m in range(degree + 1)
    ]
    b2 = _convolve(b, b, degree)
    b3 = _convolve(b2, b, degree)

    @lru_cache(None)
    def laplacian_at_origin(k: int, m: int, p: int) -> int:
        # Delta^k (|t|^(2m) t_0^p) at t=0 in R^4.
        if k == 0:
            return 1 if (m, p) == (0, 0) else 0
        value = 0
        if p >= 2:
            value += p * (p - 1) * laplacian_at_origin(k - 1, m, p - 2)
        if m >= 1:
            value += 2 * m * (2 * m + 2 + 2 * p) * laplacian_at_origin(k - 1, m - 1, p)
        return value

    radial_moments = []
    for k in range(degree + 1):
        moment = Fraction(0)
        for m in range(k + 1):
            p = 2 * k - 2 * m
            moment += b3[m] * Fraction((-1) ** m * laplacian_at_origin(k, m, p),
                                       math.factorial(p))
        radial_moments.append(moment)

    result = sum(
        Fraction(math.comb(degree, k) * (-4) ** (degree - k)) * radial_moments[k]
        for k in range(degree + 1)
    )
    if result.denominator != 1:
        raise ArithmeticError(f"exact Haar moment is not integral: {result}")
    return int(result)


def gauge_covariance_error(jmax: float) -> float:
    basis, idx, L, R, U, Ud = PW.build_link(jmax)
    _, Jfund = PW.spin_mats(0.5)
    error = 0.0
    for a in range(3):
        for i, j in __import__("itertools").product(range(2), repeat=2):
            left = L[a] @ U[(i, j)] - U[(i, j)] @ L[a]
            left -= sum(Jfund[a][p, i] * U[(p, j)] for p in range(2))
            right = R[a] @ U[(i, j)] - U[(i, j)] @ R[a]
            right -= sum(U[(i, p)] * Jfund[a][p, j] for p in range(2))
            error = max(error, float(np.linalg.norm(left)), float(np.linalg.norm(right)))
    return error


def sparse_transitions_without_hidden_threshold(matrix):
    """Keep every floating nonzero; all state pruning is explicit and logged."""
    transitions = {}
    tiny_count = 0
    minimum_nonzero = math.inf
    for column in range(matrix.shape[1]):
        rows = np.flatnonzero(matrix[:, column] != 0.0)
        entries = []
        for row in rows:
            value = matrix[row, column]
            magnitude = abs(value)
            minimum_nonzero = min(minimum_nonzero, magnitude)
            if magnitude <= 1e-14:
                tiny_count += 1
            entries.append((int(row), value))
        transitions[column] = entries
    return transitions, {"nonzero_transition_count": sum(len(x) for x in transitions.values()),
                         "transition_entries_at_or_below_1e-14": tiny_count,
                         "minimum_nonzero_transition_magnitude": (
                             float(minimum_nonzero) if math.isfinite(minimum_nonzero) else None
                         )}


def evolve_detailed(jmax: float, nmax: int, prune_tol: float):
    basis, idx, _, _, U, Ud = PW.build_link(jmax)
    D = len(basis)
    Ut, transition_diag_u = {}, {}
    Dt, transition_diag_ud = {}, {}
    for key, matrix in U.items():
        Ut[key], transition_diag_u[key] = sparse_transitions_without_hidden_threshold(matrix)
    for key, matrix in Ud.items():
        Dt[key], transition_diag_ud[key] = sparse_transitions_without_hidden_threshold(matrix)
    transition_diagnostics = {
        "fundamental": {f"{key[0]},{key[1]}": value for key, value in transition_diag_u.items()},
        "adjoint": {f"{key[0]},{key[1]}": value for key, value in transition_diag_ud.items()},
        "hidden_transition_threshold": 0.0,
    }
    zero = idx[(0.0, 0.0, 0.0)]
    state = {(zero, zero, zero, zero): 1.0 + 0.0j}
    pairs = list(__import__("itertools").combinations(range(4), 2))
    rows = []
    for power in range(1, nmax + 1):
        out = {}
        for ds, amp in state.items():
            for e1, e2 in pairs:
                for a, b in __import__("itertools").product(range(2), repeat=2):
                    for o1, c1 in Ut[(a, b)][ds[e1]]:
                        for o2, c2 in Dt[(a, b)][ds[e2]]:
                            nd = list(ds)
                            nd[e1] = o1
                            nd[e2] = o2
                            nd = tuple(nd)
                            out[nd] = out.get(nd, 0.0j) + amp * c1 * c2
        finite = all(np.isfinite([z.real, z.imag]).all() for z in out.values())
        if not finite:
            raise FloatingPointError(f"non-finite amplitude at power {power}, Jmax={jmax}")
        exact_nonzero = {key: value for key, value in out.items() if value != 0.0}
        dropped = [value for value in exact_nonzero.values() if abs(value) <= prune_tol]
        state = {key: value for key, value in exact_nonzero.items() if abs(value) > prune_tol}
        norm2 = float(sum(abs(value) ** 2 for value in state.values()))
        max_spin = max((max(basis[d][0] for d in key) for key in state), default=0.0)
        rows.append({
            "power": power,
            "support_after_pruning": len(state),
            "support_before_pruning_float_nonzero": len(exact_nonzero),
            "pruned_amplitude_count": len(dropped),
            "pruned_amplitude_l2_norm": float(math.sqrt(sum(abs(x) ** 2 for x in dropped))),
            "norm_squared": norm2,
            "max_spin_reached": float(max_spin),
            "all_amplitudes_finite": True,
        })
    return {
        "link_dimension": D,
        "transition_diagnostics": transition_diagnostics,
        "rows": rows,
        "final_state": state,
        "final_support": len(state),
        "final_norm_squared": rows[-1]["norm_squared"],
        "final_max_spin": rows[-1]["max_spin_reached"],
    }


def compare_states(left, right):
    keys = set(left) | set(right)
    differences = [left.get(key, 0.0j) - right.get(key, 0.0j) for key in keys]
    return {
        "support_union": len(keys),
        "left_only_support": len(set(left) - set(right)),
        "right_only_support": len(set(right) - set(left)),
        "max_amplitude_difference": max((abs(x) for x in differences), default=0.0),
        "state_l2_difference": float(math.sqrt(sum(abs(x) ** 2 for x in differences))),
    }


def run():
    n = 6
    exact = exact_quaternion_haar_control(n)
    known = [6, 120, 4650, 254604, 16807392]
    exact_prefix = [exact_quaternion_haar_control(k) for k in range(1, 6)]
    exact_values = exact_prefix + [exact]
    if exact_prefix != known:
        raise RuntimeError(f"independent exact control failed known prefix: {exact_prefix}")

    baseline = evolve_detailed(3.0, n, 1e-13)
    tight = evolve_detailed(3.0, n, 1e-14)
    loose = evolve_detailed(3.0, n, 1e-12)
    overcut = evolve_detailed(3.5, n, 1e-13)

    cutoff_comparison = compare_states(baseline["final_state"], overcut["final_state"])
    tight_comparison = compare_states(baseline["final_state"], tight["final_state"])
    loose_comparison = compare_states(baseline["final_state"], loose["final_state"])
    exact_error = abs(baseline["final_norm_squared"] - exact)
    norm_history_errors = {
        "Jmax_3_prune_1e-13": [
            abs(row["norm_squared"] - exact_values[row["power"] - 1])
            for row in baseline["rows"]
        ],
        "Jmax_3_prune_1e-14": [
            abs(row["norm_squared"] - exact_values[row["power"] - 1])
            for row in tight["rows"]
        ],
        "Jmax_3_prune_1e-12": [
            abs(row["norm_squared"] - exact_values[row["power"] - 1])
            for row in loose["rows"]
        ],
        "Jmax_7_2_prune_1e-13": [
            abs(row["norm_squared"] - exact_values[row["power"] - 1])
            for row in overcut["rows"]
        ],
    }
    gauge_covariance = {
        "Jmax_3": gauge_covariance_error(3.0),
        "Jmax_7_2": gauge_covariance_error(3.5),
    }
    # Link basis is identical through j=3 for both cutoffs; indices above the
    # J=3 sector start at the common prefix dimension.
    _, idx3, *_ = PW.build_link(3.0)
    low_sector_dimension = len(idx3)
    high_spin_support = sum(
        1 for key in overcut["final_state"]
        if any(index >= low_sector_dimension for index in key)
    )
    checks = {
        "independent_exact_Haar_control_reproduces_first_five_values": exact_prefix == known,
        "independent_exact_sixth_control_is_integer": isinstance(exact, int),
        "Jmax_3_reaches_exact_spin_3": math.isclose(baseline["final_max_spin"], 3.0, rel_tol=0.0, abs_tol=1e-12),
        "Jmax_3_support_reaches_expected_frontier_at_every_power": all(
            math.isclose(row["max_spin_reached"], row["power"] / 2.0, rel_tol=0.0, abs_tol=1e-12)
            for row in baseline["rows"]
        ),
        "Jmax_7_2_has_no_support_above_3_at_any_power": all(
            row["max_spin_reached"] <= 3.0 + 1e-12 for row in overcut["rows"]
        ) and high_spin_support == 0,
        "Jmax_3_and_7_2_states_agree": (
            cutoff_comparison["left_only_support"] == 0
            and cutoff_comparison["right_only_support"] == 0
            and cutoff_comparison["max_amplitude_difference"] < 1e-10
            and cutoff_comparison["state_l2_difference"] < 1e-9
        ),
        "pruning_1e_14_stable_against_1e_13": (
            tight_comparison["max_amplitude_difference"] < 1e-9
            and tight_comparison["state_l2_difference"] < 1e-8
        ),
        "pruning_1e_12_stable_against_1e_13": (
            loose_comparison["max_amplitude_difference"] < 1e-9
            and loose_comparison["state_l2_difference"] < 1e-8
        ),
        "all_pruning_and_cutoff_norm_histories_match_exact_control": max(
            error for history in norm_history_errors.values() for error in history
        ) < 1e-3,
        "baseline_norm_matches_independent_exact_control": exact_error < 1e-3,
        "gauge_covariance_at_both_cutoffs": max(gauge_covariance.values()) < 1e-12,
        "all_final_amplitudes_finite": all(
            np.isfinite([z.real, z.imag]).all()
            for result_run in (baseline, tight, loose, overcut)
            for z in result_run["final_state"].values()
        ),
    }
    result = {
        "status": "sixth Peter-Weyl cutoff finite numerical audit",
        "passed": bool(all(checks.values())),
        "power": n,
        "required_cutoff_by_representation_bound": 3.0,
        "independent_exact_control": {
            "method": "SU(2)^3 Haar integral via unit-quaternion S^3 characteristic function and exact rational Laplacian recurrence",
            "known_first_five_values": known,
            "recomputed_first_five_values": exact_prefix,
            "recomputed_first_six_values": exact_values,
            "norm_history_absolute_errors": norm_history_errors,
            "sixth_exact_norm_squared": exact,
            "derivation": "W=2 sum_{i<j} q_i dot q_j=|q_1+q_2+q_3+1|^2-4; ||W^6|0>||^2=E[W^12]",
        },
        "baseline_Jmax_3_prune_1e-13": {
            key: value for key, value in baseline.items() if key != "final_state"
        },
        "tight_prune_1e-14": {
            key: value for key, value in tight.items() if key != "final_state"
        },
        "loose_prune_1e-12": {
            key: value for key, value in loose.items() if key != "final_state"
        },
        "overcutoff_Jmax_7_2_prune_1e-13": {
            key: value for key, value in overcut.items() if key != "final_state"
        },
        "cutoff_state_comparison_J3_vs_J7_2": cutoff_comparison,
        "pruning_state_comparison_1e-13_vs_1e-14": tight_comparison,
        "pruning_state_comparison_1e-13_vs_1e-12": loose_comparison,
        "overcutoff_support_above_spin_3": high_spin_support,
        "gauge_covariance_max_errors": gauge_covariance,
        "exact_norm_absolute_error": exact_error,
        "checks": checks,
        "claim_boundary": "Finite sixth-order Peter-Weyl cutoff calculation and an independent exact SU(2)^3 Haar moment. This is not a continuum-limit proof or a physical theorem.",
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    rendered = json.dumps(result, indent=2)
    print(rendered)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
