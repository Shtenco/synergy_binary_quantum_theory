#!/usr/bin/env python3
"""Regression gate for a discrete null-optical/Raychaudhuri operator.

Scope: finite-dimensional screen-space Jacobi/Riccati discretisation only.
This does NOT derive the optical tidal operator from the full BQG constraints.
Requires NumPy; the reference integrator is implemented here (no SciPy).
"""
from __future__ import annotations

import math
import numpy as np


def _assert_symmetric(name: str, matrix: np.ndarray, tol: float = 1e-12) -> None:
    defect = float(np.linalg.norm(matrix - matrix.T, ord="fro"))
    scale = max(1.0, float(np.linalg.norm(matrix, ord="fro")))
    assert defect <= tol * scale, (name, defect, scale)


def riccati_step(B: np.ndarray, R: np.ndarray, h: float) -> np.ndarray:
    """First-order discrete Riccati map induced by the Jacobi recurrence."""
    if h <= 0:
        raise ValueError("step size must be positive")
    _assert_symmetric("B", B)
    _assert_symmetric("R", R)
    A = np.eye(B.shape[0]) + h * B
    singular_values = np.linalg.svd(A, compute_uv=False)
    if singular_values[-1] < 1e-12 * max(1.0, singular_values[0]):
        raise FloatingPointError("Riccati chart reached a caustic/singular step")
    out = np.linalg.solve(A, B) - h * R
    _assert_symmetric("updated B", out, tol=5e-11)
    return out


def integrate_riccati(B0: np.ndarray, R: np.ndarray, T: float, h: float):
    n = int(round(T / h))
    if n < 1 or not math.isclose(n * h, T, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("T must be a positive integer multiple of h")
    B = np.array(B0, dtype=float)
    history = [B.copy()]
    for _ in range(n):
        B = riccati_step(B, R, h)
        history.append(B.copy())
    return history


def integrate_jacobi(B0: np.ndarray, R: np.ndarray, T: float, h: float):
    """Second-difference Jacobi recurrence; forward optical B is first order."""
    n = int(round(T / h))
    if n < 1 or not math.isclose(n * h, T, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("T must be a positive integer multiple of h")
    d = B0.shape[0]
    eye = np.eye(d)
    J_prev = eye.copy()
    J_cur = eye + h * B0
    Js = [J_prev.copy(), J_cur.copy()]
    Bs = [np.array(B0, dtype=float)]
    for _ in range(n):
        J_next = (2.0 * eye - h * h * R) @ J_cur - J_prev
        delta = (J_next - J_cur) / h
        B_next = np.linalg.solve(J_cur.T, delta.T).T
        Js.append(J_next.copy())
        Bs.append(B_next.copy())
        J_prev, J_cur = J_cur, J_next
    return Js, Bs


def integrate_jacobi_midpoint(B0: np.ndarray, R: np.ndarray, T: float, h: float):
    """Second-order midpoint optical matrix for constant R.

    J(h) is Taylor-started with J(0)=I, J'(0)=B0, J''(0)=-R.
    B_{n+1/2} = 2/h (J_{n+1}-J_n)(J_{n+1}+J_n)^(-1).
    """
    n = int(round(T / h))
    if n < 1 or not math.isclose(n * h, T, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("T must be a positive integer multiple of h")
    d = B0.shape[0]
    eye = np.eye(d)
    J_prev = eye.copy()
    J_cur = eye + h * B0 - 0.5 * h * h * R
    midpoints = []
    area_trace_defect = 0.0
    for j in range(n):
        S = J_cur + J_prev
        sv = np.linalg.svd(S, compute_uv=False)
        if sv[-1] < 1e-12 * max(1.0, sv[0]):
            raise FloatingPointError("midpoint Jacobi chart is singular")
        delta = J_cur - J_prev
        B_mid = (2.0 / h) * np.linalg.solve(S.T, delta.T).T
        _assert_symmetric("midpoint B", B_mid, tol=5e-9)
        logdet_cur = np.linalg.slogdet(J_cur)
        logdet_prev = np.linalg.slogdet(J_prev)
        if logdet_cur[0] <= 0 or logdet_prev[0] <= 0:
            raise FloatingPointError("Jacobi area crossed a caustic")
        area_rate = (logdet_cur[1] - logdet_prev[1]) / h
        area_trace_defect = max(
            area_trace_defect, abs(float(np.trace(B_mid)) - area_rate)
        )
        midpoints.append(B_mid)
        if j < n - 1:
            J_next = (2.0 * eye - h * h * R) @ J_cur - J_prev
            J_prev, J_cur = J_cur, J_next
    return midpoints, area_trace_defect


def rk4_reference(B0: np.ndarray, R: np.ndarray, T: float, steps: int = 40000):
    """High-resolution RK4 reference for B' = -B^2 - R."""
    h = T / steps
    B = np.array(B0, dtype=float)

    def f(X):
        return -X @ X - R

    for _ in range(steps):
        k1 = f(B)
        k2 = f(B + 0.5 * h * k1)
        k3 = f(B + 0.5 * h * k2)
        k4 = f(B + h * k3)
        B = B + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    return B


def exact_zero_tidal(B0: np.ndarray, T: float) -> np.ndarray:
    eye = np.eye(B0.shape[0])
    return np.linalg.solve((eye + T * B0).T, B0.T).T


def run() -> None:
    B0 = np.array([[0.70, 0.15], [0.15, 0.30]], dtype=float)
    # Vacuum screen tidal matrix: symmetric and trace-free (R_ab k^a k^b = 0).
    R = np.array([[0.20, 0.04], [0.04, -0.20]], dtype=float)
    T = 0.5

    # 1. Jacobi and Riccati discretisations are algebraically equivalent.
    jacobi_checks = 0
    for h in (0.1, 0.05, 0.025, 0.0125):
        Bs_ric = integrate_riccati(B0, R, T, h)
        _, Bs_jac = integrate_jacobi(B0, R, T, h)
        assert len(Bs_ric) == len(Bs_jac)
        for X, Y in zip(Bs_ric, Bs_jac):
            assert np.linalg.norm(X - Y, ord="fro") < 2e-9
            jacobi_checks += 1

    # 2. Symmetry is preserved.
    symmetry_defect = 0.0
    for h in (0.1, 0.05, 0.025, 0.0125):
        for B in integrate_riccati(B0, R, T, h):
            symmetry_defect = max(
                symmetry_defect, float(np.linalg.norm(B - B.T, ord="fro"))
            )
    assert symmetry_defect < 2e-12

    # 3. First-order convergence of the forward Riccati variable.
    reference_T = rk4_reference(B0, R, T)
    errors_first = []
    for h in (0.1, 0.05, 0.025, 0.0125):
        numerical = integrate_riccati(B0, R, T, h)[-1]
        errors_first.append(
            float(np.linalg.norm(numerical - reference_T, ord="fro"))
        )
    ratios_first = [
        errors_first[i] / errors_first[i+1] for i in range(len(errors_first)-1)
    ]
    assert all(1.85 < ratio < 2.15 for ratio in ratios_first), (
        errors_first, ratios_first
    )

    # 4. Second-order midpoint extraction from the Jacobi field.
    errors_mid = []
    area_defects = []
    for h in (0.1, 0.05, 0.025, 0.0125):
        midpoint_values, area_defect = integrate_jacobi_midpoint(B0, R, T, h)
        B_mid = midpoint_values[-1]
        reference_mid = rk4_reference(B0, R, T - 0.5*h)
        errors_mid.append(
            float(np.linalg.norm(B_mid - reference_mid, ord="fro"))
        )
        area_defects.append(area_defect)
    ratios_mid = [
        errors_mid[i] / errors_mid[i+1] for i in range(len(errors_mid)-1)
    ]
    assert all(3.5 < ratio < 4.5 for ratio in ratios_mid), (
        errors_mid, ratios_mid
    )
    area_ratios = [
        area_defects[i] / area_defects[i+1]
        for i in range(len(area_defects)-1)
    ]
    assert all(3.5 < ratio < 4.5 for ratio in area_ratios), (
        area_defects, area_ratios
    )

    # 5. Zero-tidal solution is exact up to floating-point roundoff.
    R0 = np.zeros((2, 2))
    zero_tidal_error = 0.0
    for h in (0.25, 0.1, 0.05, 0.025):
        got = integrate_riccati(B0, R0, T, h)[-1]
        exact = exact_zero_tidal(B0, T)
        zero_tidal_error = max(
            zero_tidal_error, float(np.linalg.norm(got - exact, ord="fro"))
        )
    assert zero_tidal_error < 2e-12

    # 6. Vacuum Raychaudhuri trace identity on a two-dimensional screen.
    theta = float(np.trace(B0))
    sigma = B0 - 0.5 * theta * np.eye(2)
    lhs = float(np.trace(B0 @ B0))
    rhs = 0.5 * theta * theta + float(np.trace(sigma @ sigma))
    trace_identity_error = abs(lhs - rhs)
    vacuum_trace_error = abs(float(np.trace(R)))
    dB = -B0 @ B0 - R
    dtheta = float(np.trace(dB))
    raychaudhuri_rhs = -0.5 * theta * theta - float(np.trace(sigma @ sigma))
    vacuum_raychaudhuri_error = abs(dtheta - raychaudhuri_rhs)
    assert trace_identity_error < 1e-14
    assert vacuum_trace_error < 1e-14
    assert vacuum_raychaudhuri_error < 1e-14

    # 7. Detect the Riccati chart singularity at a conjugate point.
    caustic_detected = False
    try:
        integrate_riccati(-np.eye(2), R0, 1.0, 0.1)
    except FloatingPointError:
        caustic_detected = True
    assert caustic_detected

    print("Jacobi/Riccati equivalence checks =", jacobi_checks)
    print("maximum symmetry defect =", f"{symmetry_defect:.3e}")
    print("vacuum R_kk trace error =", f"{vacuum_trace_error:.3e}")
    print("vacuum Raychaudhuri RHS error =", f"{vacuum_raychaudhuri_error:.3e}")
    print("first-order errors =", [f"{x:.12e}" for x in errors_first])
    print("first-order halving ratios =", [f"{x:.8f}" for x in ratios_first])
    print("midpoint errors =", [f"{x:.12e}" for x in errors_mid])
    print("midpoint halving ratios =", [f"{x:.8f}" for x in ratios_mid])
    print("midpoint area/trace defects =", [f"{x:.12e}" for x in area_defects])
    print("midpoint area/trace halving ratios =", [f"{x:.8f}" for x in area_ratios])
    print("zero-tidal exact-solution error =", f"{zero_tidal_error:.3e}")
    print("2D Raychaudhuri trace-identity error =", f"{trace_identity_error:.3e}")
    print("caustic chart singularity detection =", caustic_detected)
    print("DISCRETE OPTICAL RICCATI GATE: PASS")
    print("SCOPE: finite-dimensional test model; no full-BQG derivation claimed")


if __name__ == "__main__":
    run()
