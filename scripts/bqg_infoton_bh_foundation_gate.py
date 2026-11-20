#!/usr/bin/env python3
"""Finite BQG bridge controls: spin-2/TT selection, Floquet control, unitary Page toy.

This is deliberately NOT a physical BQG black-hole or infoton simulation.
It tests finite algebra and conditional toy dynamics while recording the
missing physical coupling and physical-history inputs as null.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


def kron_on(op: np.ndarray, site: int, n: int = 4) -> np.ndarray:
    out = np.array([[1.0 + 0.0j]])
    eye = np.eye(2, dtype=complex)
    for i in range(n):
        out = np.kron(out, op if i == site else eye)
    return out


def spin2_gate() -> dict:
    sx = np.array([[0, 1], [1, 0]], complex) / 2
    sy = np.array([[0, -1j], [1j, 0]], complex) / 2
    sz = np.array([[1, 0], [0, -1]], complex) / 2
    J = [sum(kron_on(op, i) for i in range(4)) for op in (sx, sy, sz)]
    J2 = sum(x @ x for x in J)
    vals = np.linalg.eigvalsh(J2)
    mult = {str(x): int(np.count_nonzero(np.isclose(vals, x, atol=1e-10))) for x in (0, 2, 6)}
    plus, minus = np.zeros(16, complex), np.zeros(16, complex)
    plus[0], minus[-1] = 1, 1
    leak = []
    for state in (plus, minus):
        projector = np.linalg.eigh(J2)
        mask = np.isclose(projector[0], 6, atol=1e-10)
        P2 = projector[1][:, mask] @ projector[1][:, mask].conj().T
        leak.append(float(np.linalg.norm((np.eye(16) - P2) @ state)))
    mvals = [float(np.vdot(s, J[2] @ s).real) for s in (plus, minus)]
    return {
        "passed": mult == {"0": 2, "2": 9, "6": 5} and max(leak) < 1e-12 and np.allclose(mvals, [2, -2]),
        "J2_multiplicities": mult,
        "extremal_Jz": mvals,
        "projector_leakage": leak,
        "interpretation": "finite kinematic j=2 carrier; not a physical infoton pole",
    }


def tt_projector(k: np.ndarray) -> np.ndarray:
    k = k / np.linalg.norm(k)
    P = np.eye(3) - np.outer(k, k)
    out = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for a in range(3):
                for b in range(3):
                    out[i, j, a, b] = 0.5 * (P[i, a] * P[j, b] + P[i, b] * P[j, a] - P[i, j] * P[a, b])
    return out


def tt_selection_gate() -> dict:
    Pi = tt_projector(np.array([0.0, 0.0, 1.0]))
    scalar = np.eye(3)
    vector_source = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    tensor = np.diag([1.0, -1.0, 0.0])
    project = lambda T: np.einsum("ijab,ab->ij", Pi, T)
    norms = {
        "scalar_TT_norm": float(np.linalg.norm(project(scalar))),
        "longitudinal_vector_TT_norm": float(np.linalg.norm(project(vector_source))),
        "transverse_spin2_TT_norm": float(np.linalg.norm(project(tensor))),
    }
    passed = norms["scalar_TT_norm"] < 1e-12 and norms["longitudinal_vector_TT_norm"] < 1e-12 and norms["transverse_spin2_TT_norm"] > 1.0
    return {"passed": passed, **norms, "scope": "kinematic representation selection only"}


def rk4_monodromy(omega: float, drive: float, modulation: float, steps: int) -> np.ndarray:
    period = 2.0 * math.pi / drive
    dt = period / steps
    M = np.eye(2)

    def rhs(t: float, X: np.ndarray) -> np.ndarray:
        A = np.array([[0.0, 1.0], [-omega**2 * (1.0 + modulation * math.cos(drive * t)), 0.0]])
        return A @ X

    for i in range(steps):
        t = i * dt
        k1 = rhs(t, M)
        k2 = rhs(t + dt / 2, M + dt * k1 / 2)
        k3 = rhs(t + dt / 2, M + dt * k2 / 2)
        k4 = rhs(t + dt, M + dt * k3)
        M = M + dt * (k1 + 2*k2 + 2*k3 + k4) / 6
    return M


def floquet_mu(omega: float, drive: float, modulation: float, steps: int) -> float:
    M = rk4_monodromy(omega, drive, modulation, steps)
    period = 2.0 * math.pi / drive
    eig = np.linalg.eigvals(M)
    return max(0.0, float(np.max(np.log(np.abs(eig))) / period))


def floquet_gate() -> dict:
    omega, eps = 1.0, 0.02
    coarse = floquet_mu(omega, 2 * omega, eps, 2000)
    fine = floquet_mu(omega, 2 * omega, eps, 4000)
    leading = abs(eps) * omega / 4
    off = [floquet_mu(omega, x * omega, eps, 4000) for x in (1.8, 2.2)]
    rel = abs(fine - leading) / leading
    return {
        "passed": rel < 0.002 and abs(fine - coarse) / leading < 0.002 and max(off) < 1e-7,
        "toy_modulation_epsilon": eps,
        "xi_derived_from_BQG": None,
        "mu_coarse": coarse,
        "mu_fine": fine,
        "mu_leading": leading,
        "relative_error": rel,
        "off_resonance_mu": off,
        "scope": "conditional Mathieu/Floquet control; epsilon is imposed, not derived",
    }


def haar_unitary4(rng: np.random.Generator) -> np.ndarray:
    z = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    q, r = np.linalg.qr(z)
    phases = np.diag(r)
    phases = np.where(np.abs(phases) > 0, phases / np.abs(phases), 1.0)
    return q @ np.diag(np.conj(phases))


def apply_two_qubit(state: np.ndarray, U: np.ndarray, a: int, b: int, n: int) -> np.ndarray:
    axes = [a, b] + [i for i in range(n) if i not in (a, b)]
    inv = np.argsort(axes)
    tensor = state.reshape([2] * n).transpose(axes).reshape(4, -1)
    return (U @ tensor).reshape([2] * n).transpose(inv).reshape(-1)


def entropy_of_radiation(state: np.ndarray, n_bh: int, emitted: int) -> float:
    n_total = 2 * n_bh
    rad_axes = list(range(n_bh, n_bh + emitted))
    rest_axes = [i for i in range(n_total) if i not in rad_axes]
    psi = state.reshape([2] * n_total).transpose(rad_axes + rest_axes).reshape(2**emitted, -1)
    rho = psi @ psi.conj().T
    vals = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    vals = vals[vals > 1e-13]
    return float(-np.sum(vals * np.log(vals)) / math.log(2))


def unitary_page_toy(n_bh: int = 5, seed: int = 20261120) -> dict:
    rng = np.random.default_rng(seed)
    bh = rng.normal(size=2**n_bh) + 1j * rng.normal(size=2**n_bh)
    bh /= np.linalg.norm(bh)
    state = np.zeros(2 ** (2 * n_bh), complex)
    state[: 2**n_bh] = bh
    entropies = [0.0]
    for emitted in range(n_bh):
        active = n_bh - emitted
        for a in range(active):
            for b in range(a + 1, active):
                state = apply_two_qubit(state, haar_unitary4(rng), a, b, 2 * n_bh)
        state = apply_two_qubit(state, np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], complex), active - 1, n_bh + emitted, 2 * n_bh)
        entropies.append(entropy_of_radiation(state, n_bh, emitted + 1))
    norm_error = abs(float(np.vdot(state, state).real) - 1.0)
    passed = norm_error < 1e-10 and entropies[0] < 1e-10 and entropies[-1] < 1e-9 and max(entropies) > 0.5
    return {
        "passed": passed,
        "n_initial_BH_qubits": n_bh,
        "radiation_entropy_bits_by_emission": entropies,
        "global_norm_error": norm_error,
        "final_radiation_purity_entropy_bits": entropies[-1],
        "scope": "finite unitary scrambling-and-release control; not derived from BQG amplitudes or Hawking dynamics",
    }


def run() -> dict:
    spin = spin2_gate()
    tt = tt_selection_gate()
    floquet = floquet_gate()
    page = unitary_page_toy()
    return {
        "kind": "finite_foundation_controls_not_physical_BQG_completion",
        "passed": all(x["passed"] for x in (spin, tt, floquet, page)),
        "graviton_infoton": {"spin2": spin, "TT_selection": tt, "conditional_Floquet": floquet},
        "black_hole": {"unitary_Page_curve_control": page},
        "physical_inputs_still_missing": [
            "constraint-projected gauge-invariant infoton operator",
            "positive physical infoton spectral pole",
            "derived h-II vertex and xi",
            "full Lorentzian physical-history measure and metric",
            "dynamical trapping horizon and BQG area operator",
            "derived horizon entropy and coefficient 1/4",
            "Hawking Bogoliubov spectrum and greybody factors from BQG",
            "energy-conserving evaporation from the same amplitudes",
            "radiation density matrix from the same BQG evolution",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
