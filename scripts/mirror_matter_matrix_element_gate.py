#!/usr/bin/env python3
"""Microscopic mirror-matter matrix-element gate.

This gate separates two questions that were previously conflated:

1. Does the already proposed axial coupling sigma J5^0 provide a static
   charge for cold massive matter?
2. What is the minimal positive-mass mirror-doublet matter law that does
   provide the required static source Q_sigma = chi beta_m m?

The answer to (1) is no for a free on-shell massive Dirac particle at rest:
J5^0/J^0 = h |p|/E and vanishes at p=0; the unpolarized average vanishes
for every p. The diagonal pseudoscalar bilinear ubar i gamma5 u also vanishes.

For (2), let q=+/-1 be an internal mirror label with
    (sigma,q) -> (-sigma,-q).
Demand a positive rest mass and a constant universal logarithmic response
    d ln m_q / d sigma = q beta.
The unique solution is
    m_q(sigma) = m_* exp(q beta sigma).
For mirror-related aligned states (sigma,q)=(chi v,chi), both branches have
the same positive mass but opposite static scalar charge:
    d m/d sigma = chi beta m.
Thus beta_m=beta exactly in this minimal matter extension.

This is a finite/QM gate for the matrix-element logic. It does not derive the
numerical beta from the existing geometry-only theory; that requires an actual
microscopic matter Hamiltonian/spectrum.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


TOL = 5e-12


def gamma_matrices():
    I2 = np.eye(2, dtype=complex)
    Z2 = np.zeros((2, 2), dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)

    g0 = np.block([[I2, Z2], [Z2, -I2]])

    def gi(s):
        return np.block([[Z2, s], [-s, Z2]])

    g1, g2, g3 = gi(sx), gi(sy), gi(sz)
    g5 = 1j * g0 @ g1 @ g2 @ g3
    return (g0, g1, g2, g3), g5, (sx, sy, sz)


GAMMA, G5, PAULI = gamma_matrices()
G0 = GAMMA[0]
ALPHA = tuple(G0 @ GAMMA[i] for i in range(1, 4))
BETA_D = G0


def helicity_two_spinor(pvec, helicity):
    p = np.asarray(pvec, dtype=float)
    pnorm = float(np.linalg.norm(p))
    if pnorm < 1e-15:
        return np.array([1.0, 0.0], dtype=complex) if helicity > 0 else np.array([0.0, 1.0], dtype=complex)
    Hh = sum((p[i] / pnorm) * PAULI[i] for i in range(3))
    vals, vecs = np.linalg.eigh(Hh)
    idx = int(np.argmax(vals) if helicity > 0 else np.argmin(vals))
    return vecs[:, idx]


def dirac_u(pvec, mass, helicity):
    p = np.asarray(pvec, dtype=float)
    E = math.sqrt(mass * mass + float(p @ p))
    chi = helicity_two_spinor(p, helicity)
    sigma_p = sum(p[i] * PAULI[i] for i in range(3))
    pref = math.sqrt(E + mass)
    lower = sigma_p @ chi / (E + mass)
    u = pref * np.concatenate([chi, lower])
    return u, E


def dirac_bilinears(pvec, mass, helicity):
    u, E = dirac_u(pvec, mass, helicity)
    ubar = u.conj().T @ G0
    J0 = float(np.real(ubar @ G0 @ u))
    A0 = float(np.real(ubar @ G0 @ G5 @ u))
    scalar = float(np.real(ubar @ u))
    pseudoscalar = complex(ubar @ (1j * G5) @ u)
    return {
        "E": E,
        "J0": J0,
        "A0": A0,
        "A0_over_J0": A0 / J0,
        "scalar": scalar,
        "pseudoscalar_abs": float(abs(pseudoscalar)),
    }


def axial_gate(mass=2.0):
    rows = []
    max_formula_error = 0.0
    max_pseudoscalar = 0.0
    for p in (0.0, 0.1, 0.3, 1.0, 3.0):
        by_h = {}
        for h in (+1, -1):
            b = dirac_bilinears((0.0, 0.0, p), mass, h)
            target = h * p / b["E"]
            max_formula_error = max(max_formula_error, abs(b["A0_over_J0"] - target))
            max_pseudoscalar = max(max_pseudoscalar, b["pseudoscalar_abs"])
            by_h[str(h)] = b
        unpolarized_A0 = 0.5 * (by_h["1"]["A0"] + by_h["-1"]["A0"])
        rows.append({
            "p_over_m": p / mass,
            "helicities": by_h,
            "unpolarized_A0": unpolarized_A0,
        })

    rest = rows[0]
    rest_zero = all(abs(rest["helicities"][str(h)]["A0"]) < TOL for h in (+1, -1))
    unpolarized_zero = all(abs(r["unpolarized_A0"]) < TOL for r in rows)
    return {
        "rows": rows,
        "identity": "J5^0/J^0 = helicity*|p|/E",
        "rest_axial_density_zero": rest_zero,
        "unpolarized_axial_density_zero": unpolarized_zero,
        "max_identity_error": max_formula_error,
        "max_diagonal_pseudoscalar_abs": max_pseudoscalar,
        "static_mass_charge_from_axial_channel": 0.0,
    }


def dirac_hamiltonian(pvec, mass):
    p = np.asarray(pvec, dtype=float)
    return sum(p[i] * ALPHA[i] for i in range(3)) + mass * BETA_D


def positive_energy_state(pvec, mass):
    H = dirac_hamiltonian(pvec, mass)
    vals, vecs = np.linalg.eigh(H)
    pos = np.flatnonzero(vals > 0)
    idx = int(pos[0])
    return float(vals[idx]), vecs[:, idx]


def mirror_mass(mstar, beta, sigma, q):
    return mstar * math.exp(q * beta * sigma)


def mirror_doublet_gate(mstar=2.1, beta=0.37, v=0.8):
    aligned = {}
    source_errors = []
    mirror_masses = []
    mirror_sources = []

    for chi in (+1, -1):
        q = chi
        sigma = chi * v
        mass = mirror_mass(mstar, beta, sigma, q)
        E0, state = positive_energy_state((0.0, 0.0, 0.0), mass)
        dm_dsigma = q * beta * mass
        dH_dsigma = dm_dsigma * BETA_D
        source = float(np.real(np.vdot(state, dH_dsigma @ state)))
        target = chi * beta * mass
        source_errors.append(abs(source - target))
        mirror_masses.append(mass)
        mirror_sources.append(source)
        aligned[str(chi)] = {
            "chi": chi,
            "q": q,
            "sigma": sigma,
            "mass": mass,
            "positive_energy_rest_E": E0,
            "hellmann_feynman_source": source,
            "target_chi_beta_m": target,
            "beta_m_recovered": source / (chi * mass),
        }

    moving = []
    max_moving_error = 0.0
    q = +1
    sigma = +v
    mass = mirror_mass(mstar, beta, sigma, q)
    for p in (0.0, 0.2, 0.7, 2.0):
        E, state = positive_energy_state((0.0, 0.0, p), mass)
        dm_dsigma = q * beta * mass
        source = float(np.real(np.vdot(state, (dm_dsigma * BETA_D) @ state)))
        target = q * beta * mass * mass / E
        max_moving_error = max(max_moving_error, abs(source - target))
        moving.append({
            "p": p,
            "E": E,
            "source": source,
            "target_q_beta_m2_over_E": target,
            "source_over_rest_charge": source / (q * beta * mass),
        })

    sigma_probe = 0.63
    mplus = mirror_mass(mstar, beta, sigma_probe, +1)
    mminus = mirror_mass(mstar, beta, sigma_probe, -1)
    beta_from_mass_ratio = math.log(mplus / mminus) / (2.0 * sigma_probe)

    positive_scan = min(
        mirror_mass(mstar, beta, s, q)
        for s in np.linspace(-4.0, 4.0, 81)
        for q in (+1, -1)
    )

    return {
        "mirror_rule": "(sigma,q)->(-sigma,-q)",
        "mass_law": "m_q(sigma)=m_* exp(q beta sigma)",
        "constant_log_response": "d ln m_q/d sigma = q beta",
        "aligned_mirror_states": aligned,
        "mirror_partner_mass_difference": abs(mirror_masses[0] - mirror_masses[1]),
        "mirror_partner_source_sum": abs(mirror_sources[0] + mirror_sources[1]),
        "max_rest_source_error": max(source_errors),
        "moving_source_control": moving,
        "max_moving_source_error": max_moving_error,
        "beta_input": beta,
        "beta_from_mass_ratio": beta_from_mass_ratio,
        "beta_recovery_error": abs(beta_from_mass_ratio - beta),
        "minimum_mass_on_sigma_scan": positive_scan,
        "derived_static_beta_m": beta,
        "force_strength_formula": "alpha=beta^2/(4*pi*G*Z_sigma)",
        "regular_seed_formula": "alpha=3*beta^2*ell/(8*sqrt(2)*pi*G*J)",
    }


def run():
    axial = axial_gate()
    mirror = mirror_doublet_gate()

    passed = (
        axial["rest_axial_density_zero"]
        and axial["unpolarized_axial_density_zero"]
        and axial["max_identity_error"] < 1e-12
        and axial["max_diagonal_pseudoscalar_abs"] < 1e-12
        and mirror["mirror_partner_mass_difference"] < 1e-12
        and mirror["mirror_partner_source_sum"] < 1e-12
        and mirror["max_rest_source_error"] < 1e-12
        and mirror["max_moving_source_error"] < 1e-12
        and mirror["beta_recovery_error"] < 1e-12
        and mirror["minimum_mass_on_sigma_scan"] > 0.0
    )

    return {
        "status": "mirror matter matrix-element gate",
        "passed": bool(passed),
        "axial_channel": axial,
        "positive_mass_mirror_doublet": mirror,
        "main_result": (
            "The previously proposed axial density is not a static charge for cold unpolarized massive matter. "
            "A minimal mirror doublet with positive mass law m_q=m_* exp(q beta sigma) has equal mirror-partner "
            "masses and opposite rest-frame Hellmann-Feynman charges Q_sigma=chi beta m, so beta_m=beta."
        ),
        "falsifier": (
            "If the actual microscopic matter Hamiltonian has zero mirror-odd logarithmic mass derivative "
            "d ln m/d sigma for every physical state, then beta_m=0 and the static mirror-force branch fails."
        ),
        "scope": (
            "This gate derives/extracts the matrix-element form but does not determine beta numerically from the "
            "current geometry-only repository. That requires a realistic microscopic matter Hamiltonian."
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = run()
    text = json.dumps(out, indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
