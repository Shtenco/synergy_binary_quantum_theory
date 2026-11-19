#!/usr/bin/env python3
"""Mirror-covariant Wilson-Dirac matter carrier.

This finite gate supplies a concrete microscopic H_m(sigma) candidate rather
than only an abstract rest-mass function. On a three-dimensional spatial
lattice,

  H_q(k,sigma) = sum_i alpha_i sin(k_i)
               + beta_D [m_q(sigma) + r_W sum_i (1-cos(k_i))].

The mirror mass law obeys m_q(sigma)=m_-q(-sigma). The gate checks:
- exact mirror-partner spectral equality;
- positive aligned rest masses;
- opposite Hellmann-Feynman rest sources;
- the moving-state derivative identity;
- the standard massless Wilson corner control: one k=0 zero and no corner
  doublers.

The dimensionless mirror coupling beta remains an input. This gate shows that
the required source structure is compatible with a concrete Wilson-Dirac
carrier; it does not derive realistic gauge/chiral matter or beta numerically.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


TOL = 2e-12


def gamma_data():
    I = np.eye(2, dtype=complex)
    Z = np.zeros((2, 2), dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    g0 = np.block([[I, Z], [Z, -I]])

    def gi(s):
        return np.block([[Z, s], [-s, Z]])

    g = (g0, gi(sx), gi(sy), gi(sz))
    alpha = tuple(g0 @ g[i] for i in range(1, 4))
    return alpha, g0


ALPHA, BETA_D = gamma_data()


def mirror_mass(mstar, beta, sigma, q):
    return mstar * math.exp(q * beta * sigma)


def hamiltonian(k, mass, wilson_r=1.0):
    M = mass + wilson_r * sum(1.0 - math.cos(x) for x in k)
    H = sum(math.sin(k[i]) * ALPHA[i] for i in range(3)) + M * BETA_D
    return H, M


def positive_state(H):
    vals, vecs = np.linalg.eigh(H)
    pos = np.flatnonzero(vals > 0)
    if not len(pos):
        raise RuntimeError("no positive-energy state")
    i = int(pos[0])
    return float(vals[i]), vecs[:, i], vals


def aligned_mirror_gate(mstar=0.4, beta=0.37, v=0.8, wilson_r=1.0):
    branches = {}
    spectra = []
    sources = []
    max_spectrum_difference = 0.0

    momenta = [
        (0.0, 0.0, 0.0),
        (0.17, 0.31, 0.23),
        (0.7, 0.2, 0.4),
        (1.1, 0.8, 0.3),
    ]

    for chi in (+1, -1):
        q = chi
        sigma = chi * v
        mass = mirror_mass(mstar, beta, sigma, q)
        rows = []
        for k in momenta:
            H, M = hamiltonian(k, mass, wilson_r)
            E, state, eig = positive_state(H)
            dm = q * beta * mass
            source = float(np.real(np.vdot(state, (dm * BETA_D) @ state)))
            target = dm * M / E
            rows.append({
                "k": list(k),
                "wilson_effective_mass": M,
                "positive_energy": E,
                "source": source,
                "target_dm_times_M_over_E": target,
                "source_error": abs(source - target),
                "spectrum": eig.tolist(),
            })
        spectra.append(rows)

        H0, _ = hamiltonian((0.0, 0.0, 0.0), mass, wilson_r)
        E0, state0, _ = positive_state(H0)
        dm = q * beta * mass
        source0 = float(np.real(np.vdot(state0, (dm * BETA_D) @ state0)))
        sources.append(source0)
        branches[str(chi)] = {
            "chi": chi,
            "q": q,
            "sigma": sigma,
            "rest_mass": mass,
            "rest_positive_energy": E0,
            "rest_source": source0,
            "target_chi_beta_m": chi * beta * mass,
            "beta_m_recovered": source0 / (chi * mass),
            "momenta": rows,
        }

    for rp, rm in zip(spectra[0], spectra[1]):
        max_spectrum_difference = max(
            max_spectrum_difference,
            float(np.max(np.abs(np.asarray(rp["spectrum"]) - np.asarray(rm["spectrum"])))),
        )

    max_source_error = max(
        row["source_error"] for branch in branches.values() for row in branch["momenta"]
    )

    return {
        "branches": branches,
        "mirror_rest_mass_difference": abs(branches["1"]["rest_mass"] - branches["-1"]["rest_mass"]),
        "mirror_rest_source_sum": abs(sources[0] + sources[1]),
        "max_mirror_spectrum_difference": max_spectrum_difference,
        "max_moving_source_error": max_source_error,
    }


def wilson_corner_control(wilson_r=1.0):
    corners = list(itertools.product((0.0, math.pi), repeat=3))
    rows = []
    zeros = []
    for k in corners:
        H, M = hamiltonian(k, 0.0, wilson_r)
        eig = np.linalg.eigvalsh(H)
        min_abs = float(np.min(np.abs(eig)))
        if min_abs < 1e-12:
            zeros.append(k)
        rows.append({
            "k": list(k),
            "wilson_mass": M,
            "min_abs_energy": min_abs,
        })
    return {
        "corner_count": len(corners),
        "zero_corners": [list(k) for k in zeros],
        "zero_corner_count": len(zeros),
        "only_physical_corner_zero": zeros == [(0.0, 0.0, 0.0)],
        "rows": rows,
    }


def run():
    mirror = aligned_mirror_gate()
    wilson = wilson_corner_control()
    plus = mirror["branches"]["1"]
    minus = mirror["branches"]["-1"]

    passed = (
        mirror["mirror_rest_mass_difference"] < TOL
        and mirror["mirror_rest_source_sum"] < TOL
        and mirror["max_mirror_spectrum_difference"] < TOL
        and mirror["max_moving_source_error"] < TOL
        and plus["rest_mass"] > 0.0
        and minus["rest_mass"] > 0.0
        and abs(plus["beta_m_recovered"] - 0.37) < TOL
        and abs(minus["beta_m_recovered"] - 0.37) < TOL
        and wilson["only_physical_corner_zero"]
    )

    return {
        "status": "mirror-covariant Wilson-Dirac matter gate",
        "passed": bool(passed),
        "hamiltonian": (
            "H_q(k,sigma)=sum_i alpha_i sin(k_i)+beta_D[m_q(sigma)+r_W sum_i(1-cos(k_i))]"
        ),
        "mirror_mass_law": "m_q(sigma)=m_* exp(q beta sigma)",
        "mirror": mirror,
        "wilson_massless_corner_control": wilson,
        "main_result": (
            "A concrete positive-energy Wilson-Dirac carrier can realize equal mirror-partner spectra and "
            "opposite static sigma charges while retaining the Wilson removal of Brillouin-corner doublers."
        ),
        "scope": (
            "beta is still a microscopic input, Wilson fermions break naive chiral symmetry at finite spacing, "
            "and this regular-lattice carrier is not yet the realistic irregular PL/Peter-Weyl matter sector."
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
