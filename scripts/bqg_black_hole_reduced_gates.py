#!/usr/bin/env python3
"""BH-1..BH-7 finite/reduced controls for the BQG black-hole programme.

Every imported semiclassical ingredient is explicitly marked. This script
validates numerical machinery and consistency relations; it does NOT derive
a black hole, Hawking radiation, or a Page curve from physical BQG amplitudes.
Requires only NumPy.
"""
from __future__ import annotations

import json
import math
import numpy as np


def trapping_horizon_gate() -> dict:
    """Locate theta_out=0 in Schwarzschild Painleve-Gullstrand slices."""
    mass = 1.0
    rows = []
    for cells in (32, 64, 128, 256, 512):
        r = np.linspace(1.5 * mass, 3.0 * mass, cells + 1)
        flow = np.sqrt(2.0 * mass / r)
        theta_out = 2.0 * (1.0 - flow) / r
        crossings = np.flatnonzero(theta_out[:-1] * theta_out[1:] <= 0.0)
        if len(crossings) != 1:
            raise RuntimeError("expected exactly one marginal-surface crossing")
        i = int(crossings[0])
        rh = r[i] - theta_out[i] * (r[i + 1] - r[i]) / (theta_out[i + 1] - theta_out[i])
        flow_h = math.sqrt(2.0 * mass / rh)
        theta_in = -2.0 * (1.0 + flow_h) / rh
        rows.append({"radial_cells": cells, "r_h": float(rh),
                     "relative_radius_error": abs(float(rh) - 2 * mass) / (2 * mass),
                     "theta_in_at_root": theta_in})
    errors = [x["relative_radius_error"] for x in rows]
    order = math.log(errors[-2] / errors[-1], 2)
    passed = all(x["theta_in_at_root"] < 0 for x in rows) and errors[-1] < 2e-6 and 1.9 < order < 2.1
    return {
        "passed": passed, "benchmark": "Schwarzschild in Painleve-Gullstrand coordinates",
        "null_expansions": "theta_out=2(1-sqrt(2M/r))/r; theta_in=-2(1+sqrt(2M/r))/r",
        "refinement": rows, "observed_order_last_pair": order,
        "scope": "trapping-horizon finder validated on an inserted GR benchmark; not a BQG emergent horizon",
    }


def icosphere(level: int, radius: float = 1.0) -> tuple[int, int, float]:
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    vertices = np.array([
        (-1, phi, 0), (1, phi, 0), (-1, -phi, 0), (1, -phi, 0),
        (0, -1, phi), (0, 1, phi), (0, -1, -phi), (0, 1, -phi),
        (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1), (-phi, 0, 1),
    ], dtype=float)
    vertices = vertices / np.linalg.norm(vertices, axis=1)[:, None] * radius
    faces = np.array([
        [0,11,5],[0,5,1],[0,1,7],[0,7,10],[0,10,11],
        [1,5,9],[5,11,4],[11,10,2],[10,7,6],[7,1,8],
        [3,9,4],[3,4,2],[3,2,6],[3,6,8],[3,8,9],
        [4,9,5],[2,4,11],[6,2,10],[8,6,7],[9,8,1],
    ], dtype=int)
    for _ in range(level):
        midpoint: dict[tuple[int, int], int] = {}
        new_vertices = list(vertices)

        def mid(a: int, b: int) -> int:
            key = tuple(sorted((int(a), int(b))))
            if key not in midpoint:
                point = (vertices[a] + vertices[b]) / 2.0
                point = point / np.linalg.norm(point) * radius
                midpoint[key] = len(new_vertices)
                new_vertices.append(point)
            return midpoint[key]

        new_faces = []
        for a, b, c in faces:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            new_faces.extend(([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]))
        vertices, faces = np.asarray(new_vertices), np.asarray(new_faces, dtype=int)
    tri = vertices[faces]
    areas = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2.0
    return len(vertices), len(faces), float(np.sum(areas))


def area_operator_gate() -> dict:
    """Separate Regge triangle area from the SU(2) Casimir spectrum."""
    radius = 1.0
    exact = 4.0 * math.pi * radius**2
    refinements = []
    for level in range(5):
        nv, nf, area = icosphere(level, radius)
        refinements.append({"level": level, "vertices": nv, "triangles": nf,
                            "regge_area": area, "relative_error": abs(area - exact) / exact})
    errors = [x["relative_error"] for x in refinements]
    order = math.log(errors[-2] / errors[-1], 2)
    j = 0.5
    casimir = math.sqrt(j * (j + 1.0))
    area_scale_over_lP2 = None
    puncture_area = "a_* sqrt(j(j+1)); a_* is not fixed by the current BQG physicalization"
    passed = errors[-1] < 0.002 and 1.9 < order < 2.1 and abs(casimir - math.sqrt(0.75)) < 1e-14
    return {
        "passed": passed, "exact_sphere_area": exact, "triangulation_refinement": refinements,
        "observed_area_error_order": order, "SU2_area_spectrum": {
            "operator_ansatz": "A_H = a_* sum_p sqrt(j_p(j_p+1))",
            "j_half_Casimir_sqrt": casimir,
            "area_scale_in_Planck_units": area_scale_over_lP2,
            "status": "Casimir dependence follows SU(2); a_* and horizon puncture map remain undetermined",
        },
        "scope": "Regge surface-area convergence plus kinematic SU(2) Casimir control; no BQG horizon area operator derived",
    }


def entropy_coefficient_gate() -> dict:
    """Exact fixed-j counting exposes the unresolved area normalization."""
    j = 0.5
    degeneracy = int(2 * j + 1)
    casimir = math.sqrt(j * (j + 1.0))
    entropy_per_puncture = math.log(degeneracy)
    dimensionless_slope = entropy_per_puncture / casimir
    required_a_over_lP2 = 4.0 * dimensionless_slope
    rows = []
    for punctures in (8, 16, 32, 64, 128, 256):
        S = punctures * entropy_per_puncture
        Abar = punctures * casimir
        rows.append({"punctures": punctures, "log_state_count": S,
                     "dimensionless_area": Abar, "S_over_Abar": S / Abar})
    passed = max(abs(x["S_over_Abar"] - dimensionless_slope) for x in rows) < 1e-13
    return {
        "passed": passed, "fixed_spin": j, "degeneracy_per_puncture": degeneracy,
        "entropy_per_puncture_nats": entropy_per_puncture,
        "dimensionless_S_over_A_slope": dimensionless_slope,
        "a_star_over_lP2_required_to_match_BH_1_over_4": required_a_over_lP2,
        "refinement_in_puncture_number": rows,
        "physical_Bekenstein_Hawking_coefficient_derived": False,
        "reason": "The counting gives S/Abar exactly, but converting Abar to physical area requires a_*; choosing a_* to hit 1/4 would be calibration, not derivation.",
        "scope": "exact fixed-j microstate count and explicit no-fit obstruction",
    }


def hawking_spectrum_gate() -> dict:
    """Check Planck spectrum from the standard exponential peeling map."""
    # In units hbar=kB=c=G=1: U=U_H-A exp(-kappa u) gives T=kappa/(2 pi).
    mass = 1.0
    kappa = 1.0 / (4.0 * mass)
    temperature = kappa / (2.0 * math.pi)
    exact_integral = math.pi**2 / 6.0
    rows = []
    for n in (100, 200, 400, 800, 1600):
        x = np.linspace(1e-6, 40.0, n + 1)
        integrand = x / np.expm1(x)
        integral = float(np.sum((integrand[1:] + integrand[:-1]) * np.diff(x) / 2.0))
        flux = temperature**2 * integral / (2.0 * math.pi)
        rows.append({"quadrature_points": n + 1, "dimensionless_integral": integral,
                     "relative_integral_error": abs(integral - exact_integral) / exact_integral,
                     "chiral_1plus1_flux": flux})
    exact_flux = math.pi * temperature**2 / 12.0
    final_error = abs(rows[-1]["chiral_1plus1_flux"] - exact_flux) / exact_flux
    passed = abs(temperature - 1.0 / (8.0 * math.pi * mass)) < 1e-14 and final_error < 2e-5
    return {
        "passed": passed, "units": "G=hbar=c=kB=1", "surface_gravity": kappa,
        "temperature": temperature, "expected_schwarzschild_temperature": 1.0 / (8.0 * math.pi * mass),
        "spectrum": "n_omega = 1/(exp(omega/T_H)-1)",
        "quadrature": rows, "exact_chiral_1plus1_flux": exact_flux,
        "greybody": "Gamma_omega=1 in a single chiral 1+1D channel; 4D greybody factors not computed",
        "scope": "standard semiclassical exponential-ray-tracing control; not derived from BQG amplitudes",
    }


def evaporation_gate() -> dict:
    """Integrate the 1+1D benchmark mass-loss ODE with energy bookkeeping."""
    alpha = 1.0 / (768.0 * math.pi)
    m0, mf = 1.0, 0.2
    exact_time = (m0**3 - mf**3) / (3.0 * alpha)
    results = []
    for steps in (100, 500, 1000, 2000):
        dt = exact_time / steps
        mass, radiated = m0, 0.0

        def dm(m: float) -> float:
            return -alpha / (m * m)

        for _ in range(steps):
            k1 = dm(mass)
            m2 = mass + dt * k1 / 2.0
            k2 = dm(m2)
            m3 = mass + dt * k2 / 2.0
            k3 = dm(m3)
            m4 = mass + dt * k3
            k4 = dm(m4)
            p1, p2, p3, p4 = alpha / mass**2, alpha / m2**2, alpha / m3**2, alpha / m4**2
            radiated += dt * (p1 + 2*p2 + 2*p3 + p4) / 6.0
            mass += dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0
        results.append({"rk4_steps": steps, "final_mass": mass,
                        "radiated_energy": radiated,
                        "energy_balance_error": abs((m0 - mass) - radiated),
                        "final_mass_error": abs(mass - mf)})
    passed = results[-1]["energy_balance_error"] < 1e-8 and results[-1]["final_mass_error"] < 1e-8
    return {
        "passed": passed, "mass_loss_law": "dM/dt=-alpha/M^2",
        "alpha_from_chiral_1plus1_flux": alpha, "analytic_lifetime_to_zero": m0**3 / (3.0 * alpha),
        "integrated_interval": {"M_initial": m0, "M_final_target": mf, "time": exact_time},
        "convergence": results,
        "scope": "energy-conserving integration of the reduced flux law; coefficient is 1+1D toy input, not 4D BQG evaporation",
    }


def haar_unitary4(rng: np.random.Generator) -> np.ndarray:
    z = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    q, r = np.linalg.qr(z)
    d = np.diag(r)
    phases = np.where(np.abs(d) > 0.0, d / np.abs(d), 1.0)
    return q @ np.diag(np.conj(phases))


def apply_two_qubit(state: np.ndarray, unitary: np.ndarray, a: int, b: int, n: int) -> np.ndarray:
    axes = [a, b] + [i for i in range(n) if i not in (a, b)]
    inverse = np.argsort(axes)
    tensor = state.reshape([2] * n).transpose(axes).reshape(4, -1)
    return (unitary @ tensor).reshape([2] * n).transpose(inverse).reshape(-1)


def radiation_entropy(state: np.ndarray, n_bh: int, emitted: int) -> float:
    radiation_axes = list(range(n_bh, n_bh + emitted))
    rest_axes = [i for i in range(2*n_bh) if i not in radiation_axes]
    psi = state.reshape([2] * (2*n_bh)).transpose(radiation_axes + rest_axes).reshape(2**emitted, -1)
    rho = psi @ psi.conj().T
    eigenvalues = np.linalg.eigvalsh((rho + rho.conj().T) / 2.0)
    eigenvalues = eigenvalues[eigenvalues > 1e-13]
    return float(-np.sum(eigenvalues * np.log2(eigenvalues)))


def page_curve_once(n: int, seed: int) -> tuple[list[float], float]:
    rng = np.random.default_rng(seed)
    bh = rng.normal(size=2**n) + 1j * rng.normal(size=2**n)
    bh /= np.linalg.norm(bh)
    state = np.zeros(2**(2*n), dtype=complex)
    state[:2**n] = bh
    entropies = [0.0]
    swap = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], dtype=complex)
    for emitted in range(n):
        active = n - emitted
        for a in range(active):
            for b in range(a + 1, active):
                state = apply_two_qubit(state, haar_unitary4(rng), a, b, 2*n)
        state = apply_two_qubit(state, swap, active - 1, n + emitted, 2*n)
        entropies.append(radiation_entropy(state, n, emitted + 1))
    norm_error = abs(float(np.vdot(state, state).real) - 1.0)
    return entropies, norm_error


def page_curve_gate() -> dict:
    rows = []
    for n, seed in ((4, 4101), (5, 5101), (6, 6101), (6, 6102), (6, 6103)):
        values, norm_error = page_curve_once(n, seed)
        peak = int(np.argmax(values))
        rows.append({"initial_BH_qubits": n, "seed": seed, "radiation_entropy_bits": values,
                     "peak_emission_step": peak, "peak_entropy_bits": max(values),
                     "final_entropy_bits": values[-1], "global_norm_error": norm_error,
                     "peak_is_interior": 0 < peak < n})
    passed = all(x["peak_is_interior"] and x["final_entropy_bits"] < 1e-9 and
                 x["global_norm_error"] < 1e-10 for x in rows)
    return {
        "passed": passed, "ensemble": "finite Haar two-qubit scrambling followed by sequential swap-out",
        "runs": rows,
        "scope": "unitary Page-curve control with exact global norm and final radiation purity; no BQG dynamics or Hawking channel",
    }


def refinement_stability_gate(parts: dict) -> dict:
    horizon_errors = [x["relative_radius_error"] for x in parts["BH1"]["refinement"]]
    area_errors = [x["relative_error"] for x in parts["BH2"]["triangulation_refinement"]]
    flux_errors = [x["relative_integral_error"] for x in parts["BH4"]["quadrature"]]
    evaporation_errors = [x["energy_balance_error"] for x in parts["BH5"]["convergence"]]
    passed = (horizon_errors[-1] < horizon_errors[0] / 100.0 and
              area_errors[-1] < area_errors[0] / 100.0 and
              flux_errors[-1] < flux_errors[0] / 100.0 and
              evaporation_errors[-1] < evaporation_errors[0] / 100.0)
    return {
        "passed": passed,
        "horizon_error_reduction": horizon_errors[0] / horizon_errors[-1],
        "area_error_reduction": area_errors[0] / area_errors[-1],
        "flux_quadrature_error_reduction": flux_errors[0] / flux_errors[-1],
        "evaporation_energy_error_reduction": evaporation_errors[0] / evaporation_errors[-1],
        "scope": "numerical regulator/refinement stability of the reduced controls only; no physical BQG regulator independence",
    }


def run() -> dict:
    parts = {
        "BH1": trapping_horizon_gate(),
        "BH2": area_operator_gate(),
        "BH3": entropy_coefficient_gate(),
        "BH4": hawking_spectrum_gate(),
        "BH5": evaporation_gate(),
        "BH6": page_curve_gate(),
    }
    parts["BH7"] = refinement_stability_gate(parts)
    return {
        "kind": "BQG_black_hole_reduced_controls_not_physical_completion",
        "passed_reduced_controls": all(x["passed"] for x in parts.values()),
        "gates": parts,
        "physical_BQG_gates_still_open": [
            "derive Lorentzian metric and physical history from BQG constraints",
            "derive marginal-surface expansion from BQG metric observables",
            "derive horizon puncture/intersection map and fix area scale a_* internally",
            "derive horizon microstate Hilbert space and obtain coefficient 1/4 without calibration",
            "derive Bogoliubov map and 4D greybody transmission from BQG dynamics",
            "derive luminosity and backreaction from the same physical state",
            "derive radiation subsystem factorization/algebra and density matrix from BQG",
            "prove refinement and regulator independence of all physical observables",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
