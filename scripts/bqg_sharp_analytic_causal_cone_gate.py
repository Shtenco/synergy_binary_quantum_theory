#!/usr/bin/env python3
"""Regression gate for the sharp analytic BQG propagation cone."""
from __future__ import annotations

import math
import numpy as np


def vmax(m: float) -> float:
    # Rationalized form avoids catastrophic cancellation for large |m|.
    a = 2.0 * abs(m)
    return 9.0 / (math.hypot(a, 3.0) + a)


def mu_c(m: float) -> float:
    return math.asinh(2*abs(m)/3.0)


def imag_energy_sup_grid(m: float, mu: float, n: int = 200001) -> float:
    k = np.linspace(-math.pi, math.pi, n)
    z = 4*m*m + 9*np.cos(k + 1j*mu)**2
    e = np.sqrt(z)
    return float(np.max(np.abs(e.imag)))


def imag_energy_sup_exact(m: float, mu: float) -> float:
    return 0.5 * vmax(m) * math.sinh(2*mu)


def critical_cone_velocity(mu: float) -> float:
    if mu == 0.0:
        return 3.0
    return 3.0 * math.sinh(mu) / mu


def gapped_cone_velocity(m: float, mu: float) -> float:
    if mu == 0.0:
        return vmax(m)
    return vmax(m) * math.sinh(2*mu) / (2*mu)


def run():
    masses = [0.2, 0.5, 1.0, 2.0, 5.0]
    checks = 0

    for m in masses:
        muc = mu_c(m)
        for frac in (0.1, 0.3, 0.5, 0.7, 0.9):
            mu = frac * muc
            num = imag_energy_sup_grid(m, mu)
            exact = imag_energy_sup_exact(m, mu)
            assert abs(num - exact) < 5e-7 * max(1.0, abs(exact)), (
                m, frac, num, exact
            )

            vm = vmax(m)
            vmu = gapped_cone_velocity(m, mu)
            assert vmu > vm
            checks += 1

        xi = 1.0/muc
        assert abs(vmax(m)/3.0 - math.exp(-1.0/xi)) < 1e-13

    for mu in (0.01, 0.05, 0.1, 0.2):
        assert critical_cone_velocity(mu) > 3.0
    assert critical_cone_velocity(0.0) == 3.0
    assert abs(critical_cone_velocity(1e-6) - 3.0) < 1e-9
    for m in masses:
        assert gapped_cone_velocity(m, 0.0) == vmax(m)
    assert vmax(1e16) > 0.0

    print("gapped contour checks =", checks)
    print("critical front =", 3.0)
    for m in masses:
        print(
            f"m={m:g}",
            "mu_c=", mu_c(m),
            "xi=", 1.0/mu_c(m),
            "vmax=", vmax(m),
        )
    print("IMAGINARY-ENERGY SUPREMUM IDENTITY PASS")
    print("ANALYTIC STRIP == INVERSE CORRELATION LENGTH PASS")
    print("SHARP EXPONENTIAL CONE -> VMAX PASS")
    print("PASS")


if __name__ == "__main__":
    run()
