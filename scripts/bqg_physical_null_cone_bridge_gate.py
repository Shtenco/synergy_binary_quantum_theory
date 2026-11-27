#!/usr/bin/env python3
"""Fail-closed algebraic checks for the BQG physical null-cone bridge.

This verifies identities already stated in repository documents. It does not
reconstruct a physical Lorentzian metric or claim a full BQG light cone.
"""
from __future__ import annotations

import math
import numpy as np


def null_coordinate_velocity(q: np.ndarray, lapse: float, shift: np.ndarray,
                             direction: np.ndarray) -> np.ndarray:
    """Return u=-shift+N*s, with q(s,s)=1."""
    norm2 = float(direction @ q @ direction)
    if not np.isfinite(norm2) or norm2 <= 0:
        raise ValueError("direction must have positive q-norm")
    s = direction / math.sqrt(norm2)
    return -shift + lapse * s


def check_adm_null_identity() -> int:
    rng = np.random.default_rng(20261127)
    checks = 0
    for _ in range(500):
        X = rng.normal(size=(3, 3))
        q = X.T @ X + 0.2 * np.eye(3)
        lapse = float(rng.uniform(0.1, 3.0))
        shift = rng.normal(size=3)
        direction = rng.normal(size=3)
        u = null_coordinate_velocity(q, lapse, shift, direction)
        residual = float((u + shift) @ q @ (u + shift) - lapse**2)
        assert abs(residual) < 2e-12 * max(1.0, lapse**2)
        checks += 1
    return checks


def check_hda_tt_speed() -> int:
    # Hamiltonian equations imply c_T^2=A*B. Standard HDA has AB=1.
    checks = 0
    for A in (0.125, 0.5, 1.0, 2.0, 8.0, 17.0):
        B = 1.0 / A
        assert abs(A * B - 1.0) < 2e-15
        checks += 1
    for A, B in ((2.0, 1.0), (1.0, 2.0), (0.4, 0.7)):
        assert abs(A * B - 1.0) > 1e-6
        checks += 1
    return checks


def spectral_front(m: float) -> float:
    return math.sqrt(4.0 * m * m + 9.0) - 2.0 * abs(m)


def check_reduced_spectral_identity() -> int:
    checks = 0
    for m in (0.0, 1e-8, 0.02, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0):
        mu = math.asinh(2.0 * abs(m) / 3.0)
        assert abs(spectral_front(m) / 3.0 - math.exp(-mu)) < 3e-14
        assert spectral_front(m) > 0.0
        checks += 1
    # The reduced model does not itself fix its speed to the HDA value 1.
    assert abs(spectral_front(0.0) - 1.0) > 1.0
    assert abs(spectral_front(0.5) - 1.0) > 0.5
    return checks + 2


def check_graph_support_order() -> int:
    checks = 0
    for distance in range(1001):
        n = (distance + 1) // 2
        assert distance <= 2 * n
        if n:
            assert distance > 2 * (n - 1)
        checks += 1
    return checks


def main() -> None:
    counts = {
        "ADM null-coordinate identity": check_adm_null_identity(),
        "HDA/TT speed normalization": check_hda_tt_speed(),
        "reduced spectral-front identity": check_reduced_spectral_identity(),
        "graph support-order ceiling": check_graph_support_order(),
    }
    for name, count in counts.items():
        print(f"PASS {name}: {count} checks")
    print(f"TOTAL PASS: {sum(counts.values())} checks")
    print("SCOPE: algebraic bridge regression only; no physical BQG null cone derived")


if __name__ == "__main__":
    main()
