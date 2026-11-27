#!/usr/bin/env python3
"""Regression gate for the existing BQG reduced light cone and its IR null limit.

This gate deliberately distinguishes:
  (1) the exact finite-lattice band-front speed v_front(m),
  (2) the Lorentzian speed c_IR=3 in the joint critical scaling limit,
  (3) the graph-sector support cone, whose master parameter is not physical time.

NumPy is the only dependency.
"""
from __future__ import annotations

import math
import numpy as np


def energy(q: float, m: float) -> float:
    # q = k - pi/2, so cos(k)^2 = sin(q)^2.
    return math.sqrt(4.0 * m * m + 9.0 * math.sin(q) ** 2)


def group_speed(q: float, m: float) -> float:
    e = energy(q, m)
    if e == 0.0:
        return 3.0
    return abs(9.0 * math.sin(q) * math.cos(q) / e)


def lattice_front(m: float) -> float:
    return math.sqrt(4.0 * m * m + 9.0) - 2.0 * abs(m)


def exact_front_maximizer_x(m: float) -> float:
    # x = sin(q)^2; stationarity gives 9 x^2 + 8 m^2 x - 4 m^2 = 0.
    a = abs(m)
    return 2.0 * a * lattice_front(m) / 9.0


def continuum_energy(p: float, mu: float) -> float:
    return math.sqrt(4.0 * mu * mu + 9.0 * p * p)


def continuum_group_speed(p: float, mu: float) -> float:
    den = continuum_energy(p, mu)
    return 0.0 if den == 0.0 else 9.0 * abs(p) / den


def run() -> None:
    checks = 0

    # Exact finite-lattice front: analytic maximizer and closed-form speed agree.
    for m in (0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0):
        vf = lattice_front(m)
        x = exact_front_maximizer_x(m)
        assert 0.0 < x < 1.0
        qstar = math.asin(math.sqrt(x))
        vg = group_speed(qstar, m)
        assert abs(vg - vf) < 2e-12
        assert 0.0 < vf < 3.0
        checks += 1

    # Critical endpoint: the finite-lattice front reaches the IR null speed.
    assert lattice_front(0.0) == 3.0
    checks += 1

    # Joint scaling m=eps*mu, q=eps*p:
    # E/eps -> sqrt(4 mu^2 + 9 p^2), with O(eps^2) scaled-energy error.
    mus = (0.2, 0.7, 1.5)
    momenta = (0.3, 1.0, 2.0)
    epsilons = (0.2, 0.1, 0.05, 0.025)
    max_energy_errors = []
    max_speed_errors = []

    for eps in epsilons:
        e_errors = []
        v_errors = []
        for mu in mus:
            for p in momenta:
                m = eps * mu
                q = eps * p
                scaled_e = energy(q, m) / eps
                target_e = continuum_energy(p, mu)
                e_errors.append(abs(scaled_e - target_e))

                scaled_v = group_speed(q, m)
                target_v = continuum_group_speed(p, mu)
                v_errors.append(abs(scaled_v - target_v))
                checks += 2
        max_energy_errors.append(max(e_errors))
        max_speed_errors.append(max(v_errors))

    # Taylor errors are quadratic in eps; halving eps reduces them by ~4.
    for errors in (max_energy_errors, max_speed_errors):
        ratios = [errors[i] / errors[i + 1] for i in range(len(errors) - 1)]
        assert all(3.7 < ratio < 4.3 for ratio in ratios), ratios
        checks += len(ratios)

    # The limiting continuum mass shell is Lorentzian with c_IR=3:
    # omega^2 - 9 p^2 = (2 mu)^2.
    c_ir = 3.0
    for mu in mus:
        for p in momenta:
            omega = continuum_energy(p, mu)
            assert abs(omega * omega - c_ir * c_ir * p * p - 4.0 * mu * mu) < 1e-13
            checks += 1

    # Local 1+1 null vectors in ds^2 = -c_IR^2 dt^2 + dx^2.
    for sign in (-1.0, 1.0):
        k_t, k_x = 1.0, sign * c_ir
        null_norm = -(c_ir ** 2) * k_t * k_t + k_x * k_x
        assert abs(null_norm) < 1e-14
        checks += 1

    # The finite-m lattice front approaches 3 in the critical scaling m=eps*mu.
    for mu in mus:
        vals = [lattice_front(eps * mu) for eps in epsilons]
        assert all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))
        assert abs(vals[-1] - 3.0) < abs(vals[0] - 3.0)
        checks += 2

    print("checks =", checks)
    print("finite-lattice front maximizer identity PASS")
    print("critical front c_IR = 3 PASS")
    print("joint (m,q)->0 Lorentzian mass-shell scaling PASS")
    print("scaled energy and group-speed errors are O(epsilon^2) PASS")
    print("local null-vector norm PASS")
    print("finite-m front -> critical null speed PASS")
    print("IMPORTANT: master-flow T is not identified with physical time")
    print("IMPORTANT: full BQG gravity/matter cone universality remains OPEN")
    print("PASS")


if __name__ == "__main__":
    run()
