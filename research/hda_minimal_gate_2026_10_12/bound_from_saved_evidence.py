#!/usr/bin/env python3
"""Projector lower bound from the saved Euclidean two-hit source-column-0 data.
This does NOT regenerate SU(2) amplitudes or test the full Lorentzian HDA.
Original: research/audit_2026_10_10/results.json
"""
from math import isclose

norm = 2.879453814704954
p_norm = 0.5503455636886702
q_norm = 2.826371177242493
upstream_q_fraction = 0.9634699130467531

assert isclose(norm**2, p_norm**2 + q_norm**2, rel_tol=0, abs_tol=1e-12)
assert isclose(q_norm**2 / norm**2, upstream_q_fraction, rel_tol=0, abs_tol=1e-12)
lower_bound = q_norm / norm
assert lower_bound > 0.98
print(f"PASS: any D with Dpsi inside original 32D code has relative defect >= {lower_bound:.12f}")
print("SCOPE: one source column; Euclidean H0,H1; not full Lorentzian BQG.")
