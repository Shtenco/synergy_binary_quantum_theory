from __future__ import annotations

import math
import numpy as np

# Frozen regulator-safe K5 commutator data already produced by the repository.
TOTAL_NORM = 1.681559985798016
FRACTIONS = np.array([
    0.55596684,
    0.18246774,
    0.15668486,
    0.08535889,
    0.01602488,
    0.00262260,
    0.00087420,
], dtype=float)

reported_graph_changed = 0.44403316047555935
if not __debug__:
    raise RuntimeError("Refusing optimized execution: this gate uses assert-based scientific checks")

assert abs(FRACTIONS.sum() - 1.0) < 2e-8
assert abs(FRACTIONS[1:].sum() - reported_graph_changed) < 2e-8
channel_norms = TOTAL_NORM * np.sqrt(FRACTIONS)
assert np.all(channel_norms > 0)

# Finite operator-algebra selftest of the graph-sector heat-kernel bridge.
# Basis states 0,1 belong to graph sector G; 2,3 to graph sector G'.
# Constraints are chosen so the common physical kernel mixes the two sectors.
# This verifies the structural formula, not the numerical BQG gravity projector.

# Physical kernel basis vectors: one same-graph superposition and one cross-graph superposition.
v0 = np.array([1, 1, 0, 0], dtype=complex) / math.sqrt(2)
v1 = np.array([0, 1, 1, 0], dtype=complex) / math.sqrt(2)
# Orthonormalize the intended physical subspace.
Q, _ = np.linalg.qr(np.column_stack([v0, v1]))
P_phys_expected = Q @ Q.conj().T
R = np.eye(4, dtype=complex) - P_phys_expected

# Two positive constraints sharing exactly the same kernel.
# C_a = R H_a R with H_a positive on the complement.
H1 = np.diag([2.0, 3.0, 4.0, 5.0]).astype(complex)
H2 = np.array([
    [3.0, .2, 0, 0],
    [.2, 4.0, .1, 0],
    [0, .1, 5.0, .3],
    [0, 0, .3, 6.0],
], dtype=complex)
C1 = R @ H1 @ R
C2 = R @ H2 @ R
constraints = [C1, C2]
G = np.array([[1.2, 0.15], [0.15, 0.9]], dtype=complex)

M = np.zeros((4, 4), dtype=complex)
for a, Ca in enumerate(constraints):
    for b, Cb in enumerate(constraints):
        M += Ca.conj().T @ (G[a, b] * Cb)
M = (M + M.conj().T) / 2

ev, U = np.linalg.eigh(M)
scale = max(float(np.max(np.abs(ev))), 1.0)
zero = ev < 1e-11 * scale
assert int(np.sum(zero)) == 2
P0 = (U[:, zero]) @ (U[:, zero]).conj().T
assert np.linalg.norm(P0 - P_phys_expected) < 1e-9

# Graph-sector projectors.
Pi_G = np.diag([1, 1, 0, 0]).astype(complex)
Pi_Gp = np.diag([0, 0, 1, 1]).astype(complex)

# Exact physical projector has a nonzero cross-graph block in this selftest.
cross_phys = Pi_Gp @ P0 @ Pi_G
assert np.linalg.norm(cross_phys) > 1e-6

# Heat kernel and semigroup / convergence checks.
def heat(T: float) -> np.ndarray:
    return (U * np.exp(-T * ev)) @ U.conj().T

T1, T2 = 0.37, 0.61
assert np.linalg.norm(heat(T1 + T2) - heat(T1) @ heat(T2)) < 1e-11

# Short-T derivative: (e^{-TM}-I)/T -> -M.
for T in (1e-4, 5e-5):
    err = np.linalg.norm((heat(T) - np.eye(4)) / T + M)
    assert err < 5e-3

positive = ev[~zero]
gap = float(np.min(positive))
for T in (1.0, 2.0, 4.0):
    obs = np.linalg.norm(heat(T) - P0, 2)
    pred = math.exp(-T * gap)
    assert abs(obs - pred) < 2e-10

# Microstate matrix-element identity for the master block:
# <beta|M|alpha> = sum_ab <C_a beta| G_ab |C_b alpha>.
alpha = np.array([1, 0, 0, 0], dtype=complex)
beta = np.array([0, 0, 1, 0], dtype=complex)
lhs = np.vdot(beta, M @ alpha)
rhs = 0j
for a, Ca in enumerate(constraints):
    for b, Cb in enumerate(constraints):
        rhs += G[a, b] * np.vdot(Ca @ beta, Cb @ alpha)
assert abs(lhs - rhs) < 1e-11

print("K5 graph-change support fraction =", FRACTIONS[1:].sum())
for n0, (f, a) in enumerate(zip(FRACTIONS, channel_norms)):
    print(f"n_zero={n0}: norm2_fraction={f:.8f} channel_norm={a:.12f}")
print("finite master gap =", gap)
print("physical cross-graph block norm (selftest) =", np.linalg.norm(cross_phys))
print("DERIVED: K_T(G',G)=Pi_G' exp(-T M) Pi_G")
print("DERIVED: A_phys=<G',beta|1_{0}(M)|G,alpha>")
print("IMPORTANT: frozen K5 44.4% is constraint-commutator support, not P_phys probability")
print("PASS")
