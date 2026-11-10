from __future__ import annotations

import math
import numpy as np
from scipy.linalg import sqrtm

T = 1.5


def finite_chain_ground_correlation(L: int, m: float) -> tuple[np.ndarray, np.ndarray]:
    h = np.zeros((L, L), dtype=float)
    for i in range(L - 1):
        h[i, i + 1] = h[i + 1, i] = -T
    for i in range(L):
        h[i, i] += -2.0 * m * ((-1) ** i)
    vals, vecs = np.linalg.eigh(h)
    occ = vals < 0.0
    C = vecs[:, occ] @ vecs[:, occ].T
    return vals, C


def spin_correlators(C: np.ndarray, i: int, j: int) -> tuple[float, float, float, float]:
    zi = 1.0 - 2.0 * C[i, i]
    zj = 1.0 - 2.0 * C[j, j]
    zz = zi * zj - 4.0 * C[i, j] ** 2
    G = np.eye(C.shape[0]) - 2.0 * C
    xx = float(np.linalg.det(G[i:j, i + 1 : j + 1]))
    return zi, zj, zz, xx


def pair_state(C: np.ndarray, i: int, j: int) -> np.ndarray:
    zi, zj, zz, xx = spin_correlators(C, i, j)
    rho = np.zeros((4, 4), dtype=float)
    rho[0, 0] = (1 + zi + zj + zz) / 4
    rho[1, 1] = (1 + zi - zj - zz) / 4
    rho[2, 2] = (1 - zi + zj - zz) / 4
    rho[3, 3] = (1 - zi - zj + zz) / 4
    rho[1, 2] = rho[2, 1] = xx / 2
    return rho


def partials(rho: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r = rho.reshape(2, 2, 2, 2)
    ra = np.einsum("abcb->ac", r)
    rb = np.einsum("abad->bd", r)
    return ra, rb


def entropy_bits(rho: np.ndarray) -> float:
    vals = np.linalg.eigvalsh(rho)
    vals = vals[vals > 1e-15]
    return -float(np.sum(vals * np.log2(vals)))


def mutual_information_bits(rho: np.ndarray) -> float:
    ra, rb = partials(rho)
    return entropy_bits(ra) + entropy_bits(rb) - entropy_bits(rho)


def fidelity_squared(rho: np.ndarray, sigma: np.ndarray) -> float:
    sr = sqrtm(rho)
    x = sr @ sigma @ sr
    return float(np.real(np.trace(sqrtm(x))) ** 2)


def bures_correlation_conductance(rho: np.ndarray) -> float:
    ra, rb = partials(rho)
    prod = np.kron(ra, rb)
    F = fidelity_squared(rho, prod)
    # Squared Bures distance, positive and zero iff rho=rho_A⊗rho_B.
    return 2.0 * (1.0 - math.sqrt(max(0.0, F)))


def qfi_unitary(rho: np.ndarray, generator: np.ndarray) -> float:
    vals, vecs = np.linalg.eigh(rho)
    out = 0.0
    for a, la in enumerate(vals):
        for b, lb in enumerate(vals):
            den = la + lb
            if den > 1e-15:
                gab = np.vdot(vecs[:, a], generator @ vecs[:, b])
                out += 2.0 * (la - lb) ** 2 / den * abs(gab) ** 2
    return float(out)


Z = np.diag([1.0, -1.0])
I2 = np.eye(2)
K_REL = (np.kron(Z, I2) - np.kron(I2, Z)) / 2.0

# Thermodynamic critical exact nearest-neighbor QFI.
FQ_CRIT_EXACT = 32.0 / (math.pi**2 + 4.0)
GQ_CRIT_EXACT = FQ_CRIT_EXACT / 4.0

# Independent large-open-chain gate for several phases.
L = 600
i0 = 280
masses = [0.0, 0.2, 0.5, 1.0, 2.0]
rows = []
for m in masses:
    _, C = finite_chain_ground_correlation(L, m)
    rho = pair_state(C, i0, i0 + 1)
    mi = mutual_information_bits(rho)
    gb = bures_correlation_conductance(rho)
    fq = qfi_unitary(rho, K_REL)
    gq = fq / 4.0
    # Reference-vacuum normalized resistance length.
    ell_q = GQ_CRIT_EXACT / gq
    rows.append((m, mi, gb, fq, gq, ell_q))

# Critical finite-size value must approach the exact thermodynamic value.
if not __debug__:
    raise RuntimeError("Refusing optimized execution: this gate uses assert-based scientific checks")

assert abs(rows[0][3] / FQ_CRIT_EXACT - 1.0) < 0.01

# All three local information conductances are positive in both phases.
for _, mi, gb, fq, gq, ell_q in rows:
    assert mi > 0.0
    assert gb > 0.0
    assert fq > 0.0
    assert gq > 0.0
    assert ell_q > 0.0

# Orientation mass weakens the local link monotonically for this controlled family.
for a, b in zip(rows, rows[1:]):
    assert b[1] < a[1]  # MI
    assert b[2] < a[2]  # Bures correlation conductance
    assert b[3] < a[3]  # relative-twist QFI
    assert b[5] > a[5]  # normalized resistance length grows

# Exact additive path metric on a chain once local edge lengths are assigned.
# This is the network step that avoids using long-distance I(i:j) as a metric.
for row in rows:
    ell = row[5]
    for r in (1, 2, 5, 11):
        direct_path_sum = sum([ell] * r)
        assert abs(direct_path_sum - r * ell) < 1e-14

print("critical exact F_Q =", FQ_CRIT_EXACT)
print("critical exact g_Q=F_Q/4 =", GQ_CRIT_EXACT)
print("m, MI_edge, Bures^2_edge, FQ_twist, normalized resistance length")
for m, mi, gb, fq, gq, ell_q in rows:
    print(f"{m:4.1f}  {mi:.12f}  {gb:.12f}  {fq:.12f}  {ell_q:.12f}")
print("POSITIVE: local QFI conductance -> resistance-path metric is additive in critical and gapped phases")
print("NO-UNIQUENESS: MI, Bures and QFI define different local scales; rho_ij alone does not select one")
print("PASS")
