from __future__ import annotations

import numpy as np
import sympy as sp

I2 = sp.eye(2)
X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
Z = sp.Matrix([[1, 0], [0, -1]])
I8 = sp.eye(8)


def kron3(a: sp.Matrix, b: sp.Matrix, c: sp.Matrix) -> sp.Matrix:
    return sp.kronecker_product(a, b, c)


def entropy_bits(rho: np.ndarray) -> float:
    vals = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    vals = vals[vals > 1e-12]
    return float(-np.sum(vals * np.log2(vals)))


def partial_trace_three(rho: np.ndarray, keep: tuple[int, ...]) -> np.ndarray:
    dims = [2, 2, 2]
    arr = rho.reshape(dims + dims)
    trace_out = [i for i in range(3) if i not in keep]
    for i in sorted(trace_out, reverse=True):
        ncur = arr.ndim // 2
        arr = np.trace(arr, axis1=i, axis2=i + ncur)
    dim = 2 ** len(keep)
    return arr.reshape((dim, dim))


# Two-node exact reduced gluing operator embedded on links AB and BC.
H_AB = sp.Rational(3, 2) * I8 - sp.Rational(3, 4) * (
    kron3(X, X, I2) + kron3(Z, Z, I2)
)
H_BC = sp.Rational(3, 2) * I8 - sp.Rational(3, 4) * (
    kron3(I2, X, X) + kron3(I2, Z, Z)
)
H3 = sp.simplify(H_AB + H_BC)

lam0 = 3 - 3 * sp.sqrt(2) / 2
lam2 = 3 + 3 * sp.sqrt(2) / 2
expected = {lam0: 2, sp.Integer(3): 4, lam2: 2}
if not __debug__:
    raise RuntimeError("Refusing optimized execution: this gate uses assert-based scientific checks")

assert H3.eigenvals() == expected
assert len((H3 - lam0 * I8).nullspace()) == 2

# Symmetry-neutral ground state = normalized projector onto the two-dimensional
# ground eigenspace. This avoids arbitrary choice inside the exact degeneracy.
Hn = np.array(H3.evalf(), dtype=complex)
w, v = np.linalg.eigh(Hn)
mask = np.isclose(w, float(sp.N(lam0)))
V0 = v[:, mask]
rho0 = (V0 @ V0.conj().T) / V0.shape[1]

# Single logical nodes are maximally mixed.
for i in range(3):
    ri = partial_trace_three(rho0, (i,))
    assert np.allclose(ri, np.eye(2) / 2, atol=1e-11)

# Pair spectra: nearest neighbors versus next-nearest neighbors.
rho_ab = partial_trace_three(rho0, (0, 1))
rho_bc = partial_trace_three(rho0, (1, 2))
rho_ac = partial_trace_three(rho0, (0, 2))

spec_ab = np.linalg.eigvalsh(rho_ab)
spec_bc = np.linalg.eigvalsh(rho_bc)
spec_ac = np.linalg.eigvalsh(rho_ac)

sqrt2 = np.sqrt(2.0)
expected_nn = np.array([(3 - 2 * sqrt2) / 8, 1 / 8, 1 / 8, (3 + 2 * sqrt2) / 8])
expected_nnn = np.array([0, 1 / 4, 1 / 4, 1 / 2])
assert np.allclose(spec_ab, expected_nn, atol=1e-11)
assert np.allclose(spec_bc, expected_nn, atol=1e-11)
assert np.allclose(spec_ac, expected_nnn, atol=1e-11)

I_ab = 2.0 - entropy_bits(rho_ab)
I_bc = 2.0 - entropy_bits(rho_bc)
I_ac = 2.0 - entropy_bits(rho_ac)
assert np.isclose(I_ab, I_bc, atol=1e-11)
assert I_ab > I_ac
assert np.isclose(I_ac, 0.5, atol=1e-11)

# Oriented-volume correlations Q=(sqrt(3)/4)Y.
q2 = sp.Rational(3, 16)
YY_ab = q2 * kron3(Y, Y, I2)
YY_bc = q2 * kron3(I2, Y, Y)
YY_ac = q2 * kron3(Y, I2, Y)

corr_ab = float(np.trace(rho0 @ np.array(YY_ab.evalf(), complex)).real)
corr_bc = float(np.trace(rho0 @ np.array(YY_bc.evalf(), complex)).real)
corr_ac = float(np.trace(rho0 @ np.array(YY_ac.evalf(), complex)).real)

assert np.isclose(corr_ab, -3 / 32, atol=1e-11)
assert np.isclose(corr_bc, -3 / 32, atol=1e-11)
assert np.isclose(corr_ac, 0.0, atol=1e-11)

print("H3 spectrum:")
print("  3 - 3*sqrt(2)/2 : degeneracy 2")
print("  3                 : degeneracy 4")
print("  3 + 3*sqrt(2)/2 : degeneracy 2")
print("ground degeneracy = 2")
print("nearest-neighbor pair spectrum =", spec_ab)
print("next-nearest pair spectrum =", spec_ac)
print("I(A:B) =", I_ab)
print("I(B:C) =", I_bc)
print("I(A:C) =", I_ac)
print("<Q_A Q_B> =", corr_ab)
print("<Q_B Q_C> =", corr_bc)
print("<Q_A Q_C> =", corr_ac)
print("PASS")
