from __future__ import annotations

import sympy as sp

I = sp.I
SQRT2 = sp.sqrt(2)
SQRT3 = sp.sqrt(3)

sx = sp.Matrix([[0, 1], [1, 0]]) / 2
sy = sp.Matrix([[0, -I], [I, 0]]) / 2
sz = sp.Matrix([[1, 0], [0, -1]]) / 2
ID = sp.eye(2)
J = [sx, sy, sz]


def kron4(a, b, c, d):
    return sp.kronecker_product(a, b, c, d)


def ket(bits: str) -> sp.Matrix:
    v = sp.zeros(16, 1)
    v[int(bits, 2), 0] = 1
    return v


def jop(site: int, component: int) -> sp.Matrix:
    ops = [ID, ID, ID, ID]
    ops[site] = J[component]
    return kron4(*ops)


s = (ket("0101") - ket("0110") - ket("1001") + ket("1010")) / 2

t = (
    ket("0011")
    + ket("1100")
    - (ket("0101") + ket("0110") + ket("1001") + ket("1010")) / 2
) / SQRT3

if not __debug__:
    raise RuntimeError("Refusing optimized execution: this gate uses assert-based scientific checks")

assert sp.simplify((s.H * s)[0]) == 1
assert sp.simplify((t.H * t)[0]) == 1
assert sp.simplify((s.H * t)[0]) == 0

Q = sp.zeros(16)
for a in range(3):
    for b in range(3):
        for c in range(3):
            eps = sp.LeviCivita(a, b, c)
            if eps:
                Q += eps * jop(0, a) * jop(1, b) * jop(2, c)

D12 = sum((jop(0, a) * jop(1, a) for a in range(3)), sp.zeros(16))
D23 = sum((jop(1, a) * jop(2, a) for a in range(3)), sp.zeros(16))
assert sp.simplify(Q - I * (D12 * D23 - D23 * D12)) == sp.zeros(16)

basis = [s, t]
Qphys = sp.Matrix([
    [sp.simplify((u.H * Q * v)[0]) for v in basis]
    for u in basis
])

sigma_y = sp.Matrix([[0, -I], [I, 0]])
assert sp.simplify(Qphys - (SQRT3 / 4) * sigma_y) == sp.zeros(2)
assert sp.simplify(Qphys**2 - sp.Rational(3, 16) * sp.eye(2)) == sp.zeros(2)

psi_plus = (s + I * t) / SQRT2
psi_minus = (s - I * t) / SQRT2
assert sp.simplify(Q * psi_plus - (SQRT3 / 4) * psi_plus) == sp.zeros(16, 1)
assert sp.simplify(Q * psi_minus + (SQRT3 / 4) * psi_minus) == sp.zeros(16, 1)

alpha, phi = sp.symbols("alpha phi", real=True)
p12 = sp.cos(alpha) ** 2
p13 = sp.Rational(1, 2) - sp.cos(2 * alpha) / 4 + SQRT3 * sp.sin(2 * alpha) * sp.cos(phi) / 4
p14 = sp.Rational(1, 2) - sp.cos(2 * alpha) / 4 - SQRT3 * sp.sin(2 * alpha) * sp.cos(phi) / 4

Ap = sp.simplify(2 * sum(
    (p - sp.Rational(1, 2))**2 for p in (p12, p13, p14)
))
Qmean = SQRT3 * sp.sin(2 * alpha) * sp.sin(phi) / 4
Qvar = sp.Rational(3, 16) - Qmean**2
assert sp.trigsimp(Ap - 4 * Qvar) == 0

eps, delta = sp.symbols("eps delta", real=True)
quad_anisotropy = 3 * eps**2 + sp.Rational(3, 4) * delta**2

print("Q = i [J1.J2, J2.J3] : EXACT")
print("Q_phys =")
sp.pprint(Qphys)
print("Q_phys^2 = 3/16 I : EXACT")
print("eigenvalues:", [-SQRT3 / 4, SQRT3 / 4])
print("isotropic states: (|s> +/- i|t>)/sqrt(2)")
print("pair singlet weights:")
print("p12=p34 =", sp.simplify(p12))
print("p13=p24 =", sp.simplify(p13))
print("p14=p23 =", sp.simplify(p14))
print("isotropy locus: alpha=pi/4, phi=pi/2 or 3pi/2")
print("A_p =", sp.trigsimp(Ap))
print("<Q> =", Qmean)
print("A_p = 4 Var(Q) : EXACT")
print("local quadratic A_p =", quad_anisotropy)
print("PASS")
