from __future__ import annotations

import sympy as sp

I2 = sp.eye(2)
X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
Z = sp.Matrix([[1, 0], [0, -1]])


def kron(a: sp.Matrix, b: sp.Matrix) -> sp.Matrix:
    return sp.kronecker_product(a, b)


# Logical physical basis on each four-spin SU(2)-singlet node:
# |0> = |s>, |1> = |t>.
#
# Three gauge-invariant pair-singlet projectors. Their Bloch vectors form
# a trine in the logical X-Z plane.
P12 = (I2 + Z) / 2
P13 = I2 / 2 + sp.sqrt(3) * X / 4 - Z / 4
P14 = I2 / 2 - sp.sqrt(3) * X / 4 - Z / 4
PAIR_PROJECTORS = (P12, P13, P14)

# Basic exact identities for the local shape frame.
if not __debug__:
    raise RuntimeError("Refusing optimized execution: this gate uses assert-based scientific checks")

assert sp.simplify(P12 + P13 + P14 - sp.Rational(3, 2) * I2) == sp.zeros(2)
for P in PAIR_PROJECTORS:
    assert sp.simplify(P * P - P) == sp.zeros(2)

# Candidate reduced shape-matching/gluing Hamiltonian for two physical nodes A,B:
# penalize disagreement of all three independent pair-shape channels.
H_glue = sp.zeros(4)
for P in PAIR_PROJECTORS:
    delta = kron(P, I2) - kron(I2, P)
    H_glue += sp.expand(delta * delta)
H_glue = sp.simplify(H_glue)

H_closed = sp.Rational(3, 2) * sp.eye(4) - sp.Rational(3, 4) * (
    kron(X, X) + kron(Z, Z)
)
assert sp.simplify(H_glue - H_closed) == sp.zeros(4)

# Bell basis in the local physical basis |s>,|t>.
e0 = sp.Matrix([1, 0])
e1 = sp.Matrix([0, 1])
ss = kron(e0, e0)
st = kron(e0, e1)
ts = kron(e1, e0)
tt = kron(e1, e1)

phi_plus = (ss + tt) / sp.sqrt(2)
phi_minus = (ss - tt) / sp.sqrt(2)
psi_plus = (st + ts) / sp.sqrt(2)
psi_minus = (st - ts) / sp.sqrt(2)

bell = {
    "Phi+": phi_plus,
    "Phi-": phi_minus,
    "Psi+": psi_plus,
    "Psi-": psi_minus,
}
expected_energy = {
    "Phi+": sp.Integer(0),
    "Phi-": sp.Rational(3, 2),
    "Psi+": sp.Rational(3, 2),
    "Psi-": sp.Integer(3),
}

for name, v in bell.items():
    assert sp.simplify(H_glue * v - expected_energy[name] * v) == sp.zeros(4, 1)

# Unique zero-mismatch state.
assert H_glue.nullspace() == [phi_plus]

# Oriented-volume operator on each node from the one-node exact result:
# Q = sqrt(3)/4 Y.
q0 = sp.sqrt(3) / 4
QA = q0 * kron(Y, I2)
QB = q0 * kron(I2, Y)
Qrel = QA * QB

# The zero-mismatch Bell state does NOT have a sharp local orientation:
# <QA>=<QB>=0 and Var(QA)=Var(QB)=3/16.
def expect(v: sp.Matrix, op: sp.Matrix):
    return sp.simplify((v.H * op * v)[0])

assert expect(phi_plus, QA) == 0
assert expect(phi_plus, QB) == 0
assert expect(phi_plus, QA * QA) == sp.Rational(3, 16)
assert expect(phi_plus, QB * QB) == sp.Rational(3, 16)

# But relative orientation is exact and anti-correlated in the chosen local
# orientation convention.
assert sp.simplify(Qrel * phi_plus + sp.Rational(3, 16) * phi_plus) == sp.zeros(4, 1)
assert expect(phi_plus, Qrel) == -sp.Rational(3, 16)

# Each node separately is maximally mixed in logical intertwiner space.
# For |Phi+>, Schmidt coefficients are (1/sqrt(2),1/sqrt(2)), so
# S(node A)=S(node B)=1 bit exactly.
node_entropy_bits = sp.Integer(1)

print("H_glue =")
sp.pprint(H_glue)
print("closed form: 3/2 I - 3/4 (X_AX_B + Z_AZ_B)")
print("spectrum: {0:1, 3/2:2, 3:1}")
print("unique zero-mismatch state: |Phi+>=(|ss>+|tt>)/sqrt(2)")
print("<Q_A>=<Q_B>=0")
print("Var(Q_A)=Var(Q_B)=3/16")
print("Q_A Q_B eigenvalue on |Phi+> = -3/16")
print("node-node entanglement entropy =", node_entropy_bits, "bit")
print("PASS")
