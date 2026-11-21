#!/usr/bin/env python3
"""Exact reference algebra; NOT a calculation of BQG's physical TT."""
import sympy as s
Z=s.diag(1,-1);I2=s.eye(2);J=s.Matrix([[0,-1],[1,0]])
assert J**2==-I2
seq=[s.simplify(s.trace(J**(-t)*Z*J**t*Z)/2) for t in range(8)]
assert seq==[1,-1,1,-1,1,-1,1,-1]
assert seq[0]*seq[2]-seq[1]**2==0 and seq[1]<0
ZZ=s.kronecker_product(Z,Z);I4=s.eye(4)
C=I4-ZZ;M=C*C;P=(I4+ZZ)/2
t=s.symbols("t",nonnegative=True)
assert P*P==P and C*P==s.zeros(4)
assert P*((-t*M).exp())*P==P
O1=s.kronecker_product(Z,I2);O2=s.kronecker_product(I2,Z)
assert s.trace(P*O1*O2)/s.trace(P)==1
assert s.trace(P*O1)/s.trace(P)==0
assert s.trace(P*O2)/s.trace(P)==0
assert s.trace(I4*O1*O2)/4==0
T=s.Matrix([[s.Rational(7,10),0,s.Rational(1,10)],
            [0,s.Rational(1,2),0],
            [s.Rational(1,10),0,s.Rational(1,5)]])
V=s.Matrix([[1,0],[0,1],[0,0]])
C0=V.T*V;C1=V.T*T*V;C2=V.T*T*T*V
PV=V*C0.inv()*V.T
R=s.simplify(C2-C1*C0.inv()*C1)
assert R==s.diag(s.Rational(1,100),0)
assert R==V.T*T*(s.eye(3)-PV)*T*V
assert T.eigenvals()=={
 (s.Rational(9,10)-s.sqrt(29)/10)/2:1,
 (s.Rational(9,10)+s.sqrt(29)/10)/2:1,
 s.Rational(1,2):1
}
assert 0<(s.sqrt(29)-5)/20<s.Rational(1,10)
print("PASS: master heat time no-go, C8, connected toy, TT 2x2 Schur and Ritz bound")
print("PHYSICAL BQG: full Lorentzian HDA BLOCKED; physical TT pole BLOCKED")
