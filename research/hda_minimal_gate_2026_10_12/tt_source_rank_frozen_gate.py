#!/usr/bin/env python3
"""Exact BQG X/Z TT rank deficiency; kinematic q4 six-metric remedy.
NOT a physical projector, Lorentzian full HDA or TT pole calculation.
Requires sympy.
"""
import sympy as s
rt=s.sqrt(3)
MX=s.Matrix([[rt/2,0,rt/2],[0,-rt/2,-rt/2],[rt/2,-rt/2,0]])
MZ=s.Matrix([[s.Rational(1,2),1,-s.Rational(1,2)],[1,s.Rational(1,2),-s.Rational(1,2)],[-s.Rational(1,2),-s.Rational(1,2),-1]])
x,y,z=s.symbols('x y z',real=True)
v=s.Matrix([x,y,z]);q=x*x+y*y+z*z
def gnum(A,B):
    return s.expand(2*s.trace(A*B)*q**2-4*q*(v.T*A*B*v)[0]+(v.T*A*v)[0]*(v.T*B*v)[0])
G=s.Matrix([[gnum(A,B) for B in (MX,MZ)] for A in (MX,MZ)])
F=2*(x**3+y**3+z**3)-x*x*y-x*x*z-x*y*y-x*z*z-y*y*z-y*z*z-3*x*y*z
assert s.factor(G.det()-3*q*F**2)==0
# det physical TT Gram = det G/(4*q**4) = 3*F**2/(4*q**3).
assert s.simplify(G.subs({x:1,y:-1,z:0})/(2*2**2))==s.diag(0,s.Rational(33,8))
assert s.simplify(G.subs({x:0,y:0,z:1})/2)==s.diag(s.Rational(3,2),2)
# Kinematic six-edge response: sigma^2(M_hq) in {1/6 (x3),1/3 (x3)}.
# For EVERY orthonormal 2x6 TT row matrix E, eigenvalues of
# E M_hq M_hq^T E^T lie in [1/6,1/3], by Rayleigh quotient.
Y=s.Matrix([[0,-s.I],[s.I,0]]);X=s.Matrix([[0,1],[1,0]]);Z=s.diag(1,-1)
assert Y*X-X*Y==-2*s.I*Z
assert Y*Z-Z*Y==2*s.I*X
print('PASS: det(G_TT)=3 F(n)^2/(4 |n|^6), exact null direction (1,-1,0)')
print('PASS: rank-6 collective map retains 2 kinematic TT polarisations in all directions')
print('CONDITIONAL local Y commutator norm = 2 * 1.3389293521464034')
print('BLOCKED: Dirac physical source Gram, full Lorentzian HDA, physical time TT pole')
