#!/usr/bin/env python3
"""Exact low-cost physical TT moment theorem: reduced/free benchmark ONLY.
No microscopic Lorentzian HDA, physical rigging limit or interacting BQG TT pole is claimed.
"""
from __future__ import annotations
import math,json
import numpy as np
import sympy as s
from itertools import product

z,q=s.symbols("z q",positive=True)
w=s.symbols("w",real=True)
mass2=(1-z)**2/z
c0=z/(q*(1-z*z))
K=((1-z)**2+4*z*s.sin(w/2)**2)/(c0*(1-z*z))
assert s.simplify(K-q*(4*s.sin(w/2)**2+mass2))==0
assert s.simplify(z/(c0*(1-z*z))-q)==0
A,B,r,t=s.symbols("A B r t",real=True)
C0=A+B;C1=A*r+B*t;C2=A*r*r+B*t*t
assert s.factor(C0*C2-C1*C1)==A*B*(r-t)**2
H1=s.Matrix([[0,0,0],[0,0,1],[0,1,0]])
H2=s.diag(0,1,-1);P=s.diag(1,0,0)
assert H1*H2!=H2*H1 and (H1*H2-H2*H1)*P==s.zeros(3)
assert P*s.diag(0,1,2)*P==s.zeros(3)

R=1/math.sqrt(3);Z=1.7
cases=[]
for qmag,axis in [(0.015,[1,0,0]),(0.02,[1,2,3]),(0.055,[1,1,0])]:
    d=np.array(axis,dtype=float);d/=np.linalg.norm(d)
    om2=R*R*sum((2*np.sin(qmag*d/2))**2)
    om=math.sqrt(om2)
    rate=math.exp(-2*math.asinh(om/2))
    c=1/(Z*om*math.sqrt(om2+4))
    c1=c*rate;c2=c1*rate
    rate_hat=c1/c
    z_hat=rate_hat/(c*(1-rate_hat**2))
    om2_hat=(1-rate_hat)**2/rate_hat
    assert abs(z_hat/Z-1)<1e-12 and abs(om2_hat/om2-1)<1e-12
    assert abs(c*c2-c1*c1)<1e-10
    cases.append({"q":qmag,"C0":c,"C1":c1,"C2":c2,
                  "Z_recovered":z_hat,"Omega2_recovered":om2_hat,
                  "E_over_q":-math.log(rate_hat)/qmag})
w1,w2=1.2,.7
r1=math.exp(-R*.035);r2=math.exp(-1.4)
def det(weights):
    moments=[sum(weight*rate**n for weight,rate in zip(weights,(r1,r2))) for n in (0,1,2)]
    return moments[0]*moments[2]-moments[1]**2
assert det((w1,w2))>0 and det((w1,-w2))<0
result={"symbolic":"PASS", "reduced_free_TT":"PASS","physical_BQG_HDA":"NOT_RUN",
"physical_BQG_TT":"NOT_RUN", "three_lag_hankel_positive_mixture":det((w1,w2)),
"three_lag_ghost_negative_control":det((w1,-w2)),
"free_cases":cases,"status":"mathematical theorem and benchmark only"}
print(json.dumps(result,indent=2))
