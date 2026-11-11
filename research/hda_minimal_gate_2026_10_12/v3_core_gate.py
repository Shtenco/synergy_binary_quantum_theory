#!/usr/bin/env python3
"""K5 nodal lapse aliasing and saved Euclidean bound. Full Lorentzian HDA NOT run."""
import itertools, json
import numpy as np
import sympy as s

x,y,z=s.symbols('x y z',real=True)
coords=(x,y,z)
pts=[(0,0,0),(1,0,0),(0,1,0),(0,0,1),(s.Rational(1,3),s.Rational(1,4),s.Rational(1,5))]
def r2(u,v):return sum((a-b)**2 for a,b in zip(u,v))
N=x*s.prod(r2(coords,pts[j])/r2(pts[0],pts[j]) for j in range(1,5))
vals=[s.simplify(N.subs(dict(zip(coords,p)))) for p in pts]
grads=[[s.simplify(s.diff(N,t).subs(dict(zip(coords,p)))) for t in coords] for p in pts]
assert vals==[0]*5 and grads==[[1,0,0]]+[[0,0,0]]*4
assert -grads[0][0]==-1
edges=list(itertools.combinations(range(5),2))
n=s.symbols('n0:5');m=s.symbols('m0:5')
for i,j in edges:
    w=n[i]*m[j]-n[j]*m[i]
    mid=((n[i]+n[j])*(m[j]-m[i])-(m[i]+m[j])*(n[j]-n[i]))/2
    assert s.expand(w-mid)==0
# From research/audit_2026_10_10/results.json, hbar=1 and conditional Dpsi in P32.
full=2.879453814704954;inside=0.5503455636886702;outside=2.826371177242493
assert abs(full**2-inside**2-outside**2)<1e-12
bound=outside**2
rng=np.random.default_rng(20261012)
r=rng.normal(size=(10,7))+1j*rng.normal(size=(10,7))
g=r.conj()@r.T
eig=np.linalg.eigvalsh(g)
mx=np.max(np.sum(abs(r)**2,axis=1))
assert mx-1e-10<=eig[-1]<=10*mx+1e-10
assert abs(2.0**(-4*6)-1/16777216)<1e-20
print(json.dumps({
  'nodal_lapse_values':list(map(str,vals)),
  'nodal_gradients':[[str(t) for t in p] for p in grads],
  'independent_D_probe':-1,
  'euclidean_lambda_max_bound_conditional':bound,
  'euclidean_relative_squared_bound_conditional':bound/full**2,
  'gram_sandwich_synthetic':True,
  'scale_trap_at_n6':2.0**(-24),
  'physical_HDA':'BLOCKED — independent Lorentzian D on common habitat missing'
},indent=2))
