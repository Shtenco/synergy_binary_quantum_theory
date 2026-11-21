#!/usr/bin/env python3
"""Exact ten-pair reduction for K5 HDA. NOT a physical Lorentzian BQG calculation."""
from itertools import combinations
import numpy as np
import sympy as s

edges=list(combinations(range(5),2))
n=s.symbols("n0:5",real=True);m=s.symbols("m0:5",real=True)
w={(i,j):n[i]*m[j]-n[j]*m[i] for i,j in edges}
for i,j in edges:
    mid=((n[i]+n[j])*(m[j]-m[i])-(m[i]+m[j])*(n[j]-n[i]))/2
    assert s.expand(mid-w[i,j])==0
assert s.expand(sum(v*v for v in w.values())-
    (sum(v*v for v in n)*sum(v*v for v in m)-sum(n[i]*m[i] for i in range(5))**2))==0
for i,j,k,l in combinations(range(5),4):
    assert s.expand(w[i,j]*w[k,l]-w[i,k]*w[j,l]+w[i,l]*w[j,k])==0
# Every pointwise Leibniz derivation C^5 -> C^5 is zero:
# D(e_j)=D(e_j^2)=2 e_j D(e_j); all 25 matrix entries vanish.
assert s.diag(*[1-2*(i==j) for i in range(5) for j in range(5)]).rank()==25
a=s.symbols("a0:5");b=s.symbols("b0:5")
for i,j in edges:
    assert s.expand(a[j]*b[j]-a[i]*b[i]-
       ((a[i]+a[j])*(b[j]-b[i])+(b[i]+b[j])*(a[j]-a[i]))/2)==0
# Independent synthetic RHS verifies bilinear, lambda-polynomial and Gram identities.
rng=np.random.default_rng(20261012);d=6
def hermitian():
    t=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
    return (t+t.conj().T)/2
E=[hermitian() for _ in range(5)]
L=[hermitian() for _ in range(5)]
D={ij:hermitian() for ij in edges} # independent of all commutators
lam=1.37
H=[E[i]+lam*L[i] for i in range(5)]
comm=lambda x,y:x@y-y@x
R={ij:-1j*comm(H[ij[0]],H[ij[1]])-D[ij] for ij in edges}
V=np.stack([R[ij].reshape(-1) for ij in edges])
G=V.conj()@V.T
G=(G+G.conj().T)/2
emin,emax=np.linalg.eigvalsh(G)[[0,-1]]
assert emin>-1e-9
max_identity=0.;max_gram=0.
for t in range(250):
    N=rng.normal(size=5);M=rng.normal(size=5)
    omega=np.array([N[i]*M[j]-N[j]*M[i] for i,j in edges])
    lhs=-1j*comm(sum(N[i]*H[i] for i in range(5)),sum(M[i]*H[i] for i in range(5)))-sum(omega[z]*D[ij] for z,ij in enumerate(edges))
    rhs=sum(omega[z]*R[ij] for z,ij in enumerate(edges))
    max_identity=max(max_identity,np.linalg.norm(lhs-rhs,'fro')/max(np.linalg.norm(lhs,'fro'),1e-12))
    max_gram=max(max_gram,abs(np.linalg.norm(rhs,'fro')**2-np.real(omega@G@omega))/max(np.linalg.norm(rhs,'fro')**2,1e-12))
    assert np.linalg.norm(rhs,'fro') <= np.sqrt(emax)*np.linalg.norm(omega)+1e-10
for i,j in edges:
    polynomial=comm(E[i],E[j])+lam*(comm(E[i],L[j])+comm(L[i],E[j]))+lam**2*comm(L[i],L[j])
    assert np.linalg.norm(polynomial-comm(H[i],H[j]),'fro')<1e-10
assert max_identity<1e-12 and max_gram<1e-12
print({'status':'PASS (algebra and synthetic controls only)','ten_pairs':10,
'pointwise_derivation_dim':0,'midpoint_edge_leibniz':'exact',
'gram_eigenvalues_min_max':[float(emin),float(emax)],
'max_bilinear_relative_error':float(max_identity),'max_gram_relative_error':float(max_gram),
'physical_lorentzian_HDA_test':'NOT RUN: missing common-habitat independent microscopic D'})
