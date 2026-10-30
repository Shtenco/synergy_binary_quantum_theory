#!/usr/bin/env python3
"""Exact SO(5) q2 transporter carry gate. NOT the full Lorentzian BQG HDA/TT.
Reuses M(a,4) defined in main/scripts/su2_quantum_link_vector5_gate.py.
"""
import numpy as np
from scipy.linalg import expm

def M(a,b):
    if a>b:return -M(b,a)
    mat=np.zeros((5,5),complex)
    mat[a,b]=-1j;mat[b,a]=1j
    return mat
tol=1e-12
Pg=np.diag([1,1,1,1,0]).astype(complex)
P0=np.diag([0,0,0,0,1]).astype(complex)
U=[M(a,4) for a in range(4)]
ident=np.eye(5,dtype=complex)
for a in range(4):
    for b in range(4):
        target=np.outer(ident[a],ident[b])
        assert np.linalg.norm(Pg@U[a]@P0@U[b]@Pg-target)<tol

cycle=[0,1,3,2]  # Gray 00->01->11->10
S=np.zeros((4,4),complex)
for k in range(4):
    S[cycle[(k+1)%4],cycle[k]]=1
assert np.linalg.matrix_rank(S)==4
w=np.array([1,2,3,4],float)/np.sqrt(30)
H=sum(w[a]*U[a] for a in range(4))
L=(Pg@H@P0@H@Pg)[:4,:4]
assert np.linalg.matrix_rank(L)==1
assert np.allclose(L,np.outer(w,w))
for t in (0,.13,.67,1.7,5):
    direct=(Pg@expm(-1j*t*H)@Pg)[:4,:4]
    analytic=np.eye(4)+(np.cos(t)-1)*np.outer(w,w)
    assert np.linalg.norm(direct-analytic)<tol

Kch=np.zeros((5,5),complex)
Ks=np.zeros((5,5),complex)
for k in range(4):
    a,b=cycle[k],cycle[(k+1)%4]
    forward=Pg@U[b]@P0@U[a]@Pg
    backward=Pg@U[a]@P0@U[b]@Pg
    Kch+=(forward-backward)/(2j)
    Ks+=(forward+backward)/2
Kch=Kch[:4,:4];Ks=Ks[:4,:4]
assert np.linalg.norm(Kch-(S-S.conj().T)/(2j))<tol
assert np.linalg.norm(Ks-(S+S.conj().T)/2)<tol
assert np.linalg.norm(S-Ks-1j*Kch)<tol
assert np.linalg.norm(Kch@Ks)<tol
assert np.linalg.norm(Ks@Ks+Kch@Kch-np.eye(4))<tol
assert np.allclose(np.linalg.eigvalsh(Kch),[-1,0,0,1])
assert np.linalg.matrix_rank(Kch)==2

# S is unitary with four singular values 1.
sval=np.linalg.svd(S,compute_uv=False)
assert np.allclose(sval,np.ones(4))
rank_one_error=float(np.linalg.norm(sval[1:])/np.linalg.norm(sval))
assert abs(rank_one_error-np.sqrt(3)/2)<tol
# Positive second-order Hermitian return also cannot equal directed S.
pos=(S+S.conj().T)/2
e,v=np.linalg.eigh(pos)
closest_psd=v@np.diag(np.maximum(0,e))@v.conj().T
psd_error=float(np.linalg.norm(S-closest_psd,'fro')/np.linalg.norm(S,'fro'))
assert abs(psd_error-np.sqrt(3)/2)<tol

# Four distinct event channels are sufficient for an engineered carry lift.
T=np.zeros((8,8),complex)
for k in range(4):
    T[4+k,cycle[k]]=1
    T[cycle[(k+1)%4],4+k]=1
assert np.linalg.norm(T.conj().T@T-np.eye(8))<tol
assert np.linalg.norm(np.linalg.matrix_power(T,2)[:4,:4]-S)<tol
assert np.linalg.norm(np.linalg.matrix_power(T,8)-np.eye(8))<tol
print("PASS: micro unit factorization, rank1 no-go, exact chiral operator, engineered C8")
print("relative best rank1/PSD defect:",rank_one_error,psd_error)
print("K_chi spectrum:",np.linalg.eigvalsh(Kch))
print("PHYSICAL_BQG_HDA=OPEN PHYSICAL_BQG_TT=OPEN PHYSICAL_CARRY=OPEN")
