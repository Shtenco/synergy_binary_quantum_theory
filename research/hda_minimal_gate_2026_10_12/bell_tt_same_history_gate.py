#!/usr/bin/env python3
"""BQG Bell+shape-source history one-pass. Source Gauss singlets REAL;
clock Hamiltonian CHOSEN TOY, not full microscopic Lorentzian BQG.
"""
import numpy as np
import sympy as s
from math import sqrt,pi,exp
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Z=np.diag([1,-1]).astype(complex)
Y=np.array([[0,-1j],[1j,0]])
Id4=np.eye(4,dtype=complex)
def ket(bits):
 v=np.zeros(16,complex);v[int(bits,2)]=1;return v
a=(ket('0101')-ket('0110')-ket('1001')+ket('1010'))/2
b=(ket('0011')+ket('1100')-(ket('0101')+ket('0110')+ket('1001')+ket('1010'))/2)/sqrt(3)
K=np.column_stack((a,b))
assert np.linalg.norm(K.conj().T@K-I)<1e-12
for p in (X,Y,Z):
 J=np.zeros((16,16),complex)
 for k in range(4):
  arr=[I]*4;arr[k]=p/2
  z=arr[0]
  for y in arr[1:]:z=np.kron(z,y)
  J+=z
 assert np.linalg.norm(J@K)<1e-12
 for op in (X,Z):
  embedded=K@op@K.conj().T
  assert np.linalg.norm(J@embedded-embedded@J)<1e-12
assert np.linalg.norm(np.kron(K,K).conj().T@np.kron(K,K)-Id4)<1e-12
phi=np.array([1,0,0,1],complex)/sqrt(2)
psi=np.array([0,1,1,0],complex)/sqrt(2)
phim=np.array([1,0,0,-1],complex)/sqrt(2)
psim=np.array([0,1,-1,0],complex)/sqrt(2)
energies=[pi/4,pi/2,3*pi/4]
vectors=[psi,phim,psim]
H=sum(e*np.outer(v,v.conj()) for e,v in zip(energies,vectors))
local=[np.kron(Z,I),np.kron(X,I),
       np.kron(I,(Z+X)/sqrt(2)),np.kron(I,(Z-X)/sqrt(2))]
sources=np.column_stack((local[1]@phi,local[0]@phi))
Cs=[]
for tau in (0,1,2):
 T=np.outer(phi,phi.conj())+sum(exp(-tau*e)*np.outer(v,v.conj()) for e,v in zip(energies,vectors))
 Cs.append(sources.conj().T@T@sources)
C0,C1,C2=Cs
assert np.linalg.norm(C0-np.eye(2))<1e-12
assert np.linalg.norm(C2-C1@np.linalg.inv(C0)@C1)<1e-12
assert np.allclose(np.sort(-np.log(np.linalg.eigvalsh(C1))),energies[:2])
bell=[]
for N in (8,16,32):
 step=8/N
 R=np.outer(phi,phi.conj())+sum(np.exp(-1j*step*e)*np.outer(v,v.conj()) for e,v in zip(energies,vectors))
 assert np.linalg.norm(np.linalg.matrix_power(R,N)-Id4)<1e-12
 S=np.roll(np.eye(N,dtype=complex),1,axis=0);G=np.kron(S,R)
 V=np.vstack([np.linalg.matrix_power(R,t)/sqrt(N) for t in range(N)])
 assert np.linalg.norm(V.conj().T@V-Id4)<1e-12
 assert np.linalg.norm(G@V-V)<1e-12
 def dress(o):
  O=np.zeros((4*N,4*N),complex)
  for t in range(N):
   Rt=np.linalg.matrix_power(R,t)
   O[4*t:4*t+4,4*t:4*t+4]=Rt@o@Rt.conj().T
  assert np.linalg.norm(O@G-G@O)<1e-11
  assert np.linalg.norm(O@V-V@o)<1e-11
  return O
 A,Ap,B,Bp=map(dress,local)
 for a in (A,Ap):
  for b in (B,Bp):assert np.linalg.norm(a@b-b@a)<1e-11
 omega=V@phi
 result=float(np.vdot(omega,(A@B+A@Bp+Ap@B-Ap@Bp)@omega).real)
 assert abs(result-2*sqrt(2))<1e-11
 bell.append(result)
assert np.ptp(np.array(bell))<1e-12
# Full bare logical Pauli commutant is only scalars.
pa=[s.eye(2),s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
gen=[s.kronecker_product(pa[1],pa[0]),s.kronecker_product(pa[3],pa[0]),
     s.kronecker_product(pa[0],pa[1]),s.kronecker_product(pa[0],pa[3])]
base=[s.kronecker_product(p,q) for p in pa for q in pa]
comm=s.Matrix.hstack(*[s.Matrix.vstack(*[(g@b-b@g).reshape(16,1) for g in gen]) for b in base])
assert comm.rank()==15
# Compression of commuting bare A and B need not commute.
u=s.eye(4)[:,0];v=(s.eye(4)[:,1]+s.eye(4)[:,2])/s.sqrt(2)
W=s.Matrix.hstack(u,v);AA=s.kronecker_product(pa[3],pa[0]);BB=s.kronecker_product(pa[0],pa[1])
assert AA@BB==BB@AA
C=(W.T@AA@W)@(W.T@BB@W)-(W.T@BB@W)@(W.T@AA@W)
assert C==s.Matrix([[0,1/s.sqrt(2)],[-1/s.sqrt(2),0]])
assert abs(np.exp(-1j*pi/4)-np.exp(-1j*(pi/4+2*pi)))<1e-12
assert abs(np.exp(-pi/4)-np.exp(-pi/4-2*pi))>.45
print({'status':'PASS — CONTROL ONLY','bell_by_clock_N_8_16_32':bell,
       'C1':C1.real.tolist(),'C2':C2.real.tolist(),
       'moment_leakage_Fro':float(np.linalg.norm(C2-C1@C1)),
       'full_BQG_constraints':'UNAVAILABLE','physical_TT_pole':'OPEN','full_HDA':'OPEN'})
