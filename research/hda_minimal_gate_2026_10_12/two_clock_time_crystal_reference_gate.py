#!/usr/bin/env python3
"""Two-clock relative-rank, finite 2T no-go, phase carry and Floquet negatives.
The full matrix-free N=8,16,32 Gauss/Bell producer is in accompanying
BQG_2T_TIME_CRYSTAL_ONEPASS_20261012.zip research package.
These are reference controls only, NOT a BQG physical HDA/TT result.
"""
import math
import numpy as np
import sympy as s
I=np.eye(2);J=np.array([[0,-1],[1,0]],dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]])
Z=np.diag([1,-1])
R1=np.kron(J,I);R2=np.kron(I,J)
phi=np.array([1,0,0,1],complex)/math.sqrt(2)
assert np.allclose(R1@R2@phi,phi)
assert np.allclose(R1@phi,-R2@phi)
K=[math.pi/2*np.kron(Y,I),math.pi/2*np.kron(I,Y)]
QFI=np.array([[4*np.real(np.vdot(phi,K[a]@K[b]@phi)-np.vdot(phi,K[a]@phi)*np.vdot(phi,K[b]@phi)) for b in range(2)] for a in range(2)])
assert np.allclose(QFI,math.pi**2*np.array([[1,-1],[-1,1]]))
assert np.linalg.matrix_rank(QFI)==1
O=[np.kron(Z,I),np.kron(X,I)]
C=[]
for t in (0,1,2):
 U=np.linalg.matrix_power(R1,t); Ui=np.linalg.matrix_power(R1,-t)
 C.append(np.array([[np.vdot(phi,Ui@A@U@B@phi) for B in O] for A in O]))
assert all(np.allclose(C[t],(-1)**t*np.eye(2)) for t in range(3))
assert np.allclose(C[2]-C[1]@C[1],0)
assert np.linalg.eigvalsh(C[1]).max() < 0
def order(T):
 a=list(range(len(T)))
 for i in range(1,17):
  a=[T[v] for v in a]
  if a==list(range(len(T))):return i
ind=[2*((k+1)%4)+(1-bit) for k in range(4) for bit in range(2)]
carry=[2*k+1 if bit==0 else 2*((k+1)%4) for k in range(4) for bit in range(2)]
assert order(ind)==4 and order(carry)==8
for n in (1,8,16,32):
 eps=.08
 assert abs(math.cos(n*(math.pi+eps))-((-1)**n)*math.cos(n*eps))<1e-12
x=s.symbols("x0:6");p=s.symbols("p0:6")
g=[1,1,1,1,-1,-1]
Q11=sum(g[i]*x[i]**2 for i in range(6))
Q12=sum(x[i]*p[i] for i in range(6))
Q22=sum(g[i]*p[i]**2 for i in range(6))
def PB(a,b):return s.expand(sum(s.diff(a,x[i])*s.diff(b,p[i])-s.diff(a,p[i])*s.diff(b,x[i]) for i in range(6)))
assert PB(Q11,Q12)==2*Q11
assert PB(Q12,Q22)==2*Q22
assert PB(Q11,Q22)==4*Q12
print("PASS: 2-clock relative QFI rank1, unitary period2 not positive TT, independent clock != phase carry, external Bars Poisson algebra")
print("NOT PROVED: BQG microscopic Sp(2,R), full Lorentzian HDA, physical massless TT, robust interacting time crystal")
