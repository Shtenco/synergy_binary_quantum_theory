#!/usr/bin/env python3
"""Exact Schur/TT two-source finite MASTER gate. Mathematical toy only, not physical BQG."""
import sympy as s
import numpy as np
D=s.diag(2,3)
K=s.Matrix([[s.Rational(1,3),0,s.Rational(1,2)],[0,s.Rational(2,5),s.Rational(1,3)]])
F=s.diag(0,0,5)
A=F+K.T*D*K;B=K.T*D
M=A.row_join(B).col_join(B.T.row_join(D))
assert M.rank()==3 and A.rank()==3
assert A-B*D.inv()*B.T==F
rootD=s.diag(s.sqrt(2),s.sqrt(3))
C=s.diag(0,0,s.sqrt(5)).row_join(s.zeros(3,2)).col_join((rootD*K).row_join(rootD))
assert C.T*C==M
CP=C[:,:3];CQ=C[:,3:]
assert CP.T*(s.eye(5)-CQ*(CQ.T*CQ).inv()*CQ.T)*CP==F
C1=CQ[:,:1]
F1=CP.T*(s.eye(5)-C1*(C1.T*C1).inv()*C1.T)*CP
assert F1.rank()==2 and F.rank()==1
E=s.Matrix([[1,0],[0,1],[0,0]])
W=E.col_join(-K*E)
P=W*(W.T*W).inv()*W.T
assert P*P==P and M*P==s.zeros(5)
U=E.col_join(s.zeros(2,2))
G=s.simplify(U.T*P*U)
assert G==s.diag(s.Rational(9,10),s.Rational(25,29))
# True nonzero compression does not exclude dressed zero vectors.
assert A.det()>0
# Dropping output rows creates a SPURIOUS third zero.
cut=C[3:,:];cP=cut[:,:3];cQ=cut[:,3:]
Fcut=cP.T*(s.eye(2)-cQ*(cQ.T*cQ).inv()*cQ.T)*cP
assert Fcut==s.zeros(3)
# Shortcut uses only two complete H_v U columns.
K0=U.T*U;K1=U.T*M*U
low=s.simplify(K0-K1)  # delta=1
assert low==s.diag(s.Rational(7,9),s.Rational(13,25))
assert s.simplify(G-low)==s.diag(s.Rational(11,90),s.Rational(248,725))
eigs=np.linalg.eigvalsh(np.array(M).astype(float))
assert eigs[2]>1
print('PASS exact rational Schur, 2 physical toy source directions, output-truncation no-go')
print('coarse minimal eig',np.linalg.eigvalsh(np.array(A).astype(float)).min())
print('first positive full master eig',eigs[2])
print('full physical master TT in BQG: NOT CALCULATED; Lorentzian images and independent D missing')
