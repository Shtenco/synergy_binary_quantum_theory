#!/usr/bin/env python3
"""Exact SO5 transporter Floquet control, NOT physical BQG Lorentzian evolution."""
import json
import numpy as np
from scipy.linalg import expm
def U(a):
    z=np.zeros((5,5),complex);z[a,4]=-1j;z[4,a]=1j;return z
u=[U(a) for a in range(4)];cyc=(0,1,3,2)
edges=[(cyc[k],cyc[(k+1)%4]) for k in range(4)]
P=np.diag([1,1,1,1,0]);Q=np.eye(5)-P
def comm(a,b):return a@b-b@a
def word(eps,ed=edges):
    out=np.eye(5,dtype=complex)
    for a,b in ed:
        out=out@expm(-1j*eps*u[b])@expm(-1j*eps*u[a])@expm(1j*eps*u[b])@expm(1j*eps*u[a])
    return out
A=sum((comm(u[b],u[a]) for a,b in edges),np.zeros((5,5),complex))
K=A/(2j)
cubic=sum((comm(u[b]+u[a],comm(u[b],u[a])) for a,b in edges),np.zeros((5,5),complex))
assert np.linalg.norm(cubic)<1e-12
for a,b in edges:
    assert np.linalg.norm(comm(u[b]+u[a],comm(u[b],u[a]))-u[a]+u[b])<1e-12
assert np.allclose(np.linalg.eigvalsh(K),[-1,0,0,0,1],atol=1e-12)
e=[.4,.2,.1,.05,.025]
def log_slope(xs,ys):return float(np.polyfit(np.log(xs),np.log(ys),1)[0])
closed=[word(t) for t in e]
open_=[word(t,edges[:3]) for t in e]
AO=sum((comm(u[b],u[a]) for a,b in edges[:3]),np.zeros((5,5),complex))
errs=[np.linalg.norm(W-expm(-t*t*A),'fro') for t,W in zip(e,closed)]
errs_open=[np.linalg.norm(W-expm(-t*t*AO),'fro') for t,W in zip(e,open_)]
leaks=[np.linalg.norm(Q@W@P,'fro') for W in closed]
leaks_open=[np.linalg.norm(Q@W@P,'fro') for W in open_]
slopes=[log_slope(e[1:],data[1:]) for data in (errs,errs_open,leaks,leaks_open)]
assert abs(slopes[0]-4)<.1 and abs(slopes[1]-3)<.1
assert abs(slopes[2]-5)<.1 and abs(slopes[3]-3)<.1
S=np.zeros((5,5),complex)
for a,b in edges:S[b,a]=1
Ks=(S+S.conj().T)/2
G=(np.pi/2)*K+np.pi*(Ks@Ks-Ks)/2
assert np.linalg.norm(expm(1j*G)-S-Q)<1e-12
assert np.linalg.norm(K.conj()+K)<1e-12
n=[8,16,32,64,128,256];ref=expm(-A)
r=[np.linalg.matrix_power(word(1/np.sqrt(m)),m) for m in n]
global_err=[np.linalg.norm(v-ref,'fro') for v in r]
global_leak=[np.linalg.norm(Q@v@P,'fro') for v in r]
cont_slopes=[log_slope(n,data) for data in (global_err,global_leak)]
assert abs(cont_slopes[0]+1)<.1 and abs(cont_slopes[1]+1.5)<.1
print(json.dumps({'status':'MATHEMATICAL_CONTROL_PASS_NOT_PHYSICAL_BQG',
  'single_word_slopes_epsilon4_epsilon3_epsilon5_epsilon3':slopes,
  'n_cycle_slopes_global_error_and_Q_leak':cont_slopes,
  'continuum_n':n,'operator_errors':global_err,'Q_leakage':global_leak,
  'K_chi_spectrum':np.linalg.eigvalsh(K).tolist(),
  'physical_HDA':'OPEN','TT_pole':'OPEN','microscopic_carry_selection':'OPEN'},indent=2))
