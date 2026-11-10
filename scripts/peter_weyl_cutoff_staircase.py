#!/usr/bin/env python3
"""Exact finite SU(2) Peter-Weyl link/cutoff staircase regression test.

Constructs H_link^{Jmax}=direct_sum_{j<=Jmax} V_j^L tensor V_j^R and the
fundamental multiplication operator from Clebsch-Gordan coefficients.  It
checks exact left/right gauge covariance and evolves the gauge vacuum under a
sum of six two-link Wilson loops on a four-link dipole.  The vacuum history
W^n must stabilize exactly once Jmax >= n/2.
"""
from __future__ import annotations
import argparse,itertools,json,math
from functools import lru_cache
from pathlib import Path
import numpy as np
import sympy as sp
from sympy.physics.wigner import clebsch_gordan

@lru_cache(None)
def _cg(a,b,c,d,e,f):
    return complex(clebsch_gordan(*[sp.Rational(x,2) for x in (a,b,c,d,e,f)]).evalf())
def CG(j1,j2,J,m1,m2,M):return _cg(*(int(round(2*x)) for x in (j1,j2,J,m1,m2,M)))
def mvals(j):return np.array([j-r for r in range(int(round(2*j))+1)],float)
def spin_mats(j):
    ms=mvals(j);d=len(ms);Jz=np.diag(ms);Jp=np.zeros((d,d),complex)
    if j>0:
        for c,m in enumerate(ms):
            if c>0:Jp[c-1,c]=math.sqrt(j*(j+1)-m*(m+1))
    Jm=Jp.conj().T
    return ms,[(Jp+Jm)/2,(Jp-Jm)/(2j) if j>0 else np.zeros((1,1),complex),Jz]

def build_link(jmax):
    js=np.arange(0,int(round(2*jmax))+1)/2;basis=[]
    for j in js:
        for m in mvals(j):
            for n in mvals(j):basis.append((float(j),float(m),float(n)))
    idx={b:i for i,b in enumerate(basis)};D=len(basis)
    L=[np.zeros((D,D),complex) for _ in range(3)];R=[np.zeros((D,D),complex) for _ in range(3)]
    for j in js:
        ms,J=spin_mats(float(j))
        for mi,m in enumerate(ms):
            for ni,n in enumerate(ms):
                col=idx[(float(j),float(m),float(n))]
                for pi,p in enumerate(ms):
                    for a in range(3):
                        L[a][idx[(float(j),float(p),float(n))],col]=J[a][pi,mi]
                        R[a][idx[(float(j),float(m),float(p))],col]=J[a][pi,ni]
    fs=mvals(.5);U={};Ud={}
    for ai,a in enumerate(fs):
        for bi,b in enumerate(fs):
            M=np.zeros((D,D),complex)
            for j,m,n in basis:
                col=idx[(j,m,n)]
                for J in (j-.5,j+.5):
                    if J<0 or J>jmax+1e-12:continue
                    MM=m+a;NN=n+b
                    if abs(MM)>J+1e-12 or abs(NN)>J+1e-12:continue
                    c=math.sqrt((2*j+1)/(2*J+1))*CG(j,.5,J,m,a,MM)*CG(j,.5,J,n,b,NN)
                    if c != 0.0:M[idx[(float(J),float(MM),float(NN))],col]+=c
            U[(ai,bi)]=M;Ud[(ai,bi)]=M.conj().T
    return basis,idx,L,R,U,Ud

def sparse_transitions(M):
    return {c:[(int(r),M[r,c]) for r in np.flatnonzero(abs(M[:,c])>1e-14)] for c in range(M.shape[1])}
def evolve(jmax,nmax=5):
    basis,idx,L,R,U,Ud=build_link(jmax);D=len(basis);Ut={k:sparse_transitions(v) for k,v in U.items()};Dt={k:sparse_transitions(v) for k,v in Ud.items()}
    zero=idx[(0.,0.,0.)];state={(zero,zero,zero,zero):1+0j};pairs=list(itertools.combinations(range(4),2));rows=[]
    for n in range(1,nmax+1):
        out={}
        for ds,amp in state.items():
            for e1,e2 in pairs:
                for ab in itertools.product(range(2),repeat=2):
                    for o1,c1 in Ut[ab][ds[e1]]:
                        for o2,c2 in Dt[ab][ds[e2]]:
                            nd=list(ds);nd[e1]=o1;nd[e2]=o2;nd=tuple(nd);out[nd]=out.get(nd,0)+amp*c1*c2
        state={k:v for k,v in out.items() if abs(v)>1e-13};n2=float(sum(abs(v)**2 for v in state.values()));maxj=max(max(basis[d][0] for d in ds) for ds in state)
        rows.append({"power":n,"norm_squared":n2,"support":len(state),"max_spin_reached":maxj})
    return D,rows

def run():
    J=[.5,1,1.5,2,2.5];allrows={};cov={}
    for jm in J:
        basis,idx,L,R,U,Ud=build_link(jm);_,Jf=spin_mats(.5);e=0.0
        for a in range(3):
            for i,j in itertools.product(range(2),repeat=2):
                e=max(e,float(np.linalg.norm(L[a]@U[(i,j)]-U[(i,j)]@L[a]-sum(Jf[a][p,i]*U[(p,j)] for p in range(2)))))
                e=max(e,float(np.linalg.norm(R[a]@U[(i,j)]-U[(i,j)]@R[a]-sum(U[(i,p)]*Jf[a][p,j] for p in range(2)))))
        D,rows=evolve(jm);cov[str(jm)]=e;allrows[str(jm)]={"link_dimension":D,"rows":rows}
    stable=[6,120,4650,254604,16807392]
    checks=[]
    for n,target in enumerate(stable,1):
        jm=n/2;obs=allrows[str(jm)]["rows"][n-1]["norm_squared"];checks.append({"power":n,"minimum_exact_Jmax":jm,"observed":obs,"target_integer":target,"error":abs(obs-target)})
    passed=max(cov.values())<1e-12 and max(x["error"] for x in checks)<1e-7
    return {"status":"exact Peter-Weyl cutoff staircase","passed":bool(passed),"gauge_covariance_max_errors":cov,"cutoff_runs":allrows,"stable_vacuum_W_norm_squared":checks,"theorem":"A word of r fundamental holonomies is cutoff-exact on an input spin-j sector whenever j+r/2 <= Jmax. For the vacuum Wilson history W^n this gives the observed threshold Jmax=n/2."}

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--output",type=Path);a=ap.parse_args();o=run();t=json.dumps(o,indent=2);print(t)
    if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(t+"\n",encoding="utf-8")
    return 0 if o["passed"] else 1
if __name__=="__main__":raise SystemExit(main())
