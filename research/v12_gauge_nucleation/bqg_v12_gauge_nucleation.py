#!/usr/bin/env python3
"""BQG V12 finite gauge-dressed link nucleation candidate.

Reproduces the SO(5) vector representation of
scripts/su2_quantum_link_vector5_gate.py with two charged endpoints.
It is a GAUGE-INVARIANT TOY MODEL, not the actual K5 physical projector.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.linalg import eigh

I2 = np.eye(2, dtype=complex)
SIG = [
    np.array([[0,1],[1,0]],complex),
    np.array([[0,-1j],[1j,0]],complex),
    np.array([[1,0],[0,-1]],complex),
]

def so5_link():
    def M(a,b):
        A=np.zeros((5,5),complex)
        A[a,b]=-1j; A[b,a]=1j
        return A
    eps=np.zeros((3,3,3),int)
    from itertools import permutations
    for p in permutations(range(3)):
        inv=sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
        eps[p]=-1 if inv%2 else 1
    J=[sum((0.5*eps[a,b,c]*M(b,c) for b in range(3) for c in range(3)),
           np.zeros((5,5),complex)) for a in range(3)]
    K=[M(a,3) for a in range(3)]
    L=[0.5*(J[a]+K[a]) for a in range(3)]
    R=[0.5*(J[a]-K[a]) for a in range(3)]
    return L,R,[M(a,4) for a in range(4)]

def build(g=0.4, gap=1.2):
    I3=np.eye(3,dtype=complex)
    I5=np.eye(5,dtype=complex)
    Jm=[]
    for sig in SIG:
        j=np.zeros((3,3),dtype=complex)
        j[1:,1:]=0.5*sig
        Jm.append(j)
    L,R,U=so5_link()
    GL=[np.kron(np.kron(Jm[k],I5),I3)+np.kron(np.kron(I3,L[k]),I3) for k in range(3)]
    GR=[np.kron(np.kron(I3,R[k]),I3)+np.kron(np.kron(I3,I5),Jm[k]) for k in range(3)]
    casimir=sum((G@G for G in GL+GR),np.zeros((45,45),complex))
    vals,vecs=eigh(casimir)
    singlets=vecs[:,vals<1e-9]
    vacuum=np.zeros(45,complex)
    vacuum[np.ravel_multi_index((0,4,0),(3,5,3))]=1
    assert singlets.shape[1]==2, vals[:5]
    basis=np.eye(45,dtype=complex)[:,np.ravel_multi_index((1,0,1),(3,5,3))]
    omega=singlets@(singlets.conj().T@basis)
    omega-=vacuum*np.vdot(vacuum,omega)
    omega/=np.linalg.norm(omega)
    projector=np.outer(omega,omega.conj())
    H=gap*projector-g*(np.outer(omega,vacuum.conj())+np.outer(vacuum,omega.conj()))
    phys=singlets@singlets.conj().T
    bare=np.kron(np.kron(I3,U[0]),I3)@vacuum
    checks=dict(
        gauss=float(max(np.linalg.norm(G@omega) for G in GL+GR)),
        commutator=float(max(np.linalg.norm(H@G-G@H) for G in GL+GR)),
        bare_gauss=float(np.linalg.norm(casimir@bare)),
        bare_physical_projection=float(np.linalg.norm(phys@bare)),
        singlet_dimension=int(singlets.shape[1]),
        casimir_gap=float(vals[2]),
        link_active_probability_on_omega=float(np.linalg.norm(omega.reshape(3,5,3)[:,:4,:])**2),
    )
    return H,vacuum,omega,checks

def run(g=0.4,gap=1.2):
    H,vacuum,omega,checks=build(g,gap)
    evals,V=eigh(H)
    times=[0,0.5,1,2,3,4,6,8]
    rows=[]
    Om=math.sqrt(gap*gap+4*g*g)
    maxerr=0.
    for t in times:
        ket=V@(np.exp(-1j*evals*t)*(V.conj().T@vacuum))
        p=float(abs(np.vdot(omega,ket))**2)
        analytic=(4*g*g/Om**2)*math.sin(Om*t/2)**2
        maxerr=max(maxerr,abs(p-analytic))
        rows.append(dict(time=t,link_creation_probability=p,analytic=analytic,
            vacuum_probability=float(abs(np.vdot(vacuum,ket))**2),
            norm=float(np.vdot(ket,ket).real)))
    checks['analytic_probability_error']=maxerr
    checks['max_probability']=4*g*g/Om**2
    checks['first_maximum_time']=math.pi/Om
    checks['Hamiltonian_hermiticity_error']=float(np.linalg.norm(H-H.conj().T))
    assert checks['gauss']<1e-12
    assert checks['commutator']<1e-12
    assert checks['bare_physical_projection']<1e-12
    assert checks['bare_gauss']>0.1
    assert maxerr<1e-12
    assert all(abs(row['norm']-1)<1e-12 for row in rows)
    return {'hypothesis':'finite gauge-dressed nucleation candidate, not derived BQG gravity',
        'g':g,'gap':gap,'Hilbert_dimension':45,'checks':checks,'time_series':rows,
        'limitations':['Endpoint charges and vacuum are selected by hand',
            'Coupling g is not derived from K5 constraints',
            'Potential vertex pair is selected before link creation',
            'Time is external Hamiltonian time, not BQG physical time',
            'No 3+1D limit, graviton or Einstein action is established']}

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    out=run()
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(out,indent=2,ensure_ascii=False))
