#!/usr/bin/env python3
from __future__ import annotations
import functools, itertools, json
import numpy as np
import bqg_safe_boundary_master_gate as R
import bqg_s5_intertwiner_action_gate as S

IRREPS=('5','41','32','311','221','2111','11111')


def all_perms(): return tuple(itertools.permutations(range(5)))

def fixed_points(p): return sum(i==p[i] for i in range(5))

def perm_sign(p):
    inv=sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))
    return -1 if inv%2 else 1

def cycle_type(p):
    seen=set(); lens=[]
    for i in range(5):
        if i in seen: continue
        j=i; n=0
        while j not in seen:
            seen.add(j); n+=1; j=p[j]
        lens.append(n)
    return tuple(sorted(lens,reverse=True))

def irrep_characters(p):
    # Standard representation [4,1] is permutation rep minus the trivial line.
    chiV=fixed_points(p)-1
    p2=S.compose_perm(p,p)
    chiV2=fixed_points(p2)-1
    sym2=(chiV*chiV+chiV2)//2
    wedge2=(chiV*chiV-chiV2)//2
    chi32=sym2-1-chiV
    chi311=wedge2
    sg=perm_sign(p)
    return {
        '5':1,
        '41':chiV,
        '32':chi32,
        '311':chi311,
        '221':sg*chi32,
        '2111':sg*chiV,
        '11111':sg,
    }


def boundary_basis(): return R.basis_full_jhalf()

@functools.lru_cache(None)
def boundary_representation_matrix(p):
    p=tuple(p); basis=boundary_basis(); idx={k:i for i,k in enumerate(basis)}
    U=np.zeros((32,32),complex)
    for col,key in enumerate(basis):
        for kp,c in S.act_on_key(key,p).items(): U[idx[kp],col]+=c
    return U


def character_orthogonality():
    ps=all_perms(); C=np.array([[irrep_characters(p)[r] for r in IRREPS] for p in ps],float)
    G=(C.T@C)/len(ps)
    return {'max_offdiag':float(np.max(np.abs(G-np.eye(len(IRREPS))) * (1-np.eye(len(IRREPS))))),
            'max_diag_error':float(np.max(np.abs(np.diag(G)-1))),
            'gram':G.tolist()}


def boundary_decomposition():
    ps=all_perms(); chars=[]
    for p in ps: chars.append(float(np.trace(boundary_representation_matrix(p)).real))
    mult={}
    for r in IRREPS:
        m=sum(chars[i]*irrep_characters(p)[r] for i,p in enumerate(ps))/len(ps)
        mult[r]=int(round(m)) if abs(m-round(m))<1e-9 else float(m)
    dims={r:irrep_characters(S.identity_perm())[r] for r in IRREPS}
    recon=sum(dims[r]*mult[r] for r in IRREPS)
    by_class={}
    for p,ch in zip(ps,chars):
        ct=str(cycle_type(p)); by_class.setdefault(ct,[]).append(ch)
    class_chars={k:float(np.mean(v)) for k,v in by_class.items()}
    class_spread={k:float(np.max(v)-np.min(v)) for k,v in by_class.items()}
    return {'representation_dimension':32,'multiplicities':mult,'irrep_dimensions':dims,
            'dimension_reconstructed':int(round(recon)), 'character_by_cycle_type':class_chars,
            'max_character_spread_within_class':max(class_spread.values()),
            'scope':'exact 32D all-j=1/2 K5 boundary representation including SU(2) edge-orientation cocycle and local recoupling'}


def run():
    return {'status':'S5 irrep decomposition of the 32D K5 boundary sector',
            'orthogonality':character_orthogonality(),'decomposition':boundary_decomposition()}

if __name__=='__main__': print(json.dumps(run(),indent=2))

