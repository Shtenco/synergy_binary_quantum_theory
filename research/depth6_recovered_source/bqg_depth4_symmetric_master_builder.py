#!/usr/bin/env python3
from __future__ import annotations
import functools, itertools, math, numpy as np
import bqg_safe_boundary_master_gate as R
import bqg_s5_intertwiner_action_gate as S
import bqg_s5_irrep_decomposition_gate as D


def perm_sign(p):
    inv=sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))
    return -1 if inv%2 else 1

@functools.lru_cache(None)
def coset_reps(spins):
    spins=tuple(spins); reps={}
    for p in D.all_perms():
        t=S.permute_spins(spins,p)
        reps.setdefault(t,p)
    return tuple(reps.values())

@functools.lru_cache(None)
def _basis_labels(spins):
    spins=tuple(spins)
    opts=[R.allowed_k2_t(*R.local_spins(spins,v)) for v in R.VERT]
    return tuple(itertools.product(*opts))

@functools.lru_cache(None)
def _transform_matrix(spins,p):
    spins=tuple(spins);p=tuple(p);tsp=S.permute_spins(spins,p)
    bi=_basis_labels(spins); bo=_basis_labels(tsp); oi={k:i for i,k in enumerate(bo)}
    U=np.zeros((len(bo),len(bi)),complex)
    for j,Ks in enumerate(bi):
        for (sp2,K2),c in S.act_on_key((spins,Ks),p).items():
            assert sp2==tsp
            U[oi[K2],j]+=c
    return tsp,U

def transform_coeff(spins,p,V):
    tsp,U=_transform_matrix(tuple(spins),tuple(p))
    return tsp,U@np.asarray(V,complex)

# Fast exact factorized coefficient transport: tensor product of local recoupling maps.
def transform_coeff_fast(spins,p,V):
    spins=tuple(spins);p=tuple(p);V=np.asarray(V,complex)
    tsp=S.permute_spins(spins,p)
    old_dims=[len(R.allowed_k2_t(*R.local_spins(spins,v))) for v in R.VERT]
    new_dims=[len(R.allowed_k2_t(*R.local_spins(tsp,w))) for w in R.VERT]
    X=V.reshape(tuple(old_dims)+(V.shape[1],))
    # Apply old-vertex v map, keeping its axis position v for now.
    for v in R.VERT:
        M,_,_=S.local_recoupling_matrix(spins,v,p)
        Y=np.tensordot(M,X,axes=(1,v))
        # Y axes: [new_v_axis] + all X axes except v. Move front back to v.
        X=np.moveaxis(Y,0,v)
    pinv=S.inverse_perm(p)
    X=np.transpose(X,tuple(pinv[w] for w in R.VERT)+(5,))
    X*=S.edge_orientation_phase(spins,p)
    return tsp,X.reshape((int(np.prod(new_dims)),V.shape[1]))

transform_coeff=transform_coeff_fast
