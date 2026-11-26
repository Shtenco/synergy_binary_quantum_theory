#!/usr/bin/env python3
from __future__ import annotations
import functools, itertools
import numpy as np
import bqg_safe_boundary_master_gate as R


def allowed_pair(s1,s2):
    return tuple(range(abs(s1-s2),s1+s2+1,2))

@functools.lru_cache(None)
def multiplicity_labels(spins,J):
    s1,s2,s3,s4=spins
    out=[]
    for K in allowed_pair(s1,s2):
        for L in allowed_pair(s3,s4):
            if abs(K-L)<=J<=K+L and (J-abs(K-L))%2==0:
                out.append((K,L))
    return tuple(out)

@functools.lru_cache(None)
def allowed_total_J(spins):
    vals=set()
    s1,s2,s3,s4=spins
    for K in allowed_pair(s1,s2):
        for L in allowed_pair(s3,s4):
            vals.update(range(abs(K-L),K+L+1,2))
    return tuple(sorted(vals))

@functools.lru_cache(None)
def coupled_tensor(spins,K,L,J,M):
    spins=tuple(spins); dims=tuple(s+1 for s in spins)
    T=np.zeros(dims,complex)
    mvs=[R.m2vals_t(s) for s in spins]
    for idx in itertools.product(*[range(d) for d in dims]):
        m1,m2,m3,m4=[mvs[a][idx[a]] for a in range(4)]
        MK=m1+m2; ML=m3+m4
        if MK+ML!=M: continue
        if abs(MK)>K or abs(ML)>L: continue
        c1=R.cg2(spins[0],spins[1],K,m1,m2,MK)
        if abs(c1)<1e-15: continue
        c2=R.cg2(spins[2],spins[3],L,m3,m4,ML)
        if abs(c2)<1e-15: continue
        c3=R.cg2(K,L,J,MK,ML,M)
        if abs(c3)<1e-15: continue
        T[idx]=c1*c2*c3
    return T


def apply_q_triple_tensor(T,spins,axes):
    mats=[R.spin_mats_cached(s) for s in spins]
    out=np.zeros_like(T,dtype=complex)
    a,b,c=axes
    for i,j,k in itertools.product(range(3),repeat=3):
        e=R.EPS3[i,j,k]
        if not e: continue
        X=R.apply_axis_np(T,a,mats[a][i])
        X=R.apply_axis_np(X,b,mats[b][j])
        X=R.apply_axis_np(X,c,mats[c][k])
        out += e*X
    return out


def apply_q_full_tensor(T,spins):
    out=np.zeros_like(T,dtype=complex)
    for r in range(4):
        axes=tuple(i for i in range(4) if i!=r)
        out += ((-1)**r)*apply_q_triple_tensor(T,spins,axes)
    return out

@functools.lru_cache(None)
def reduced_q_block(spins,J):
    spins=tuple(spins); labels=multiplicity_labels(spins,J); n=len(labels)
    Q=np.zeros((n,n),complex); M=J
    basis=[coupled_tensor(spins,K,L,J,M) for K,L in labels]
    qb=[apply_q_full_tensor(B,spins) for B in basis]
    for a,A in enumerate(basis):
        for b,QB in enumerate(qb): Q[a,b]=np.vdot(A,QB)
    return 0.5*(Q+Q.conj().T)

@functools.lru_cache(None)
def reduced_volume_block(spins,J):
    Q=reduced_q_block(tuple(spins),J)
    ev,U=np.linalg.eigh(Q)
    return (U*np.sqrt(np.abs(ev)))@U.conj().T


def basis_diagnostics(spins):
    spins=tuple(spins); vecs=[]
    for J in allowed_total_J(spins):
        for M in R.m2vals_t(J):
            for K,L in multiplicity_labels(spins,J):
                vecs.append(coupled_tensor(spins,K,L,J,M).reshape(-1))
    B=np.column_stack(vecs) if vecs else np.zeros((int(np.prod([s+1 for s in spins])),0),complex)
    G=B.conj().T@B
    return {'tensor_dimension':int(np.prod([s+1 for s in spins])),'basis_dimension':B.shape[1],
            'orthonormality_error':float(np.linalg.norm(G-np.eye(B.shape[1])))}


def apply_volume_tensor(T,spins):
    spins=tuple(spins); out=np.zeros_like(T,dtype=complex)
    for J in allowed_total_J(spins):
        labels=multiplicity_labels(spins,J); V=reduced_volume_block(spins,J)
        for M in R.m2vals_t(J):
            Bs=[coupled_tensor(spins,K,L,J,M) for K,L in labels]
            c=np.array([np.vdot(B,T) for B in Bs],complex)
            d=V@c
            for a,B in enumerate(Bs): out += d[a]*B
    return out


def apply_volume_batch(T,spins):
    """Vectorized exact volume action on a batch T[B,d1,d2,d3,d4]."""
    spins=tuple(spins); Bn=T.shape[0]; dims=tuple(s+1 for s in spins)
    X=T.reshape(Bn,-1); Y=np.zeros_like(X,dtype=complex)
    for J in allowed_total_J(spins):
        labels=multiplicity_labels(spins,J); V=reduced_volume_block(spins,J)
        for M in R.m2vals_t(J):
            Bs=np.stack([coupled_tensor(spins,K,L,J,M).reshape(-1) for K,L in labels],axis=0)
            c=X@Bs.conj().T
            d=c@V.T
            Y += d@Bs
    return Y.reshape((Bn,)+dims)
