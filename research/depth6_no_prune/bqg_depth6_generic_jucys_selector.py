#!/usr/bin/env python3
from __future__ import annotations
import functools, math
import numpy as np
import bqg_depth6_41 as Z
import bqg_s5_intertwiner_action_gate as S
import bqg_depth4_symmetric_master_builder as X
import bqg_s5_irrep_decomposition_gate as D

S5=Z.S5; S3=Z.S3; G3=Z.G3

def trans(a,b):
    p=list(range(5)); p[a],p[b]=p[b],p[a]; return tuple(p)

J4_TERMS=tuple(trans(1,k) for k in (2,3,4))
J5_TERMS=tuple(trans(0,k) for k in (1,2,3,4))

# seed, J4 target/spectrum, J5 target/spectrum
CFG={
    '32':   ('triv',-1,(-1,3), 0,(-2,0,3)),
    '311':  ('triv',-1,(-1,3),-2,(-2,0,3)),
    '221':  ('sign', 1,(-3,1), 0,(-3,0,2)),
    '2111': ('sign',-3,(-3,1), 1,(-4,1)),
}

@functools.lru_cache(None)
def s3_bases_of_s5(rep):
    rep=tuple(rep)
    return tuple(sorted({Z.canon3(S.permute_spins(rep,p))[0] for p in S5}))

def seed_basis(q,seed):
    return Z.s3_triv(q) if seed=='triv' else Z.s3_sign(q)

def seed_char(p,seed):
    return 1 if seed=='triv' else Z.perm_sign(p)

@functools.lru_cache(None)
def layout(rep,seed):
    bs=s3_bases_of_s5(tuple(rep)); sl={}; off=0
    for q in bs:
        d=seed_basis(q,seed).shape[1]
        sl[q]=slice(off,off+d); off+=d
    return bs,sl,off

def expand(rep,seed,C,tol=0.0):
    bs,sl,d=layout(tuple(rep),seed); C=np.asarray(C,complex)
    if C.ndim==1: C=C[:,None]
    if C.shape[0]!=d: raise ValueError(('coord rows',C.shape,d))
    raw={}
    for q in bs:
        seg=C[sl[q],:]
        if not seg.size or not np.any(seg): continue  # exact-zero support removal only
        V=seed_basis(q,seed); B=V@seg
        n3=G3/len(Z.s3_stabilizer(q)); scale=1/math.sqrt(n3)
        seen={}
        for p in S3:
            ts,TB=X.transform_coeff(q,p,B)
            if ts not in seen:
                seen[ts]=scale*seed_char(p,seed)*TB
        if len(seen)!=int(n3):
            raise RuntimeError(('s3 orbit size',q,len(seen),n3))
        for ts,TB in seen.items():
            if ts in raw: raise RuntimeError(('s3 overlap',q,ts))
            raw[ts]=TB
    return raw

def reduce(rep,seed,raw,k):
    bs,sl,d=layout(tuple(rep),seed); Y=np.zeros((d,k),complex)
    for q in bs:
        B=raw.get(q)
        if B is None: continue
        V=seed_basis(q,seed); n3=G3/len(Z.s3_stabilizer(q))
        Y[sl[q],:]=math.sqrt(n3)*(V.conj().T@B)
    return Y

def apply_sum(rep,seed,C,terms):
    C=np.asarray(C,complex)
    if C.ndim==1: C=C[:,None]
    raw=expand(rep,seed,C); acc={}
    for p in terms:
        for sp,B in raw.items():
            ts,TB=X.transform_coeff(sp,p,B)
            acc[ts]=acc.get(ts,0)+TB
    return reduce(rep,seed,acc,C.shape[1])

def exact_mult(rep,key):
    rep=tuple(rep); H=Z.s5_stabilizer(rep); val=0j
    for p in H:
        ts,U=X._transform_matrix(rep,tuple(p))
        if ts!=rep: raise RuntimeError(('bad stabilizer',rep,p,ts))
        val += D.irrep_characters(p)[key]*np.trace(U)
    val/=len(H)
    if abs(val.imag)>1e-8 or abs(val.real-round(val.real))>1e-7:
        raise RuntimeError(('noninteger multiplicity',rep,key,val))
    return int(round(val.real))

def apply_poly(rep,seed,C,target,spectrum,terms):
    Y=np.asarray(C,complex)
    for a in spectrum:
        if a==target: continue
        Y=(apply_sum(rep,seed,Y,terms)-a*Y)/(target-a)
    return Y

def apply_selector_projector(rep,key,C):
    seed,j4,s4,j5,s5=CFG[key]
    Y=apply_poly(rep,seed,C,j4,s4,J4_TERMS)
    return apply_poly(rep,seed,Y,j5,s5,J5_TERMS)

def fast_selector(rep,key,extra=0,seednum=12345,tol=1e-10):
    rep=tuple(rep); seed=CFG[key][0]; _,_,d=layout(rep,seed); m=exact_mult(rep,key)
    if m==0: return np.zeros((d,0),complex)
    k=min(d,m+extra)
    salt={'32':1,'311':2,'221':3,'2111':4}[key]
    rng=np.random.default_rng(seednum + 17*sum((i+1)*x for i,x in enumerate(rep)) + 1000003*salt)
    C=rng.standard_normal((d,k))+1j*rng.standard_normal((d,k))
    Y=apply_selector_projector(rep,key,C)
    Q,R=np.linalg.qr(Y,mode='reduced')
    s=np.linalg.svd(R,compute_uv=False)
    rank_tol=max(Y.shape)*np.finfo(float).eps*max(float(s[0]) if len(s) else 1.,1.)*100
    rank=int((s>rank_tol).sum())
    if rank<m:
        raise RuntimeError(('selector rank',key,rep,d,m,k,rank,rank_tol))
    Q=Q[:,:m]
    seed,j4,s4,j5,s5=CFG[key]
    e4=np.linalg.norm(apply_sum(rep,seed,Q,J4_TERMS)-j4*Q)/max(np.linalg.norm(Q),1e-30)
    e5=np.linalg.norm(apply_sum(rep,seed,Q,J5_TERMS)-j5*Q)/max(np.linalg.norm(Q),1e-30)
    ep=np.linalg.norm(apply_selector_projector(rep,key,Q)-Q)/max(np.linalg.norm(Q),1e-30)
    oo=np.linalg.norm(Q.conj().T@Q-np.eye(m))
    if max(e4,e5,ep,oo)>tol:
        raise RuntimeError(('selector residual',key,rep,d,m,e4,e5,ep,oo))
    return Q
