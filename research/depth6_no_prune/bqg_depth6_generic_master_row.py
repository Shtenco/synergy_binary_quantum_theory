#!/usr/bin/env python3
from __future__ import annotations
import functools,math
import numpy as np
import bqg_depth6_generic_jucys_selector as J
import bqg_depth6_41 as Z
import bqg_block_apply_invariants as A
import bqg_depth4_symmetric_master_builder as X
import bqg_s5_intertwiner_action_gate as S

S2=tuple(p for p in Z.S3 if p[2]==2)
G2=len(S2)

def char(p,seed): return 1 if seed=='triv' else Z.perm_sign(p)
def flip(seed): return 'sign' if seed=='triv' else 'triv'
def basis3(q,seed): return Z.s3_triv(q) if seed=='triv' else Z.s3_sign(q)

@functools.lru_cache(None)
def s2_stabilizer(spins):
    spins=tuple(spins);return tuple(p for p in S2 if S.permute_spins(spins,p)==spins)
@functools.lru_cache(None)
def canon2(spins):
    spins=tuple(spins);best=None;bp=None
    for p in S2:
        t=S.permute_spins(spins,p)
        if best is None or t<best:best=t;bp=p
    return best,bp
@functools.lru_cache(None)
def s2_triv(spins):return Z._projector_basis(tuple(spins),s2_stabilizer(tuple(spins)),lambda p:1)
@functools.lru_cache(None)
def s2_sign(spins):return Z._projector_basis(tuple(spins),s2_stabilizer(tuple(spins)),Z.perm_sign)
def basis2(q,seed):return s2_triv(q) if seed=='triv' else s2_sign(q)

def reduced_s3_node(qin,Xcoord,v,seed,tol=0.0):
    if v not in (0,1):raise ValueError(v)
    qin=tuple(qin);Xcoord=np.asarray(Xcoord,complex)
    if Xcoord.ndim==1:Xcoord=Xcoord[:,None]
    Vin=basis3(qin,seed);B=Vin@Xcoord
    H=A.H_node_apply_ultrafast(qin,v,B,tol);tmp={};outseed=flip(seed)
    for so,M in H.items():
        q,p=Z.canon3(tuple(so));tq,TM=X.transform_coeff(tuple(so),tuple(p),M);assert tq==q
        TM=char(p,outseed)*TM;tmp[q]=tmp.get(q,0)+TM
    gin=len(Z.s3_stabilizer(qin));out={}
    for q,M in tmp.items():
        V=basis3(q,outseed)
        if V.shape[1]==0:continue
        gout=len(Z.s3_stabilizer(q));Y=math.sqrt(gout/gin)*(V.conj().T@M)
        if np.linalg.norm(Y)>tol:out[q]=Y
    return out

def raw_to_s2_coords(raw,seed,tol=0.0):
    qs=sorted({canon2(tuple(sp))[0] for sp in raw});out={}
    for q in qs:
        B=raw.get(q)
        if B is None:continue
        V=basis2(q,seed);n2=G2/len(s2_stabilizer(q));Y=math.sqrt(n2)*(V.conj().T@B)
        if np.linalg.norm(Y)>tol:out[q]=Y
    return out

def reduced_s2_node(qin,Xcoord,v,seed,tol=0.0):
    if v!=2:raise ValueError(v)
    qin=tuple(qin);Xcoord=np.asarray(Xcoord,complex)
    if Xcoord.ndim==1:Xcoord=Xcoord[:,None]
    Vin=basis2(qin,seed);B=Vin@Xcoord
    H=A.H_node_apply_ultrafast(qin,v,B,tol);tmp={};outseed=flip(seed)
    for so,M in H.items():
        q,p=canon2(tuple(so));tq,TM=X.transform_coeff(tuple(so),tuple(p),M);assert tq==q
        TM=char(p,outseed)*TM;tmp[q]=tmp.get(q,0)+TM
    gin=len(s2_stabilizer(qin));out={}
    for q,M in tmp.items():
        V=basis2(q,outseed)
        if V.shape[1]==0:continue
        gout=len(s2_stabilizer(q));Y=math.sqrt(gout/gin)*(V.conj().T@M)
        if np.linalg.norm(Y)>tol:out[q]=Y
    return out

def merge(acc,part,prefix,scale=1.0):
    for q,M in part.items():
        key=(prefix,tuple(q));Y=scale*M
        acc[key]=acc.get(key,0)+Y

def master_map(rep,key,W=None,tol=0.0):
    rep=tuple(rep);seed=J.CFG[key][0]
    if W is None:W=J.fast_selector(rep,key)
    m=W.shape[1]
    if m==0:return {},W
    bs,sl,d=J.layout(rep,seed);acc={}
    # v=0,1: S3 fixes both vertices; output type flips by Hamiltonian sign cocycle.
    for b in bs:
        seg=W[sl[b],:]
        if not seg.size or np.linalg.norm(seg)<tol:continue
        merge(acc,reduced_s3_node(b,seg,0,seed,tol),'v0')
        merge(acc,reduced_s3_node(b,seg,1,seed,tol),'v1')
    # v=2 represents the S3-related orbit {2,3,4}; Gram weight is exactly 3.
    raw=J.expand(rep,seed,W);c2=raw_to_s2_coords(raw,seed,tol)
    for b,seg in c2.items():
        merge(acc,reduced_s2_node(b,seg,2,seed,tol),'v2',math.sqrt(3.0))
    return {q:M for q,M in acc.items() if np.linalg.norm(M)>tol},W

def gram(mp):
    if not mp:return np.zeros((0,0),complex)
    m=next(iter(mp.values())).shape[1];G=np.zeros((m,m),complex)
    for M in mp.values():G+=M.conj().T@M
    return (G+G.conj().T)/2
