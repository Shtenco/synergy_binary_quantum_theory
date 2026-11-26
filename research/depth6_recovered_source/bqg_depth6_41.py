#!/usr/bin/env python3
from __future__ import annotations
import functools, itertools, math, json, time
from pathlib import Path
import numpy as np

import bqg_block_apply_invariants as A
import bqg_depth4_symmetric_master_builder as X
import bqg_s5_intertwiner_action_gate as S
import bqg_s5_irrep_decomposition_gate as D

S5=tuple(D.all_perms())
S4=tuple(p for p in S5 if p[0]==0)
G4=len(S4)


def perm_sign(p):
    inv=sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))
    return -1 if inv&1 else 1

@functools.lru_cache(None)
def s4_stabilizer(spins):
    spins=tuple(spins)
    return tuple(p for p in S4 if S.permute_spins(spins,p)==spins)

@functools.lru_cache(None)
def s5_stabilizer(spins):
    spins=tuple(spins)
    return tuple(p for p in S5 if S.permute_spins(spins,p)==spins)

@functools.lru_cache(None)
def canon4(spins):
    spins=tuple(spins)
    best=None;bp=None
    for p in S4:
        t=S.permute_spins(spins,p)
        if best is None or t<best:
            best=t;bp=p
    return best,bp

@functools.lru_cache(None)
def canon5(spins):
    spins=tuple(spins)
    best=None;bp=None
    for p in S5:
        t=S.permute_spins(spins,p)
        if best is None or t<best:
            best=t;bp=p
    return best,bp


def _projector_basis(spins,perms,char,tol=5e-9):
    spins=tuple(spins); n=len(A.gauss_basis_labels(spins))
    P=np.zeros((n,n),complex)
    for p in perms:
        ts,U=X._transform_matrix(spins,tuple(p))
        if ts!=spins: raise RuntimeError('bad stabilizer')
        P += char(p)*U
    P/=len(perms); P=(P+P.conj().T)/2
    w,V=np.linalg.eigh(P); B=V[:,w>0.5]
    if B.size:
        rr=np.linalg.norm(P@B-B)/max(np.linalg.norm(B),1e-30)
        oo=np.linalg.norm(B.conj().T@B-np.eye(B.shape[1]))
        if rr>tol or oo>tol: raise RuntimeError(('projector',spins,rr,oo))
    return B

@functools.lru_cache(None)
def s4_triv(spins): return _projector_basis(tuple(spins),s4_stabilizer(tuple(spins)),lambda p:1)
@functools.lru_cache(None)
def s4_sign(spins): return _projector_basis(tuple(spins),s4_stabilizer(tuple(spins)),perm_sign)
@functools.lru_cache(None)
def s5_triv(spins): return _projector_basis(tuple(spins),s5_stabilizer(tuple(spins)),lambda p:1)

@functools.lru_cache(None)
def s4_bases_of_s5(rep):
    rep=tuple(rep)
    return tuple(sorted({canon4(S.permute_spins(rep,p))[0] for p in S5}))

@functools.lru_cache(None)
def perm_rep_to(rep,target):
    rep=tuple(rep);target=tuple(target)
    for p in S5:
        if S.permute_spins(rep,p)==target:return p
    raise RuntimeError('no perm')


@functools.lru_cache(None)
def s5_trivial_embedding_in_s4(rep):
    """Return S4 bases, row slices, E5. Columns of E5 are normalized S5-trivial states
    expressed in normalized S4-trivial orbit coordinates."""
    rep=tuple(rep); bases=s4_bases_of_s5(rep); V5=s5_triv(rep)
    n5=120/len(s5_stabilizer(rep))
    parts=[];slices={};off=0
    for b in bases:
        V4=s4_triv(b); d=V4.shape[1]; slices[b]=slice(off,off+d);off+=d
        p=perm_rep_to(rep,b); tb,TV=X.transform_coeff(rep,p,V5); assert tb==b
        n4=24/len(s4_stabilizer(b))
        Y=math.sqrt(n4/n5)*(V4.conj().T@TV)
        parts.append(Y)
    E=np.vstack(parts) if parts else np.zeros((0,V5.shape[1]),complex)
    if E.size:
        err=np.linalg.norm(E.conj().T@E-np.eye(E.shape[1]))
        if err>2e-8: raise RuntimeError(('E5 orth',rep,E.shape,err))
    return bases,slices,E


def m41_for_rep(rep):
    bases,slices,E=s5_trivial_embedding_in_s4(tuple(rep))
    d4=sum(s.stop-s.start for s in slices.values()); return d4-E.shape[1]


@functools.lru_cache(None)
def complement41(rep):
    bases,slices,E=s5_trivial_embedding_in_s4(tuple(rep)); d=E.shape[0];m5=E.shape[1]
    if m5==0:return bases,slices,np.eye(d,dtype=complex)
    Q,R=np.linalg.qr(E,mode='complete')
    # fix arbitrary phases not needed; fail closed on complement and orthogonality
    C=Q[:,m5:]
    e1=np.linalg.norm(E.conj().T@C) if C.size else 0.
    e2=np.linalg.norm(C.conj().T@C-np.eye(C.shape[1])) if C.size else 0.
    if e1>2e-8 or e2>2e-8:raise RuntimeError(('C41',rep,e1,e2))
    return bases,slices,C


def reduced_s4_triv_to_sign(qin,Xcoord,tol=1e-11,selected=None):
    """H0 on normalized S4-trivial orbit coordinates Xcoord, returned in normalized
    S4-sign output orbit coordinates keyed by canonical output spin q."""
    qin=tuple(qin); Xcoord=np.asarray(Xcoord,complex)
    Vin=s4_triv(qin)
    if Xcoord.ndim==1:Xcoord=Xcoord[:,None]
    if Xcoord.shape[0]!=Vin.shape[1]:raise ValueError(('rows',qin,Xcoord.shape,Vin.shape))
    B=Vin@Xcoord
    H=A.H_node_apply_ultrafast(qin,0,B,tol)
    tmp={}; sel=set(map(tuple,selected)) if selected is not None else None
    for so,W in H.items():
        q,p=canon4(tuple(so))
        if sel is not None and q not in sel:continue
        tq,TW=X.transform_coeff(tuple(so),tuple(p),W); assert tq==q
        TW=perm_sign(p)*TW
        tmp[q]=tmp.get(q,0)+TW
    gin=len(s4_stabilizer(qin));out={}
    for q,Z in tmp.items():
        Vout=s4_sign(q)
        if Vout.shape[1]==0:continue
        gout=len(s4_stabilizer(q))
        Y=math.sqrt(gout/gin)*(Vout.conj().T@Z)
        if np.linalg.norm(Y)>tol:out[q]=Y
    return out


def explicit_s4_triv_to_sign(qin,Xcoord,tol=1e-11):
    qin=tuple(qin);Xcoord=np.asarray(Xcoord,complex)
    if Xcoord.ndim==1:Xcoord=Xcoord[:,None]
    Vin=s4_triv(qin);B=Vin@Xcoord;n4=24/len(s4_stabilizer(qin));scale=1/math.sqrt(n4)
    members={}
    for p in S4:
        ts,TB=X.transform_coeff(qin,p,B)
        if ts not in members:members[ts]=scale*TB
    full={}
    for sp,BB in members.items():
        H=A.H_node_apply_ultrafast(sp,0,BB,tol)
        for so,W in H.items():full[so]=full.get(so,0)+W
    out={}
    for so in list(full):
        q,_=canon4(tuple(so))
        if q!=tuple(so):continue
        Vout=s4_sign(q)
        if Vout.shape[1]==0:continue
        nout=24/len(s4_stabilizer(q))
        Y=math.sqrt(nout)*(Vout.conj().T@full[q])
        if np.linalg.norm(Y)>tol:out[q]=Y
    return out


def compare_s4(qin,cols=2):
    d=s4_triv(tuple(qin)).shape[1]
    k=min(cols,d);X0=np.eye(d,k,dtype=complex)
    a=reduced_s4_triv_to_sign(qin,X0);b=explicit_s4_triv_to_sign(qin,X0)
    keys=set(a)|set(b);num=den=0.;mx=0.
    for q in keys:
        aa=a.get(q);bb=b.get(q)
        if aa is None:aa=np.zeros_like(bb)
        if bb is None:bb=np.zeros_like(aa)
        dd=np.linalg.norm(aa-bb);nn=np.linalg.norm(bb)
        num+=dd*dd;den+=nn*nn
        if nn>1e-10:mx=max(mx,dd/nn)
    return {'qin':list(qin),'d':d,'k':k,'keys':len(keys),'rel':math.sqrt(num/max(den,1e-60)),'mx':mx}


def H0_41_block(rep,C=None,tol=1e-11,selected=None):
    rep=tuple(rep);bases,slices,C0=complement41(rep)
    if C is None:C=C0
    else:C=np.asarray(C,complex)
    if C.shape[0]!=C0.shape[0]:raise ValueError('C rows')
    acc={}
    for b in bases:
        seg=C[slices[b],:]
        if np.linalg.norm(seg)<tol:continue
        out=reduced_s4_triv_to_sign(b,seg,tol,selected)
        for q,Y in out.items():acc[q]=acc.get(q,0)+Y
    return {q:Y for q,Y in acc.items() if np.linalg.norm(Y)>tol}


def probe41(rep,k=2,seed=0,tol=1e-10):
    rep=tuple(rep);bases,slices,E=s5_trivial_embedding_in_s4(rep);d=E.shape[0];m=d-E.shape[1]
    if m<=0:return {},m
    rng=np.random.default_rng(seed)
    R=rng.standard_normal((d,k))+1j*rng.standard_normal((d,k))
    if E.shape[1]:R-=E@(E.conj().T@R)
    # orthonormalize probes within complement
    Q,_=np.linalg.qr(R,mode='reduced');Q=Q[:,:min(k,m)]
    acc={}
    for b in bases:
        seg=Q[slices[b],:]
        if np.linalg.norm(seg)<tol:continue
        out=reduced_s4_triv_to_sign(b,seg,tol)
        for q,Y in out.items():acc[q]=acc.get(q,0)+Y
    return {q:Y for q,Y in acc.items() if np.linalg.norm(Y)>tol},m

# --- targeted selected-q path for numeric peeling ---
@functools.lru_cache(None)
def s4_images(q):
    q=tuple(q)
    return tuple(sorted({S.permute_spins(q,p) for p in S4}))


def reduced_s4_triv_to_sign_selected(qin,Xcoord,selected,tol=1e-11):
    import bqg_targeted_h as T
    import bqg_safe_boundary_master_gate as R
    qin=tuple(qin); Xcoord=np.asarray(Xcoord,complex)
    if Xcoord.ndim==1:Xcoord=Xcoord[:,None]
    selected=tuple(map(tuple,selected)); sel=set(selected)
    Vin=s4_triv(qin)
    if Xcoord.shape[0]!=Vin.shape[1]:raise ValueError(('rows',qin,Xcoord.shape,Vin.shape))
    B=Vin@Xcoord
    # only concrete S4 images compatible with one of the six node-0 triangles
    targets=set(); pairs=tuple(itertools.combinations(R.NEIG[0],2));Jmax=max(qin)+1
    for q in selected:
        for tt in s4_images(q):
            if any(T.compatible_triangle_target(qin,tt,0,a,b,Jmax) for a,b in pairs):targets.add(tt)
    if not targets:return {}
    H=T.H_node_apply_targets(qin,0,B,targets,tol)
    tmp={}
    for so,W in H.items():
        q,p=canon4(tuple(so))
        if q not in sel:continue
        tq,TW=X.transform_coeff(tuple(so),tuple(p),W); assert tq==q
        TW=perm_sign(p)*TW
        tmp[q]=tmp.get(q,0)+TW
    gin=len(s4_stabilizer(qin));out={}
    for q,Z in tmp.items():
        Vout=s4_sign(q)
        if Vout.shape[1]==0:continue
        gout=len(s4_stabilizer(q))
        Y=math.sqrt(gout/gin)*(Vout.conj().T@Z)
        if np.linalg.norm(Y)>tol:out[q]=Y
    return out


def H0_41_block_selected(rep,selected,C=None,tol=1e-11):
    rep=tuple(rep); selected=tuple(map(tuple,selected));bases,slices,C0=complement41(rep)
    if C is None:C=C0
    else:C=np.asarray(C,complex)
    if C.shape[0]!=C0.shape[0]:raise ValueError('C rows')
    acc={}
    for b in bases:
        seg=C[slices[b],:]
        if np.linalg.norm(seg)<tol:continue
        out=reduced_s4_triv_to_sign_selected(b,seg,selected,tol)
        for q,Y in out.items():acc[q]=acc.get(q,0)+Y
    return {q:Y for q,Y in acc.items() if np.linalg.norm(Y)>tol}

# --- alternative [4,1] branch: S4 [3,1] -> choose its S3-trivial line ---
S3=tuple(p for p in S4 if p[1]==1)  # fixes vertices 0 and 1; permutes 2,3,4
G3=len(S3)
T01=(1,0,2,3,4)

@functools.lru_cache(None)
def s3_stabilizer(spins):
    spins=tuple(spins); return tuple(p for p in S3 if S.permute_spins(spins,p)==spins)
@functools.lru_cache(None)
def canon3(spins):
    spins=tuple(spins);best=None;bp=None
    for p in S3:
        t=S.permute_spins(spins,p)
        if best is None or t<best:best=t;bp=p
    return best,bp
@functools.lru_cache(None)
def s3_triv(spins):return _projector_basis(tuple(spins),s3_stabilizer(tuple(spins)),lambda p:1)
@functools.lru_cache(None)
def s3_sign(spins):return _projector_basis(tuple(spins),s3_stabilizer(tuple(spins)),perm_sign)
@functools.lru_cache(None)
def s3_images(q):
    q=tuple(q);return tuple(sorted({S.permute_spins(q,p) for p in S3}))


def expand_s4_coords(rep,C):
    rep=tuple(rep);bases,slices,E=s5_trivial_embedding_in_s4(rep);C=np.asarray(C,complex);out={}
    for b in bases:
        seg=C[slices[b],:]
        if not seg.size or np.linalg.norm(seg)<1e-13:continue
        B=s4_triv(b)@seg;n4=24/len(s4_stabilizer(b));scale=1/math.sqrt(n4)
        seen={}
        for p in S4:
            ts,TB=X.transform_coeff(b,p,B)
            if ts not in seen:seen[ts]=TB
        if len(seen)!=int(n4):raise RuntimeError(('s4 orbit size',b,len(seen),n4))
        for ts,TB in seen.items():out[ts]=out.get(ts,0)+scale*TB
    return out


@functools.lru_cache(None)
def branch41_s3_coords(rep):
    """Normalized S3-invariant line inside the S4-[3,1] branch of each [4,1] multiplicity copy.
    Returns S3 canonical bases, slices, W coordinate matrix (d3 x m41)."""
    rep=tuple(rep);bases4,slices4,C=complement41(rep);m=C.shape[1]
    if m==0:return (),{},np.zeros((0,0),complex)
    v=expand_s4_coords(rep,C);tv={}
    for sp,B in v.items():
        ts,TB=X.transform_coeff(sp,T01,B);tv[ts]=tv.get(ts,0)+TB
    keys=set(v)|set(tv);w={}
    fac=4/math.sqrt(15)
    for sp in keys:
        a=tv.get(sp);b=v.get(sp)
        if a is None:a=np.zeros_like(b)
        if b is None:b=np.zeros_like(a)
        z=fac*(a+0.25*b)
        if np.linalg.norm(z)>1e-12:w[sp]=z
    # reduce normalized invariant state to normalized S3 orbit coordinates
    bases3=tuple(sorted({canon3(sp)[0] for sp in w}));parts=[];slices={};off=0
    for q in bases3:
        Bq=w.get(q)
        if Bq is None:raise RuntimeError(('missing canonical s3',q))
        V3=s3_triv(q);n3=6/len(s3_stabilizer(q));Y=math.sqrt(n3)*(V3.conj().T@Bq)
        d=Y.shape[0];slices[q]=slice(off,off+d);off+=d;parts.append(Y)
    W=np.vstack(parts) if parts else np.zeros((0,m),complex)
    if m:
        e=np.linalg.norm(W.conj().T@W-np.eye(m))
        if e>3e-8:raise RuntimeError(('W41s3 orth',rep,W.shape,e))
    return bases3,slices,W


def reduced_s3_triv_to_sign_selected(qin,Xcoord,selected,tol=1e-11):
    import bqg_targeted_h as T
    import bqg_safe_boundary_master_gate as R
    qin=tuple(qin);Xcoord=np.asarray(Xcoord,complex)
    if Xcoord.ndim==1:Xcoord=Xcoord[:,None]
    selected=tuple(map(tuple,selected));sel=set(selected);Vin=s3_triv(qin);B=Vin@Xcoord
    pairs=tuple(itertools.combinations(R.NEIG[0],2));Jmax=max(qin)+1;targets=set()
    for q in selected:
        for tt in s3_images(q):
            if any(T.compatible_triangle_target(qin,tt,0,a,b,Jmax) for a,b in pairs):targets.add(tt)
    if not targets:return {}
    H=T.H_node_apply_targets(qin,0,B,targets,tol);tmp={}
    for so,M in H.items():
        q,p=canon3(tuple(so))
        if q not in sel:continue
        tq,TM=X.transform_coeff(tuple(so),tuple(p),M);assert tq==q
        tmp[q]=tmp.get(q,0)+perm_sign(p)*TM
    gin=len(s3_stabilizer(qin));out={}
    for q,Z in tmp.items():
        V=s3_sign(q)
        if V.shape[1]==0:continue
        gout=len(s3_stabilizer(q));Y=math.sqrt(gout/gin)*(V.conj().T@Z)
        if np.linalg.norm(Y)>tol:out[q]=Y
    return out


def H0_41_s3_selected(rep,selected,W=None,tol=1e-11):
    rep=tuple(rep);bases,slices,W0=branch41_s3_coords(rep)
    if W is None:W=W0
    acc={}
    for b in bases:
        seg=W[slices[b],:]
        if np.linalg.norm(seg)<tol:continue
        out=reduced_s3_triv_to_sign_selected(b,seg,selected,tol)
        for q,Y in out.items():acc[q]=acc.get(q,0)+Y
    return {q:Y for q,Y in acc.items() if np.linalg.norm(Y)>tol}

# --- dual S3-invariant branch pair inside S5 [4,1] ---
def _raw_invariant_to_s3_dict(raw,m,tol=1e-11):
    """Convert an explicitly expanded S3-invariant raw spin-block state to normalized S3-orbit coordinates."""
    qs=tuple(sorted({canon3(tuple(sp))[0] for sp in raw}))
    out={}
    for q in qs:
        B=raw.get(q)
        if B is None: raise RuntimeError(('missing s3 canonical raw',q))
        V=s3_triv(q); n3=G3/len(s3_stabilizer(q)); Y=math.sqrt(n3)*(V.conj().T@B)
        if np.linalg.norm(Y)>tol: out[q]=Y
    return out

@functools.lru_cache(None)
def branch41_dual_s3_dicts(rep):
    """Two orthonormal S3-invariant branch lines v,w for every [4,1] multiplicity copy.
    v is the S4-[4] line; w is the S4-[3,1] S3-trivial line."""
    rep=tuple(rep);bases,slices,C=complement41(rep);m=C.shape[1]
    if m==0:return {},{},0
    v=expand_s4_coords(rep,C)
    tv={}
    for sp,B in v.items():
        ts,TB=X.transform_coeff(sp,T01,B);tv[ts]=tv.get(ts,0)+TB
    keys=set(v)|set(tv);w={};fac=4/math.sqrt(15)
    for sp in keys:
        a=tv.get(sp);b=v.get(sp)
        if a is None:a=np.zeros_like(b)
        if b is None:b=np.zeros_like(a)
        z=fac*(a+0.25*b)
        if np.linalg.norm(z)>1e-12:w[sp]=z
    Vd=_raw_invariant_to_s3_dict(v,m); Wd=_raw_invariant_to_s3_dict(w,m)
    # Gram checks in normalized S3 coordinate direct sum.
    def gram(d1,d2):
        g=np.zeros((m,m),complex)
        for q in set(d1)|set(d2):
            a=d1.get(q);b=d2.get(q)
            if a is None:
                # infer row count from b
                a=np.zeros_like(b)
            if b is None:b=np.zeros_like(a)
            g+=a.conj().T@b
        return g
    e1=np.linalg.norm(gram(Vd,Vd)-np.eye(m));e2=np.linalg.norm(gram(Wd,Wd)-np.eye(m));ec=np.linalg.norm(gram(Vd,Wd))
    if max(e1,e2,ec)>5e-8:raise RuntimeError(('dual s3 gram',rep,m,e1,e2,ec))
    return Vd,Wd,m

def H0_41_dual_s3_selected(rep,selected,tol=1e-11):
    """Stack H0 on the two orthonormal S3-invariant [4,1] branch vectors.
    Returns q -> matrix with 2*dim(S3-sign(q)) rows and m multiplicity columns."""
    rep=tuple(rep);selected=tuple(map(tuple,selected));Vd,Wd,m=branch41_dual_s3_dicts(rep)
    accV={};accW={}
    for b in sorted(set(Vd)|set(Wd)):
        a=Vd.get(b);w=Wd.get(b)
        if a is None:a=np.zeros_like(w)
        if w is None:w=np.zeros_like(a)
        Xpair=np.hstack([a,w])
        if np.linalg.norm(Xpair)<tol:continue
        out=reduced_s3_triv_to_sign_selected(b,Xpair,selected,tol)
        for q,Y in out.items():
            yv=Y[:,:m];yw=Y[:,m:]
            accV[q]=accV.get(q,0)+yv;accW[q]=accW.get(q,0)+yw
    out={}
    for q in selected:
        a=accV.get(q);w=accW.get(q)
        if a is None and w is None:continue
        if a is None:a=np.zeros_like(w)
        if w is None:w=np.zeros_like(a)
        Y=np.vstack([a,w])
        if np.linalg.norm(Y)>tol:out[q]=Y
    return out

