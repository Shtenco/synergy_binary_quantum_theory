#!/usr/bin/env python3
"""First regulator-safe Peter-Weyl K5 HH column with genuine local volume.

This verifier implements the finite link cutoff
    H_link = direct_sum_{j<=Jmax} V_j^L tensor V_j^R
with doubled integer spins, exact SU(2) Clebsch-Gordan recoupling and the
four-valent volume V=sqrt(|J1 . (J2 x J3)|).

It evaluates one orientation-covariant Hermitian node-Hamiltonian commutator
[H0,H1] on the all-j=1/2, K_v=0 K5 boundary state at Jmax=5/2, and also checks
whether the old j=1/2 five-tetrahedron BF vertex V5 remains in the kernel.

This is a first safe-column falsifier, not the full quantum HDA closure test.
The exact SymPy/CG implementation is intentionally reference-quality rather
than fast; the full run can be slow on a cold cache.
"""
from __future__ import annotations
import argparse, functools, itertools, json, math
from pathlib import Path
import numpy as np
import sympy as sp
from sympy.physics.wigner import clebsch_gordan
import bqg_pruning_audit as PRUNE

VERT=tuple(range(5))
EDGES=tuple(itertools.combinations(VERT,2))
EIDX={e:i for i,e in enumerate(EDGES)}
NEIG={v:tuple(w for w in VERT if w!=v) for v in VERT}
LEGIDX={(v,w):NEIG[v].index(w) for v in VERT for w in NEIG[v]}

EPS3=np.zeros((3,3,3),int)
for p in itertools.permutations(range(3)):
    inv=sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
    EPS3[p]=(-1)**inv
EPS2=np.array([[0,1],[-1,0]],complex)

@PRUNE.audit_lru_cache(None)
def cg2(s1,s2,S,m1,m2,M):
    vals=[sp.Rational(x,2) for x in (s1,s2,S,m1,m2,M)]
    return float(clebsch_gordan(*vals).evalf())

@PRUNE.audit_lru_cache(None)
def m2vals_t(s):
    return tuple(range(s,-s-1,-2))

@PRUNE.audit_lru_cache(None)
def allowed_k2_t(s1,s2,s3,s4):
    return tuple(sorted(set(range(abs(s1-s2),s1+s2+1,2)) &
                        set(range(abs(s3-s4),s3+s4+1,2))))

@PRUNE.audit_lru_cache(None)
def intertwiner_tensor_cached(spins,K):
    spins=tuple(spins)
    dims=[s+1 for s in spins]
    T=np.zeros(dims,dtype=complex)
    mvs=[m2vals_t(s) for s in spins]
    for idxs in itertools.product(*[range(d) for d in dims]):
        ms=[mvs[a][idxs[a]] for a in range(4)]
        val=0.0
        for M in m2vals_t(K):
            phase=(-1)**((K-M)//2)
            val += (cg2(spins[0],spins[1],K,ms[0],ms[1],M)
                    *cg2(spins[2],spins[3],K,ms[2],ms[3],-M)
                    *phase/math.sqrt(K+1))
        T[idxs]=val
    return T

@PRUNE.audit_lru_cache(None)
def spin_mats_cached(s):
    j=s/2
    ms=np.array(m2vals_t(s),float)/2
    d=s+1
    Jz=np.diag(ms).astype(complex)
    Jp=np.zeros((d,d),complex)
    for col,m in enumerate(ms):
        mp=m+1
        arr=np.where(np.isclose(ms,mp))[0]
        if len(arr):
            Jp[arr[0],col]=math.sqrt(max(0,j*(j+1)-m*(m+1)))
    Jm=Jp.conj().T
    return ((Jp+Jm)/2,(Jp-Jm)/(2j) if j>0 else np.zeros((1,1),complex),Jz)

def apply_axis_np(T,axis,M):
    return np.moveaxis(np.tensordot(M,np.moveaxis(T,axis,0),axes=(1,0)),0,axis)

def local_spins(spins,v):
    return tuple(spins[EIDX[tuple(sorted((v,w)))]] for w in NEIG[v])

@PRUNE.audit_lru_cache(None)
def cg_leg_matrix(s_in,s_out,a2):
    Mi=m2vals_t(s_in); Mo=m2vals_t(s_out)
    M=np.zeros((s_out+1,s_in+1),complex)
    for ci,m in enumerate(Mi):
        for ro,mo in enumerate(Mo):
            M[ro,ci]=cg2(s_in,1,s_out,m,a2,mo)
    return M

def hit_mats(s_in,s_out,x,y,i,j):
    if x<y:
        Ml=cg_leg_matrix(s_in,s_out,1 if i==0 else -1)
        Mr=cg_leg_matrix(s_in,s_out,1 if j==0 else -1)
        norm=math.sqrt((s_in+1)/(s_out+1))
    else:
        Ml=cg_leg_matrix(s_out,s_in,1 if j==0 else -1).conj().T
        Mr=cg_leg_matrix(s_out,s_in,1 if i==0 else -1).conj().T
        norm=math.sqrt((s_out+1)/(s_in+1))
    return Ml,Mr,norm

@PRUNE.audit_lru_cache(None)
def epsilon_j(s):
    ms=m2vals_t(s);d=s+1
    E=np.zeros((d,d),complex)
    idx={m:i for i,m in enumerate(ms)}
    for m in ms:
        E[idx[m],idx[-m]]=(-1)**((s-m)//2)
    return E

@PRUNE.audit_lru_cache(None)
def oriented_intertwiner(v,spins_local,K):
    T=intertwiner_tensor_cached(tuple(spins_local),K).copy()
    for leg,w in enumerate(NEIG[v]):
        if w<v:
            T=apply_axis_np(T,leg,epsilon_j(spins_local[leg]))
    return T

def initial_factorized_oriented(key):
    spins,Ks=key
    return (spins,tuple(oriented_intertwiner(v,local_spins(spins,v),Ks[v]).copy()
                        for v in VERT),1+0j)

def apply_hit_branch(branch,x,y,i,j,Jmax2):
    spins,tensors,amp=branch
    ei=EIDX[tuple(sorted((x,y)))]
    s=spins[ei]
    out=[]
    for so in (s-1,s+1):
        if so<0 or so>Jmax2: continue
        Ml,Mr,norm=hit_mats(s,so,x,y,i,j)
        u,v=sorted((x,y))
        tens=list(tensors)
        tens[u]=apply_axis_np(tens[u],LEGIDX[(u,v)],Ml)
        tens[v]=apply_axis_np(tens[v],LEGIDX[(v,u)],Mr)
        spn=list(spins);spn[ei]=so;spn=tuple(spn)
        out.append((spn,tuple(tens),amp*norm))
    return out

def apply_path_branch(branch,path,iout,iin,Jmax2):
    m=len(path)-1
    outs=[]
    for mids in itertools.product(range(2),repeat=max(0,m-1)):
        cols=(iout,)+tuple(mids)+(iin,)
        branches=[branch]
        for r,(x,y) in enumerate(zip(path[:-1],path[1:])):
            nb=[]
            for br in branches:
                nb.extend(apply_hit_branch(br,x,y,cols[r],cols[r+1],Jmax2))
            branches=nb
            if not branches: break
        outs.extend(branches)
    return outs

@PRUNE.audit_lru_cache(None)
def volume123_matrix(s1,s2,s3):
    mats=[spin_mats_cached(s) for s in (s1,s2,s3)]
    d=(s1+1)*(s2+1)*(s3+1)
    Q=np.zeros((d,d),complex)
    for a,b,c in itertools.product(range(3),repeat=3):
        e=EPS3[a,b,c]
        if e:
            Q += e*np.kron(np.kron(mats[0][a],mats[1][b]),mats[2][c])
    Q=(Q+Q.conj().T)/2
    ev,U=np.linalg.eigh(Q)
    return (U*np.sqrt(np.abs(ev)))@U.conj().T

def apply_volume_tensor(T,spins_local):
    d1,d2,d3,d4=[s+1 for s in spins_local]
    V=volume123_matrix(spins_local[0],spins_local[1],spins_local[2])
    A=T.reshape(d1*d2*d3,d4)
    return (V@A).reshape(d1,d2,d3,d4)

def apply_volume_tensor_oriented(T,spins_local,v):
    X=T
    for leg,w in enumerate(NEIG[v]):
        if w<v:
            X=apply_axis_np(X,leg,epsilon_j(spins_local[leg]).conj().T)
    X=apply_volume_tensor(X,spins_local)
    for leg,w in enumerate(NEIG[v]):
        if w<v:
            X=apply_axis_np(X,leg,epsilon_j(spins_local[leg]))
    return X

def add_dict(dst,src,scale=1):
    for k,a in src.items():
        prior=dst.get(k,0j)
        incoming=scale*a
        v=prior+incoming
        PRUNE.record_merge("k5.add_dict.1e-13",k,prior,incoming,v,1e-13)
        if PRUNE.keep(v,1e-13): dst[k]=v
        elif k in dst: del dst[k]

def adjoint_path(path,iout,iin):
    return tuple(reversed(path)), iin, iout

def adjoint_sequence(seq):
    out=[]
    for op in reversed(seq):
        if op[0]=='V': out.append(op)
        else:
            p,io,ii=adjoint_path(op[1],op[2],op[3])
            out.append(('P',p,io,ii))
    return tuple(out)

def T_sequences(v,a,b,c):
    seqs=[]
    for i,j,k in itertools.product(range(2),repeat=3):
        base1=(('V',v),('P',(c,v),k,i),('P',(v,c),j,k))
        base2=(('P',(c,v),k,i),('V',v),('P',(v,c),j,k))
        pf=('P',(v,a,b,v),i,j)
        pr=('P',(v,b,a,v),i,j)
        seqs += [(+1,base1+(pf,)),(-1,base1+(pr,)),
                 (-1,base2+(pf,)),(+1,base2+(pr,))]
    return tuple(seqs)

@PRUNE.audit_lru_cache(None)
def oriented_specs(v):
    neigh=NEIG[v]
    specs=[]
    for r in range(4):
        tri=tuple(neigh[i] for i in range(4) if i!=r)
        sign=(-1)**r
        a,b,c=tri
        specs += [(sign,(v,a,b,c)),(sign,(v,b,c,a)),(sign,(v,c,a,b))]
    return tuple(specs)

def apply_primitive_oriented(input_key,seq,Jmax2):
    branches=[initial_factorized_oriented(input_key)]
    for op in seq:
        if op[0]=='V':
            v=op[1];nb=[]
            for spins,tensors,amp in branches:
                t=list(tensors)
                t[v]=apply_volume_tensor_oriented(t[v],local_spins(spins,v),v)
                nb.append((spins,tuple(t),amp))
            branches=nb
        else:
            nb=[]
            for br in branches:
                nb.extend(apply_path_branch(br,op[1],op[2],op[3],Jmax2))
            branches=nb
        if not branches: break
    out={}
    for spins,tensors,amp in branches:
        local_opts=[];ok=True
        for v in VERT:
            ls=local_spins(spins,v);opts=[]
            for K in allowed_k2_t(*ls):
                c=np.vdot(oriented_intertwiner(v,ls,K),tensors[v])
                PRUNE.record_threshold("k5.primitive.intertwiner_projection.1e-12",(spins,v,K),c,1e-12,incoming=amp)
                if PRUNE.keep(c,1e-12): opts.append((K,c))
            if not opts: ok=False;break
            local_opts.append(opts)
        if not ok: continue
        for ch in itertools.product(*local_opts):
            val=amp
            for _,c in ch: val*=c
            ko=(spins,tuple(k for k,_ in ch))
            PRUNE.record_threshold("k5.primitive.product_coefficient.1e-12",ko,val,1e-12)
            if PRUNE.keep(val,1e-12):
                out[ko]=out.get(ko,0j)+val
    for k,v in out.items():
        if abs(v)<=1e-11: PRUNE.record_threshold("k5.primitive.output.1e-11",k,v,1e-11)
    return {k:v for k,v in out.items() if PRUNE.keep(v,1e-11)}

def apply_T_oriented_basis(input_key,v,a,b,c,Jmax2,adj=False):
    out={}
    for coef,seq in T_sequences(v,a,b,c):
        if adj: seq=adjoint_sequence(seq)
        add_dict(out,apply_primitive_oriented(input_key,seq,Jmax2),coef)
    return out

@PRUNE.audit_lru_cache(None)
def T_cached(key,v,a,b,c,Jmax2,adj):
    return tuple(apply_T_oriented_basis(key,v,a,b,c,Jmax2,adj).items())

def apply_T_cached_state(state,spec,Jmax2,adj=False):
    out={}
    for key,amp0 in state.items():
        for ko,c in T_cached(key,*spec,Jmax2,adj):
            out[ko]=out.get(ko,0j)+amp0*c
    for k,v in out.items():
        if abs(v)<=1e-10: PRUNE.record_threshold("k5.cached_T_state.1e-10",k,v,1e-10)
    return {k:v for k,v in out.items() if PRUNE.keep(v,1e-10)}

def apply_H_cached_state(state,v,Jmax2):
    out={}
    for sign,spec in oriented_specs(v):
        rr=apply_T_cached_state(state,spec,Jmax2,False)
        aa=apply_T_cached_state(state,spec,Jmax2,True)
        add_dict(out,rr,0.5*sign);add_dict(out,aa,0.5*sign)
    return out

def prune_state(st,tol=1e-8):
    for k,v in st.items():
        if abs(v)<=tol: PRUNE.record_threshold("k5.prune_state",k,v,tol)
    return {k:v for k,v in st.items() if PRUNE.keep(v,tol)}

def prev_intertwiners():
    z=np.array([1,0],complex);o=np.array([0,1],complex)
    s=(np.kron(z,o)-np.kron(o,z))/math.sqrt(2)
    i0=np.kron(s,s)
    tp=np.kron(z,z);t0=(np.kron(z,o)+np.kron(o,z))/math.sqrt(2);tm=np.kron(o,o)
    i1=(np.kron(tp,tm)-np.kron(t0,t0)+np.kron(tm,tp))/math.sqrt(3)
    return [i0.reshape(2,2,2,2),i1.reshape(2,2,2,2)]

def apply_axis_simple(T,axis,M):
    A=np.moveaxis(T,axis,0);B=np.tensordot(M,A,axes=(1,0));return np.moveaxis(B,0,axis)

def node_tensor_v5(v,iota,I):
    T=I[iota].copy();neighbors=[w for w in range(5) if w!=v]
    for ax,w in enumerate(neighbors):
        if w<v: T=apply_axis_simple(T,ax,EPS2)
    return T

def v5_tensor():
    I=prev_intertwiners();V=np.zeros((2,)*5,complex)
    for io in itertools.product(range(2),repeat=5):
        T=[node_tensor_v5(v,io[v],I) for v in range(5)]
        V[io]=np.einsum('abcd,aefg,behi,cfhj,dgij->',*T,optimize=True)
    return V.reshape(-1)

def basis_full_jhalf():
    spins=(1,)*10
    return [(spins,tuple(bits)) for bits in itertools.product((0,2),repeat=5)]

def norm2_state(state):
    return float(sum(abs(v)**2 for v in state.values()))

def compose_on_sparse(state,target_v,Jmax2):
    out={}
    for key,amp in state.items():
        rr=apply_H_cached_state({key:1+0j},target_v,Jmax2)
        for ko,c in rr.items(): out[ko]=out.get(ko,0j)+amp*c
    return prune_state(out,1e-8)

def run():
    JMAX2=5
    initial=basis_full_jhalf()[0]
    psi0={initial:1+0j}
    h0=prune_state(apply_H_cached_state(psi0,0,JMAX2),1e-8)
    h1=prune_state(apply_H_cached_state(psi0,1,JMAX2),1e-8)
    h1h0=compose_on_sparse(h0,1,JMAX2)
    h0h1=compose_on_sparse(h1,0,JMAX2)
    comm={};add_dict(comm,h0h1,+1);add_dict(comm,h1h0,-1);comm=prune_state(comm,1e-8)
    total2=norm2_state(comm)
    fixed={k:v for k,v in comm.items() if k[0]==(1,)*10}
    fixed2=norm2_state(fixed)
    max_spin2=max(max(k[0]) for k in comm) if comm else 0

    V5=v5_tensor();v5norm=float(np.linalg.norm(V5));full=basis_full_jhalf()
    coeffs=[]
    for idx,val in enumerate(V5):
        coefficient=val/v5norm
        # Preserve the historical production decision on the raw V5 entry;
        # normalization is applied only to the retained coefficient.
        PRUNE.record_threshold("k5.V5_input_decomposition.1e-12",full[idx],val,1e-12)
        if PRUNE.keep(val,1e-12):
            coeffs.append((full[idx],coefficient))
    h0v={}
    for key,amp in coeffs:
        rr=apply_H_cached_state({key:1+0j},0,JMAX2)
        for ko,c in rr.items(): h0v[ko]=h0v.get(ko,0j)+amp*c
    h0v=prune_state(h0v,1e-8);h0v2=norm2_state(h0v)
    h0v_fixed={k:v for k,v in h0v.items() if k[0]==(1,)*10}
    h0v_max_spin2=max(max(k[0]) for k in h0v) if h0v else 0

    result={
      'status':'first regulator-safe Peter-Weyl K5 HH column','Jmax':2.5,
      'input':'all ten links j=1/2; all five recoupling labels K=0',
      'H0_support':len(h0),'H1_support':len(h1),'H1H0_support':len(h1h0),'H0H1_support':len(h0h1),
      'commutator_support':len(comm),'commutator_norm':math.sqrt(total2),'commutator_norm_squared':total2,
      'fixed_all_jhalf_support':len(fixed),'fixed_all_jhalf_norm':math.sqrt(fixed2),
      'fixed_all_jhalf_norm_squared_fraction':fixed2/total2 if total2 else 0.0,
      'outside_all_jhalf_norm_squared_fraction':1-fixed2/total2 if total2 else 0.0,
      'outside_all_jhalf_amplitude_fraction':math.sqrt(max(0.0,1-fixed2/total2)) if total2 else 0.0,
      'max_spin_reached_in_commutator':max_spin2/2,
      'old_V5_nonzero_components':len(coeffs),'old_V5_norm_squared':float(np.vdot(V5,V5).real),
      'safe_H0_on_old_V5_norm':math.sqrt(h0v2),'safe_H0_on_old_V5_support':len(h0v),
      'safe_H0_on_old_V5_fixed_all_jhalf_support':len(h0v_fixed),
      'safe_H0_on_old_V5_fixed_all_jhalf_norm':math.sqrt(norm2_state(h0v_fixed)),
      'safe_H0_on_old_V5_max_spin_reached':h0v_max_spin2/2,
      'scope_note':'One safe HH column only; no discrete D(k,l), simplicity or Lorentzian correction yet.'
    }
    targets={'commutator_norm':1.681559985798016,'fixed_fraction':0.29790166313739946,'H0V5_norm':1.4002194669856702}
    result['regression_errors']={
      'commutator_norm':abs(result['commutator_norm']-targets['commutator_norm']),
      'fixed_fraction':abs(result['fixed_all_jhalf_norm_squared_fraction']-targets['fixed_fraction']),
      'H0V5_norm':abs(result['safe_H0_on_old_V5_norm']-targets['H0V5_norm'])}
    result['passed_regression']=(max(result['regression_errors'].values())<5e-8 and result['safe_H0_on_old_V5_fixed_all_jhalf_support']==0)
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path);a=ap.parse_args()
    out=run();txt=json.dumps(out,indent=2);print(txt)
    if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(txt+'\n',encoding='utf-8')
    return 0 if out['passed_regression'] else 1

if __name__=='__main__':raise SystemExit(main())
