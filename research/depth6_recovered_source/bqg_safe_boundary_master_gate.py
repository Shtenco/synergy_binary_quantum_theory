#!/usr/bin/env python3
from __future__ import annotations
import functools, itertools, math, json, argparse
from pathlib import Path
import numpy as np
import sympy as sp
from sympy.physics.wigner import clebsch_gordan

# --- small tested sparse-master utilities ---
def sparse_inner(a: dict, b: dict):
    if len(a) > len(b):
        a,b=b,a
    return sum(np.conj(v)*b.get(k,0j) for k,v in a.items())

def gram_from_sparse_columns(cols):
    n=len(cols); G=np.zeros((n,n),complex)
    for i in range(n):
        G[i,i]=sparse_inner(cols[i],cols[i])
        for j in range(i+1,n):
            z=sparse_inner(cols[i],cols[j])
            G[i,j]=z; G[j,i]=np.conj(z)
    return 0.5*(G+G.conj().T)

def kernel_diagnostics(M, rtol=1e-10):
    H=0.5*(M+M.conj().T)
    ev,U=np.linalg.eigh(H)
    scale=max(float(np.max(np.abs(ev))),1.0)
    tol=rtol*scale
    zero=np.abs(ev)<=tol
    pos=ev[ev>tol]
    return {
        'rank':int(np.sum(~zero)),
        'nullity':int(np.sum(zero)),
        'tolerance':float(tol),
        'eigenvalues':[float(x) for x in ev],
        'min_eigenvalue':float(ev[0]),
        'min_positive':float(np.min(pos)) if len(pos) else None,
        'kernel_vectors':U[:,zero],
    }

# --- exact repository-safe Peter-Weyl machinery used by k5_peter_weyl_safe_hda_column.py ---
VERT=tuple(range(5))
EDGES=tuple(itertools.combinations(VERT,2))
EIDX={e:i for i,e in enumerate(EDGES)}
NEIG={v:tuple(w for w in VERT if w!=v) for v in VERT}
LEGIDX={(v,w):NEIG[v].index(w) for v in VERT for w in NEIG[v]}
EPS3=np.zeros((3,3,3),int)
for p in itertools.permutations(range(3)):
    inv=sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
    EPS3[p]=(-1)**inv

@functools.lru_cache(None)
def cg2(s1,s2,S,m1,m2,M):
    vals=[sp.Rational(x,2) for x in (s1,s2,S,m1,m2,M)]
    return float(clebsch_gordan(*vals).evalf())

@functools.lru_cache(None)
def m2vals_t(s): return tuple(range(s,-s-1,-2))

@functools.lru_cache(None)
def allowed_k2_t(s1,s2,s3,s4):
    return tuple(sorted(set(range(abs(s1-s2),s1+s2+1,2)) & set(range(abs(s3-s4),s3+s4+1,2))))

@functools.lru_cache(None)
def intertwiner_tensor_cached(spins,K):
    spins=tuple(spins); dims=[s+1 for s in spins]
    T=np.zeros(dims,dtype=complex); mvs=[m2vals_t(s) for s in spins]
    for idxs in itertools.product(*[range(d) for d in dims]):
        ms=[mvs[a][idxs[a]] for a in range(4)]; val=0.0
        for M in m2vals_t(K):
            phase=(-1)**((K-M)//2)
            val += (cg2(spins[0],spins[1],K,ms[0],ms[1],M)*
                    cg2(spins[2],spins[3],K,ms[2],ms[3],-M)*phase/math.sqrt(K+1))
        T[idxs]=val
    return T

@functools.lru_cache(None)
def spin_mats_cached(s):
    j=s/2; ms=np.array(m2vals_t(s),float)/2; d=s+1
    Jz=np.diag(ms).astype(complex); Jp=np.zeros((d,d),complex)
    for col,m in enumerate(ms):
        arr=np.where(np.isclose(ms,m+1))[0]
        if len(arr): Jp[arr[0],col]=math.sqrt(max(0,j*(j+1)-m*(m+1)))
    Jm=Jp.conj().T
    return ((Jp+Jm)/2,(Jp-Jm)/(2j) if j>0 else np.zeros((1,1),complex),Jz)

def apply_axis_np(T,axis,M):
    return np.moveaxis(np.tensordot(M,np.moveaxis(T,axis,0),axes=(1,0)),0,axis)

def local_spins(spins,v):
    return tuple(spins[EIDX[tuple(sorted((v,w)))]] for w in NEIG[v])

@functools.lru_cache(None)
def cg_leg_matrix(s_in,s_out,a2):
    Mi=m2vals_t(s_in); Mo=m2vals_t(s_out); M=np.zeros((s_out+1,s_in+1),complex)
    for ci,m in enumerate(Mi):
        for ro,mo in enumerate(Mo): M[ro,ci]=cg2(s_in,1,s_out,m,a2,mo)
    return M

def hit_mats(s_in,s_out,x,y,i,j):
    if x<y:
        Ml=cg_leg_matrix(s_in,s_out,1 if i==0 else -1); Mr=cg_leg_matrix(s_in,s_out,1 if j==0 else -1)
        norm=math.sqrt((s_in+1)/(s_out+1))
    else:
        Ml=cg_leg_matrix(s_out,s_in,1 if j==0 else -1).conj().T; Mr=cg_leg_matrix(s_out,s_in,1 if i==0 else -1).conj().T
        norm=math.sqrt((s_out+1)/(s_in+1))
    return Ml,Mr,norm

@functools.lru_cache(None)
def epsilon_j(s):
    ms=m2vals_t(s); d=s+1; E=np.zeros((d,d),complex); idx={m:i for i,m in enumerate(ms)}
    for m in ms: E[idx[m],idx[-m]]=(-1)**((s-m)//2)
    return E

@functools.lru_cache(None)
def oriented_intertwiner(v,spins_local,K):
    T=intertwiner_tensor_cached(tuple(spins_local),K).copy()
    for leg,w in enumerate(NEIG[v]):
        if w<v: T=apply_axis_np(T,leg,epsilon_j(spins_local[leg]))
    return T

def initial_factorized_oriented(key):
    spins,Ks=key
    return (spins,tuple(oriented_intertwiner(v,local_spins(spins,v),Ks[v]).copy() for v in VERT),1+0j)

def apply_hit_branch(branch,x,y,i,j,Jmax2):
    spins,tensors,amp=branch; ei=EIDX[tuple(sorted((x,y)))]; s=spins[ei]; out=[]
    for so in (s-1,s+1):
        if so<0 or so>Jmax2: continue
        Ml,Mr,norm=hit_mats(s,so,x,y,i,j); u,v=sorted((x,y)); tens=list(tensors)
        tens[u]=apply_axis_np(tens[u],LEGIDX[(u,v)],Ml); tens[v]=apply_axis_np(tens[v],LEGIDX[(v,u)],Mr)
        spn=list(spins); spn[ei]=so
        out.append((tuple(spn),tuple(tens),amp*norm))
    return out

def apply_path_branch(branch,path,iout,iin,Jmax2):
    m=len(path)-1; outs=[]
    for mids in itertools.product(range(2),repeat=max(0,m-1)):
        cols=(iout,)+tuple(mids)+(iin,); branches=[branch]
        for r,(x,y) in enumerate(zip(path[:-1],path[1:])):
            nb=[]
            for br in branches: nb.extend(apply_hit_branch(br,x,y,cols[r],cols[r+1],Jmax2))
            branches=nb
            if not branches: break
        outs.extend(branches)
    return outs

@functools.lru_cache(None)
def volume123_matrix(s1,s2,s3):
    mats=[spin_mats_cached(s) for s in (s1,s2,s3)]; d=(s1+1)*(s2+1)*(s3+1); Q=np.zeros((d,d),complex)
    for a,b,c in itertools.product(range(3),repeat=3):
        e=EPS3[a,b,c]
        if e: Q += e*np.kron(np.kron(mats[0][a],mats[1][b]),mats[2][c])
    Q=(Q+Q.conj().T)/2; ev,U=np.linalg.eigh(Q)
    return (U*np.sqrt(np.abs(ev)))@U.conj().T

def apply_volume_tensor(T,spins_local):
    d1,d2,d3,d4=[s+1 for s in spins_local]; V=volume123_matrix(spins_local[0],spins_local[1],spins_local[2])
    return (V@T.reshape(d1*d2*d3,d4)).reshape(d1,d2,d3,d4)

def apply_volume_tensor_oriented(T,spins_local,v):
    X=T
    for leg,w in enumerate(NEIG[v]):
        if w<v: X=apply_axis_np(X,leg,epsilon_j(spins_local[leg]).conj().T)
    X=apply_volume_tensor(X,spins_local)
    for leg,w in enumerate(NEIG[v]):
        if w<v: X=apply_axis_np(X,leg,epsilon_j(spins_local[leg]))
    return X

def add_dict(dst,src,scale=1):
    for k,a in src.items():
        v=dst.get(k,0j)+scale*a
        if abs(v)>1e-13: dst[k]=v
        elif k in dst: del dst[k]

def adjoint_sequence(seq):
    out=[]
    for op in reversed(seq):
        if op[0]=='V': out.append(op)
        else: out.append(('P',tuple(reversed(op[1])),op[3],op[2]))
    return tuple(out)

def T_sequences(v,a,b,c):
    seqs=[]
    for i,j,k in itertools.product(range(2),repeat=3):
        base1=(('V',v),('P',(c,v),k,i),('P',(v,c),j,k)); base2=(('P',(c,v),k,i),('V',v),('P',(v,c),j,k))
        pf=('P',(v,a,b,v),i,j); pr=('P',(v,b,a,v),i,j)
        seqs += [(+1,base1+(pf,)),(-1,base1+(pr,)),(-1,base2+(pf,)),(+1,base2+(pr,))]
    return tuple(seqs)

@functools.lru_cache(None)
def oriented_specs(v):
    neigh=NEIG[v]; specs=[]
    for r in range(4):
        tri=tuple(neigh[i] for i in range(4) if i!=r); sign=(-1)**r; a,b,c=tri
        specs += [(sign,(v,a,b,c)),(sign,(v,b,c,a)),(sign,(v,c,a,b))]
    return tuple(specs)

def apply_primitive_oriented(input_key,seq,Jmax2):
    branches=[initial_factorized_oriented(input_key)]
    for op in seq:
        if op[0]=='V':
            vv=op[1]; nb=[]
            for spins,tensors,amp in branches:
                t=list(tensors); t[vv]=apply_volume_tensor_oriented(t[vv],local_spins(spins,vv),vv); nb.append((spins,tuple(t),amp))
            branches=nb
        else:
            nb=[]
            for br in branches: nb.extend(apply_path_branch(br,op[1],op[2],op[3],Jmax2))
            branches=nb
        if not branches: break
    out={}
    for spins,tensors,amp in branches:
        local_opts=[]; ok=True
        for vv in VERT:
            ls=local_spins(spins,vv); opts=[]
            for K in allowed_k2_t(*ls):
                c=np.vdot(oriented_intertwiner(vv,ls,K),tensors[vv])
                if abs(c)>1e-12: opts.append((K,c))
            if not opts: ok=False; break
            local_opts.append(opts)
        if not ok: continue
        for ch in itertools.product(*local_opts):
            val=amp
            for _,c in ch: val*=c
            if abs(val)>1e-12:
                ko=(spins,tuple(k for k,_ in ch)); out[ko]=out.get(ko,0j)+val
    return {k:v for k,v in out.items() if abs(v)>1e-11}

def apply_T_oriented_basis(input_key,v,a,b,c,Jmax2,adj=False):
    out={}
    for coef,seq in T_sequences(v,a,b,c):
        if adj: seq=adjoint_sequence(seq)
        add_dict(out,apply_primitive_oriented(input_key,seq,Jmax2),coef)
    return out

@functools.lru_cache(None)
def T_cached(key,v,a,b,c,Jmax2,adj): return tuple(apply_T_oriented_basis(key,v,a,b,c,Jmax2,adj).items())

def apply_T_cached_state(state,spec,Jmax2,adj=False):
    out={}
    for key,amp0 in state.items():
        for ko,c in T_cached(key,*spec,Jmax2,adj): out[ko]=out.get(ko,0j)+amp0*c
    return {k:v for k,v in out.items() if abs(v)>1e-10}

def apply_H_cached_state(state,v,Jmax2):
    out={}
    for sign,spec in oriented_specs(v):
        rr=apply_T_cached_state(state,spec,Jmax2,False); aa=apply_T_cached_state(state,spec,Jmax2,True)
        add_dict(out,rr,0.5*sign); add_dict(out,aa,0.5*sign)
    return out

def basis_full_jhalf():
    spins=(1,)*10
    return [(spins,tuple(bits)) for bits in itertools.product((0,2),repeat=5)]

def safe_boundary_master(Jmax2=5):
    basis=basis_full_jhalf(); n=len(basis); M=np.zeros((n,n),complex)
    node_rows=[]; support=set()
    all_cols=[]
    for v in VERT:
        cols=[]
        for key in basis:
            col=apply_H_cached_state({key:1+0j},v,Jmax2)
            cols.append(col); support.update(col.keys())
        G=gram_from_sparse_columns(cols); M += G
        d=kernel_diagnostics(G)
        node_rows.append({'node':v,'rank':d['rank'],'nullity':d['nullity'],'min_positive':d['min_positive'], 'output_union':len(set().union(*(set(c) for c in cols)))})
        all_cols.append(cols)
    D=kernel_diagnostics(M,rtol=1e-9)
    return M,D,node_rows,support,all_cols

def run():
    M,D,rows,support,_=safe_boundary_master(5)
    return {
        'status':'regulator-safe genuine-volume boundary-sector master constraint',
        'domain':'32 all-j=1/2 K5 Gauss-intertwiner boundary states',
        'Jmax':'5/2',
        'master_definition':'M_boundary=sum_v (H_v^safe)^dagger H_v^safe restricted only in the domain, with full safe outputs retained',
        'output_union_size':len(support),
        'node_diagnostics':rows,
        'master_rank':D['rank'],
        'master_nullity':D['nullity'],
        'master_min_eigenvalue':D['min_eigenvalue'],
        'master_min_positive':D['min_positive'],
        'master_eigenvalues':D['eigenvalues'],
        'interpretation':('nullity=0 proves no physical solution lies entirely in the old all-j=1/2 32D boundary sector under the genuine-volume safe Euclidean node constraints; it does not determine the kernel on the full 816/24364 reachable habitat.'),
    }

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path); a=ap.parse_args(); out=run(); text=json.dumps(out,indent=2)
    print(text)
    if a.output: a.output.write_text(text+'\n')

if __name__=='__main__': main()

