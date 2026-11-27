#!/usr/bin/env python3
from __future__ import annotations
import functools, itertools, math
import numpy as np
import bqg_safe_boundary_master_gate as R
import bqg_covariant_volume_recoupling as Rec

# Reduced local basis label = (K,L,J,M) in doubled-spin convention.

def add_state(dst,src,scale=1.0,tol=1e-13):
    for k,a in src.items():
        z=dst.get(k,0j)+scale*a
        if abs(z)>tol: dst[k]=z
        elif k in dst: del dst[k]

@functools.lru_cache(None)
def oriented_coupled_tensor(v,spins_local,K,L,J,M):
    T=Rec.coupled_tensor(tuple(spins_local),K,L,J,M).copy()
    for leg,w in enumerate(R.NEIG[v]):
        if w<v:
            T=R.apply_axis_np(T,leg,R.epsilon_j(spins_local[leg]))
    return T

@functools.lru_cache(None)
def local_hit_outputs(v,spins_local,leg,s_out,kind,color,label):
    """Exact one-leg reduced transfer in the oriented coupled basis.

    kind='cg'  : A = <s_in,1/2|s_out> component color
    kind='adj' : A = corresponding reversed-edge adjoint component
    Returns tuple((K,L,J,M), amplitude).
    """
    spins_local=tuple(spins_local); K,L,J,M=label; s_in=spins_local[leg]
    a2=1 if color==0 else -1
    if kind=='cg':
        A=R.cg_leg_matrix(s_in,s_out,a2)
    elif kind=='adj':
        A=R.cg_leg_matrix(s_out,s_in,a2).conj().T
    else:
        raise ValueError(kind)
    spo=list(spins_local); spo[leg]=s_out; spo=tuple(spo)
    Bin=oriented_coupled_tensor(v,spins_local,K,L,J,M)
    Y=R.apply_axis_np(Bin,leg,A)
    outs=[]
    # Rank-1/2 operator allows only J' = J +/- 1 (doubled notation),
    # but we deliberately retain an exact overlap filter for convention safety.
    for Jo in Rec.allowed_total_J(spo):
        if abs(Jo-J)!=1: continue
        for Mo in R.m2vals_t(Jo):
            for Ko,Lo in Rec.multiplicity_labels(spo,Jo):
                Bout=oriented_coupled_tensor(v,spo,Ko,Lo,Jo,Mo)
                a=np.vdot(Bout,Y)
                if abs(a)>1e-12: outs.append(((Ko,Lo,Jo,Mo),complex(a)))
    return tuple(outs)


def gauss_key_to_reduced(key):
    spins,Ks=key
    return (tuple(spins),tuple((K,K,0,0) for K in Ks))


def reduced_key_to_gauss(key):
    spins,locs=key; Ks=[]
    for K,L,J,M in locs:
        if J!=0 or M!=0 or K!=L: return None
        Ks.append(K)
    return (tuple(spins),tuple(Ks))


def apply_hit_reduced(state,x,y,i,j,Jmax2):
    out={}; ei=R.EIDX[tuple(sorted((x,y)))]
    u,w=sorted((x,y))
    for (spins,locs),amp in state.items():
        s=spins[ei]
        for so in (s-1,s+1):
            if so<0 or so>Jmax2: continue
            if x<y:
                kind_u,col_u='cg',i; kind_w,col_w='cg',j
                norm=math.sqrt((s+1)/(so+1))
            else:
                kind_u,col_u='adj',j; kind_w,col_w='adj',i
                norm=math.sqrt((so+1)/(s+1))
            lsu=R.local_spins(spins,u); lsw=R.local_spins(spins,w)
            legu=R.LEGIDX[(u,w)]; legw=R.LEGIDX[(w,u)]
            ou=local_hit_outputs(u,lsu,legu,so,kind_u,col_u,locs[u])
            ow=local_hit_outputs(w,lsw,legw,so,kind_w,col_w,locs[w])
            if not ou or not ow: continue
            spn=list(spins); spn[ei]=so; spn=tuple(spn)
            for lu,au in ou:
                for lw,aw in ow:
                    ll=list(locs); ll[u]=lu; ll[w]=lw
                    key=(spn,tuple(ll)); z=amp*norm*au*aw
                    out[key]=out.get(key,0j)+z
    return {k:v for k,v in out.items() if abs(v)>1e-12}


def apply_path_reduced(state,path,iout,iin,Jmax2):
    m=len(path)-1; out={}
    for mids in itertools.product(range(2),repeat=max(0,m-1)):
        cols=(iout,)+tuple(mids)+(iin,); cur=state
        for r,(x,y) in enumerate(zip(path[:-1],path[1:])):
            cur=apply_hit_reduced(cur,x,y,cols[r],cols[r+1],Jmax2)
            if not cur: break
        add_state(out,cur)
    return out


def apply_volume_reduced(state,v):
    out={}
    for (spins,locs),amp in state.items():
        K,L,J,M=locs[v]; ls=R.local_spins(spins,v)
        labels=Rec.multiplicity_labels(tuple(ls),J)
        try: idx=labels.index((K,L))
        except ValueError: continue
        V=Rec.reduced_volume_block(tuple(ls),J)
        for a,(Ko,Lo) in enumerate(labels):
            z=V[a,idx]
            if abs(z)>1e-12:
                ll=list(locs); ll[v]=(Ko,Lo,J,M)
                key=(spins,tuple(ll)); out[key]=out.get(key,0j)+amp*z
    return out


def apply_sequence_reduced(state,seq,Jmax2):
    cur=state
    for op in seq:
        if op[0]=='V': cur=apply_volume_reduced(cur,op[1])
        else: cur=apply_path_reduced(cur,op[1],op[2],op[3],Jmax2)
        if not cur: break
    return cur


def apply_primitive_gauss(input_state,seq,Jmax2):
    # input_state uses ordinary Gauss keys.
    red={gauss_key_to_reduced(k):a for k,a in input_state.items()}
    rr=apply_sequence_reduced(red,seq,Jmax2); out={}
    for k,a in rr.items():
        g=reduced_key_to_gauss(k)
        if g is not None and abs(a)>1e-11: out[g]=out.get(g,0j)+a
    return out

@functools.lru_cache(None)
def primitive_on_key(key,seq,Jmax2):
    return tuple(apply_primitive_gauss({key:1+0j},seq,Jmax2).items())


def apply_H_gauss_state(state,v,Jmax2):
    out={}
    for sign,spec in R.oriented_specs(v):
        vv,a,b,c=spec
        for coef,seq0 in R.T_sequences(vv,a,b,c):
            for adj in (False,True):
                seq=R.adjoint_sequence(seq0) if adj else seq0
                overall=0.5*sign*coef
                for key,amp in state.items():
                    for ko,z in primitive_on_key(key,seq,Jmax2):
                        out[ko]=out.get(ko,0j)+overall*amp*z
    return {k:v for k,v in out.items() if abs(v)>1e-10}

if __name__=='__main__':
    key=R.basis_full_jhalf()[0]
    o=apply_H_gauss_state({key:1+0j},0,3)
    print({'output_keys':len(o),'norm2':sum(abs(x)**2 for x in o.values())})


def apply_H_gauss_state_collective(state,v,Jmax2):
    """Exact H action, applying each primitive to the whole superposition at once.
    This merges identical reduced intermediate states after every edge hit.
    """
    out={}
    for sign,spec in R.oriented_specs(v):
        vv,a,b,c=spec
        for coef,seq0 in R.T_sequences(vv,a,b,c):
            for adj in (False,True):
                seq=R.adjoint_sequence(seq0) if adj else seq0
                rr=apply_primitive_gauss(state,seq,Jmax2)
                add_state(out,rr,0.5*sign*coef,tol=1e-11)
    return {k:z for k,z in out.items() if abs(z)>1e-10}


def apply_oriented_spec_gauss_state(state,sign,spec,Jmax2):
    """Contribution of one oriented (v,a,b,c) term to the Hermitian H_v."""
    vv,a,b,c=spec; out={}
    for coef,seq0 in R.T_sequences(vv,a,b,c):
        for adj in (False,True):
            seq=R.adjoint_sequence(seq0) if adj else seq0
            rr=apply_primitive_gauss(state,seq,Jmax2)
            add_state(out,rr,0.5*sign*coef,tol=1e-11)
    return {k:z for k,z in out.items() if abs(z)>1e-10}
