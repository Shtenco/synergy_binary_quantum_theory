#!/usr/bin/env python3
from __future__ import annotations
import functools, itertools, json
import numpy as np
import bqg_safe_boundary_master_gate as R


def identity_perm():
    return (0,1,2,3,4)


def compose_perm(q,p):
    """Composition q after p for old->new tuples: (q o p)[i] = q[p[i]]."""
    return tuple(q[p[i]] for i in range(5))


def inverse_perm(p):
    inv=[0]*5
    for i,j in enumerate(p): inv[j]=i
    return tuple(inv)


def permute_spins(spins,p):
    out=[None]*len(R.EDGES)
    for old,(u,v) in enumerate(R.EDGES):
        tgt=R.EIDX[tuple(sorted((p[u],p[v])))]
        out[tgt]=spins[old]
    return tuple(out)


def edge_orientation_phase(spins,p):
    """SU(2) dualization phase from edges whose canonical orientation flips.

    Edges are canonically oriented from the lower to the higher K5 vertex label.
    Under a vertex permutation p, an old edge (u<v) is reversed exactly when
    p[u] > p[v]. Re-expressing the image in the canonical target orientation
    contributes the Frobenius-Schur factor (-1)^(2j)=(-1)^s for that edge.
    """
    exponent=0
    for e,(u,v) in enumerate(R.EDGES):
        if p[u] > p[v]: exponent += int(spins[e])
    return -1 if exponent%2 else 1


@functools.lru_cache(None)
def _local_recoupling_cached(spins,v,p):
    spins=tuple(spins); p=tuple(p)
    vp=p[v]
    target_spins=permute_spins(spins,p)
    old_ls=R.local_spins(spins,v)
    new_ls=R.local_spins(target_spins,vp)
    oldKs=R.allowed_k2_t(*old_ls)
    newKs=R.allowed_k2_t(*new_ls)
    pinv=inverse_perm(p)
    old_neigh=R.NEIG[v]
    new_neigh=R.NEIG[vp]
    axes=tuple(old_neigh.index(pinv[wp]) for wp in new_neigh)
    M=np.zeros((len(newKs),len(oldKs)),dtype=complex)
    for j,K in enumerate(oldKs):
        T=np.transpose(R.intertwiner_tensor_cached(old_ls,K),axes)
        for i,Kp in enumerate(newKs):
            M[i,j]=np.vdot(R.intertwiner_tensor_cached(new_ls,Kp),T)
    return M,oldKs,newKs


def local_recoupling_matrix(spins,v,p):
    M,oldKs,newKs=_local_recoupling_cached(tuple(spins),int(v),tuple(p))
    return M.copy(),tuple(oldKs),tuple(newKs)


def act_on_key(key,p,tol=0.0):
    spins,Ks=key; p=tuple(p)
    tsp=permute_spins(spins,p)
    per_target=[None]*5
    for v in R.VERT:
        M,oldKs,newKs=_local_recoupling_cached(tuple(spins),v,p)
        try: col=oldKs.index(Ks[v])
        except ValueError: return {}
        opts=[]
        for row,Kp in enumerate(newKs):
            c=M[row,col]
            if abs(c)>tol: opts.append((Kp,c))
        if not opts: return {}
        per_target[p[v]]=opts
    out={}
    orient=edge_orientation_phase(spins,p)
    for choices in itertools.product(*per_target):
        newKs=tuple(k for k,_ in choices)
        c=complex(orient)
        for _,a in choices: c*=a
        if abs(c)>tol: out[(tsp,newKs)]=out.get((tsp,newKs),0j)+c
    return {k:v for k,v in out.items() if abs(v)>tol}


def act_on_state(state,p,tol=0.0):
    out={}
    for key,amp in state.items():
        for kp,c in act_on_key(key,p,tol).items():
            z=out.get(kp,0j)+amp*c
            if abs(z)>tol: out[kp]=z
            elif kp in out: del out[kp]
    return out


def state_norm2(state):
    return float(sum(abs(v)**2 for v in state.values()))


def state_distance(a,b):
    keys=set(a)|set(b)
    return float(np.sqrt(sum(abs(a.get(k,0j)-b.get(k,0j))**2 for k in keys)))


def run():
    p=(1,0,2,3,4); q=(0,2,1,3,4)
    key=((1,)*10,(0,2,2,0,2))
    group_err=state_distance(act_on_state(act_on_state({key:1+0j},p),q),act_on_state({key:1+0j},compose_perm(q,p)))
    unit=[]
    for v in R.VERT:
        M,_,_=local_recoupling_matrix((1,)*10,v,p)
        unit.append(float(np.linalg.norm(M.conj().T@M-np.eye(M.shape[1]))))
    return {
        'status':'exact S5 vertex-permutation action on K5 spin-network recoupling basis',
        'group_law_error':group_err,
        'max_local_unitarity_error':max(unit),
        'scope':'basis-recoupling action; covariance with the regulator-safe Hamiltonian is a separate gate',
    }

if __name__=='__main__':
    print(json.dumps(run(),indent=2))
