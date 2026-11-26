#!/usr/bin/env python3
from __future__ import annotations
import itertools,json
V=range(5); EDGES=list(itertools.combinations(V,2)); EI={e:i for i,e in enumerate(EDGES)}

def specs(v):
    neigh=sorted(w for w in V if w!=v);out=[]
    for r in range(4):
        tri=tuple(neigh[i] for i in range(4) if i!=r); a,b,c=tri
        out.extend([(v,a,b,c),(v,b,c,a),(v,c,a,b)])
    return out

def mult4(ss):
    s1,s2,s3,s4=ss
    return len(set(range(abs(s1-s2),s1+s2+1,2)) & set(range(abs(s3-s4),s3+s4+1,2)))

def gauss_mult(spins):
    m=1
    for v in V:
        ss=[spins[EI[tuple(sorted((v,w)))]] for w in V if w!=v]
        m*=mult4(ss)
    return m

def triangle_hit_states(spins,spec):
    v,a,b,c=spec
    tri=[EI[tuple(sorted((v,a)))],EI[tuple(sorted((a,b)))],EI[tuple(sorted((b,v)))]]
    states={tuple(spins)}
    for e in tri:
        nxt=set()
        for st in states:
            s=st[e]
            if s>0:
                q=list(st);q[e]=s-1;nxt.add(tuple(q))
            q=list(st);q[e]=s+1;nxt.add(tuple(q))
        states=nxt
    return states

def apply_H(spinset,v):
    out=set()
    for s in spinset:
        for sp in specs(v):out |= triangle_hit_states(s,sp)
    return out

def summarize(sset):
    valid=[s for s in sset if gauss_mult(s)>0]
    return {'raw_spin_assignments':len(sset),'Gauss_admissible_spin_assignments':len(valid),
            'spin_network_basis_dimension':sum(gauss_mult(s) for s in valid),
            'max_final_spin':max(max(s) for s in valid)/2 if valid else 0.0},set(valid)

def run():
    start={(1,)*10}; one={}; valid1={}
    for v in V:
        sm,vs=summarize(apply_H(start,v));one[str(v)]=sm;valid1[v]=vs
    pairs={}
    for v,w in itertools.combinations(V,2):
        both=apply_H(valid1[v],w)|apply_H(valid1[w],v)
        sm,_=summarize(both);pairs[f'{v}-{w}']=sm
    # Global nested shells under all five effective final-support H moves.
    shells=[];cur=start
    for depth in range(4):
        sm,valid=summarize(cur);sm['depth']=depth;shells.append(sm)
        if depth<3:
            nxt=set()
            for v in V:nxt |= apply_H(valid,v)
            cur=nxt
    # Safe intermediate wall: each H has one radial intermediate +/-1 and
    # final triangle hits +/-1. Starting from s=1, depth d final max <=1+d;
    # an H applied at depth d-1 may need one extra radial intermediate, but
    # that equals the same 1+d wall.
    return {'status':'K5 reachability with exact radial-return selection rule',
            'radial_final_rule':'s_radial,out = s_radial,in exactly; only the three triangle links change by +/-1 in final support',
            'one_H_by_vertex':one,'HH_by_vertex_pair':pairs,'global_shells':shells,
            'safe_wall_claim':'For d successive H actions from all s=1, Jmax_safe=(1+d)/2 suffices for both final triangle support and the one-step radial intermediate. Thus HH requires Jmax=3/2, not 5/2.',
            'scope':'Representation-selection after exact CG radial-return. Amplitude cancellations can reduce support further.'}
if __name__=='__main__':print(json.dumps(run(),indent=2))
