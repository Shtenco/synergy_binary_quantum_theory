#!/usr/bin/env python3
from __future__ import annotations
import numpy as np
import bqg_block_apply_invariants as A
import bqg_safe_boundary_master_gate as R


def compatible_triangle_target(spins,target,v,a,b,Jmax2):
    spins=tuple(spins); target=tuple(target)
    es=(R.EIDX[tuple(sorted((v,a)))],R.EIDX[tuple(sorted((a,b)))],R.EIDX[tuple(sorted((b,v)))])
    eset=set(es)
    for e,(x,y) in enumerate(zip(spins,target)):
        if e in eset:
            if y not in (x-1,x+1) or y<0 or y>Jmax2:return False
        elif y!=x:return False
    return True


def apply_triangle_path_collapsed_targets(state,path,iout,iin,Jmax2,targets,tol=1e-12):
    v,a,b,v2=path; assert v==v2 and len({v,a,b})==3
    e1=R.EIDX[tuple(sorted((v,a)))]; e2=R.EIDX[tuple(sorted((a,b)))]; e3=R.EIDX[tuple(sorted((b,v)))]
    out={}
    for (spins,locs),amp in state.items():
        s1,s2,s3=spins[e1],spins[e2],spins[e3]
        for target in targets:
            if not compatible_triangle_target(spins,target,v,a,b,Jmax2):continue
            o1,o2,o3=target[e1],target[e2],target[e3]
            nrm=A.edge_norm(s1,o1,v,a)*A.edge_norm(s2,o2,a,b)*A.edge_norm(s3,o3,b,v)
            lsv=R.local_spins(spins,v); lv1=R.LEGIDX[(v,a)]; lv3=R.LEGIDX[(v,b)]
            k1,c1=A.endpoint_kind_color(v,a,iout,0,v); k3,c3=A.endpoint_kind_color(b,v,0,iin,v)
            ov=A.local_two_hit_outputs(v,lsv,lv1,o1,k1,c1,lv3,o3,k3,c3,locs[v])
            if not ov:continue
            lsa=R.local_spins(spins,a); la1=R.LEGIDX[(a,v)]; la2=R.LEGIDX[(a,b)]
            ka1,_=A.endpoint_kind_color(v,a,0,0,a); ka2,_=A.endpoint_kind_color(a,b,0,0,a)
            oa=A.local_two_hit_shared_outputs(a,lsa,la1,o1,ka1,la2,o2,ka2,locs[a])
            if not oa:continue
            lsb=R.local_spins(spins,b); lb2=R.LEGIDX[(b,a)]; lb3=R.LEGIDX[(b,v)]
            kb2,_=A.endpoint_kind_color(a,b,0,0,b); kb3,_=A.endpoint_kind_color(b,v,0,0,b)
            ob=A.local_two_hit_shared_outputs(b,lsb,lb2,o2,kb2,lb3,o3,kb3,locs[b])
            if not ob:continue
            spn=tuple(target)
            for lvo,av in ov:
                for lao,aa in oa:
                    for lbo,ab in ob:
                        ll=list(locs);ll[v]=lvo;ll[a]=lao;ll[b]=lbo
                        A.add_vec(out,(spn,tuple(ll)),amp*(nrm*av*aa*ab),tol)
    return out


def triangle_forward_targets(state,v,a,b,i,j,Jmax2,targets,tol=1e-12):
    out={}
    pf=apply_triangle_path_collapsed_targets(state,(v,a,b,v),i,j,Jmax2,targets,tol)
    pr=apply_triangle_path_collapsed_targets(state,(v,b,a,v),i,j,Jmax2,targets,tol)
    A.apply_scaled_merge(out,pf,+1,tol);A.apply_scaled_merge(out,pr,-1,tol);return out

def triangle_adjoint_targets(state,v,a,b,i,j,Jmax2,targets,tol=1e-12):
    out={}
    pfd=apply_triangle_path_collapsed_targets(state,(v,b,a,v),j,i,Jmax2,targets,tol)
    prd=apply_triangle_path_collapsed_targets(state,(v,a,b,v),j,i,Jmax2,targets,tol)
    A.apply_scaled_merge(out,pfd,+1,tol);A.apply_scaled_merge(out,prd,-1,tol);return out


def oriented_T_targets(init,v,a,b,c,Jmax2,targets,tol=1e-12):
    out={}
    for i in range(2):
        for j in range(2):
            r=A.apply_radial_collapsed_block(init,v,c,i,j,False,tol)
            t=triangle_forward_targets(r,v,a,b,i,j,Jmax2,targets,tol)
            A.apply_scaled_merge(out,t,+0.5,tol)
            ld=triangle_adjoint_targets(init,v,a,b,i,j,Jmax2,targets,tol)
            rd=A.apply_radial_collapsed_block(ld,v,c,i,j,True,tol)
            A.apply_scaled_merge(out,rd,+0.5,tol)
    return out


def H_node_apply_targets(spins,v,B,targets,tol=1e-10,Jmax2=None):
    spins=tuple(spins);B=np.asarray(B,complex);targets=tuple(sorted(set(map(tuple,targets))))
    if Jmax2 is None:Jmax2=max(spins)+1
    # filter impossible globally once
    init=A.initial_reduced_block(spins,B);acc={}
    for sign,spec in R.oriented_specs(v):
        vv,a,b,c=spec
        ts=tuple(t for t in targets if compatible_triangle_target(spins,t,vv,a,b,Jmax2))
        if not ts:continue
        rr=oriented_T_targets(init,vv,a,b,c,Jmax2,ts,1e-12)
        A.apply_scaled_merge(acc,rr,sign,1e-11)
    return A.project_gauss_block(acc,tol)

