#!/usr/bin/env python3
from __future__ import annotations
import itertools, numpy as np
import bqg_safe_boundary_master_gate as R
import bqg_reduced_transfer_hamiltonian as H
import bqg_covariant_volume_recoupling as Rec

TOL=1e-12

def add_vec(dst,key,val,tol=TOL):
    if key in dst: z=dst[key]+val
    else: z=np.array(val,complex,copy=True)
    if np.linalg.norm(z)>tol: dst[key]=z
    elif key in dst: del dst[key]

def gauss_basis_labels(spins):
    opts=[R.allowed_k2_t(*R.local_spins(spins,v)) for v in R.VERT]
    return list(itertools.product(*opts))

def initial_reduced_block(spins,B,tol=TOL):
    spins=tuple(spins); B=np.asarray(B,complex)
    ks=gauss_basis_labels(spins)
    assert B.shape[0]==len(ks)
    out={}
    for i,Ks in enumerate(ks):
        row=B[i]
        if np.linalg.norm(row)<=tol: continue
        locs=tuple((K,K,0,0) for K in Ks)
        out[(spins,locs)]=row.copy()
    return out

def apply_hit_block(state,x,y,i,j,Jmax2,tol=TOL):
    out={}; ei=R.EIDX[tuple(sorted((x,y)))]; u,w=sorted((x,y))
    for (spins,locs),amp in state.items():
        s=spins[ei]
        for so in (s-1,s+1):
            if so<0 or so>Jmax2: continue
            if x<y:
                kind_u,col_u='cg',i; kind_w,col_w='cg',j
                norm=np.sqrt((s+1)/(so+1))
            else:
                kind_u,col_u='adj',j; kind_w,col_w='adj',i
                norm=np.sqrt((so+1)/(s+1))
            lsu=R.local_spins(spins,u); lsw=R.local_spins(spins,w)
            legu=R.LEGIDX[(u,w)]; legw=R.LEGIDX[(w,u)]
            ou=H.local_hit_outputs(u,lsu,legu,so,kind_u,col_u,locs[u],tol)
            ow=H.local_hit_outputs(w,lsw,legw,so,kind_w,col_w,locs[w],tol)
            if not ou or not ow: continue
            spn=list(spins); spn[ei]=so; spn=tuple(spn)
            for lu,au in ou:
                for lw,aw in ow:
                    ll=list(locs); ll[u]=lu; ll[w]=lw
                    add_vec(out,(spn,tuple(ll)),amp*(norm*au*aw),tol)
    return out

def apply_path_block(state,path,iout,iin,Jmax2,tol=TOL):
    m=len(path)-1; out={}
    for mids in itertools.product(range(2),repeat=max(0,m-1)):
        cols=(iout,)+tuple(mids)+(iin,); cur=state
        for r,(x,y) in enumerate(zip(path[:-1],path[1:])):
            cur=apply_hit_block(cur,x,y,cols[r],cols[r+1],Jmax2,tol)
            if not cur: break
        for k,a in cur.items(): add_vec(out,k,a,tol)
    return out

def apply_volume_block(state,v,tol=TOL):
    out={}
    for (spins,locs),amp in state.items():
        K,L,J,M=locs[v]; ls=R.local_spins(spins,v)
        labels=Rec.multiplicity_labels(tuple(ls),J)
        try: idx=labels.index((K,L))
        except ValueError: continue
        V=Rec.reduced_volume_block(tuple(ls),J)
        for a,(Ko,Lo) in enumerate(labels):
            z=V[a,idx]
            if abs(z)>tol:
                ll=list(locs); ll[v]=(Ko,Lo,J,M)
                add_vec(out,(spins,tuple(ll)),amp*z,tol)
    return out

def apply_sequence_block(state,seq,Jmax2,tol=TOL):
    cur=state
    for op in seq:
        if op[0]=='V': cur=apply_volume_block(cur,op[1],tol)
        else: cur=apply_path_block(cur,op[1],op[2],op[3],Jmax2,tol)
        if not cur: break
    return cur

def project_gauss_block(red,tol=1e-10):
    # dict output spin -> dense gauss-row x ncol matrix
    tmp={}
    for (spins,locs),vec in red.items():
        Ks=[]; ok=True
        for K,L,J,M in locs:
            if J!=0 or M!=0 or K!=L: ok=False; break
            Ks.append(K)
        if not ok: continue
        so=tuple(spins); key=tuple(Ks)
        tmp.setdefault(so,{})[key]=tmp.setdefault(so,{}).get(key,0)+vec
    out={}
    for so,d in tmp.items():
        labels=gauss_basis_labels(so); idx={k:i for i,k in enumerate(labels)}
        ncol=len(next(iter(d.values()))) if d else 0
        M=np.zeros((len(labels),ncol),complex)
        for k,v in d.items(): M[idx[k]]+=v
        M[np.abs(M)<tol]=0
        if np.linalg.norm(M)>tol: out[so]=M
    return out

def apply_scaled_merge(dst,src,scale=1.0,tol=1e-12):
    for k,a in src.items(): add_vec(dst,k,scale*a,tol)

def radial_forward_block(init,v,c,i,j,Jmax2,tol=1e-12):
    # sum_k P(v,c)_{j,k} P(c,v)_{k,i} V - P(v,c)_{j,k} V P(c,v)_{k,i}
    out={}
    for k in range(2):
        a=apply_volume_block(init,v,tol)
        a=apply_path_block(a,(c,v),k,i,Jmax2,tol)
        a=apply_path_block(a,(v,c),j,k,Jmax2,tol)
        apply_scaled_merge(out,a,+1,tol)
        b=apply_path_block(init,(c,v),k,i,Jmax2,tol)
        b=apply_volume_block(b,v,tol)
        b=apply_path_block(b,(v,c),j,k,Jmax2,tol)
        apply_scaled_merge(out,b,-1,tol)
    return out

def triangle_forward_block(state,v,a,b,i,j,Jmax2,tol=1e-12):
    out={}
    pf=apply_path_block(state,(v,a,b,v),i,j,Jmax2,tol)
    pr=apply_path_block(state,(v,b,a,v),i,j,Jmax2,tol)
    apply_scaled_merge(out,pf,+1,tol); apply_scaled_merge(out,pr,-1,tol)
    return out

def triangle_adjoint_block(state,v,a,b,i,j,Jmax2,tol=1e-12):
    # pf^dag-pr^dag = P(v,b,a,v)_{j,i} - P(v,a,b,v)_{j,i}
    out={}
    pfd=apply_path_block(state,(v,b,a,v),j,i,Jmax2,tol)
    prd=apply_path_block(state,(v,a,b,v),j,i,Jmax2,tol)
    apply_scaled_merge(out,pfd,+1,tol); apply_scaled_merge(out,prd,-1,tol)
    return out

def radial_adjoint_block(state,v,c,i,j,Jmax2,tol=1e-12):
    # adjoint of radial_forward(i,j): sum_k [V P(v,c)_{i,k} P(c,v)_{k,j} - P(v,c)_{i,k} V P(c,v)_{k,j}]
    # In application order: P(c,v)_{k,j}, P(v,c)_{i,k}, V  minus P(c,v)_{k,j}, V, P(v,c)_{i,k}
    out={}
    for k in range(2):
        a=apply_path_block(state,(c,v),k,j,Jmax2,tol)
        a=apply_path_block(a,(v,c),i,k,Jmax2,tol)
        a=apply_volume_block(a,v,tol)
        apply_scaled_merge(out,a,+1,tol)
        b=apply_path_block(state,(c,v),k,j,Jmax2,tol)
        b=apply_volume_block(b,v,tol)
        b=apply_path_block(b,(v,c),i,k,Jmax2,tol)
        apply_scaled_merge(out,b,-1,tol)
    return out

def oriented_T_factorized(init,v,a,b,c,Jmax2,tol=1e-12):
    out={}
    for i in range(2):
        for j in range(2):
            r=radial_forward_block(init,v,c,i,j,Jmax2,tol)
            t=triangle_forward_block(r,v,a,b,i,j,Jmax2,tol)
            apply_scaled_merge(out,t,+1,tol)
            # T^dag = (radial)^dag (triangle)^dag; apply loop adjoint first then radial adjoint
            ld=triangle_adjoint_block(init,v,a,b,i,j,Jmax2,tol)
            rd=radial_adjoint_block(ld,v,c,i,j,Jmax2,tol)
            apply_scaled_merge(out,rd,+1,tol)
    # Hermitian completion 1/2(T+Tdag)
    return {k:0.5*a for k,a in out.items() if np.linalg.norm(a)>tol}

def H_node_apply(spins,v,B,tol=1e-10,Jmax2=None):
    spins=tuple(spins); B=np.asarray(B,complex)
    if Jmax2 is None: Jmax2=max(spins)+1
    init=initial_reduced_block(spins,B)
    acc={}
    for sign,spec in R.oriented_specs(v):
        vv,a,b,c=spec
        rr=oriented_T_factorized(init,vv,a,b,c,Jmax2,1e-12)
        apply_scaled_merge(acc,rr,sign,1e-11)
    return project_gauss_block(acc,tol)

if __name__=='__main__':
    sp=(1,)*10; B=np.eye(1,dtype=complex)
    # one basis row only: embed first gauss basis vector
    BB=np.zeros((32,1),complex);BB[0,0]=1
    o=H_node_apply(sp,0,BB)
    print({'spin_blocks':len(o),'out_dim':sum(M.shape[0] for M in o.values()),'nnz':sum(np.count_nonzero(M) for M in o.values())})

# ---- exact radial-collapsed backend ----
import functools

@functools.lru_cache(None)
def all_local_labels(spins_local):
    sp=tuple(spins_local);out=[]
    for J in Rec.allowed_total_J(sp):
        for M in R.m2vals_t(J):
            for K,L in Rec.multiplicity_labels(sp,J):out.append((K,L,J,M))
    return tuple(out)

@functools.lru_cache(None)
def far_radial_weight(s,so,v,c):
    # Exact contraction of the remote endpoint c, including both link normalization factors.
    W=None
    for k in range(2):
        # first c->v component (k,0); second v->c component (0,k)
        Ml1,Mr1,n1=R.hit_mats(s,so,c,v,k,0)
        u1,w1=sorted((c,v)); A = Ml1 if c==u1 else Mr1
        Ml2,Mr2,n2=R.hit_mats(so,s,v,c,0,k)
        u2,w2=sorted((v,c)); B = Ml2 if c==u2 else Mr2
        Z=(n1*n2)*(B@A); W=Z if W is None else W+Z
    scalar=np.trace(W)/(s+1)
    err=np.linalg.norm(W-scalar*np.eye(s+1))
    if err>1e-10: raise RuntimeError(('far radial not scalar',s,so,v,c,err))
    return complex(scalar)

def _volume_local_outputs(spins_local,label,tol=1e-12):
    K,L,J,M=label; labels=Rec.multiplicity_labels(tuple(spins_local),J)
    try: idx=labels.index((K,L))
    except ValueError:return ()
    V=Rec.reduced_volume_block(tuple(spins_local),J);out=[]
    for a,(Ko,Lo) in enumerate(labels):
        z=V[a,idx]
        if abs(z)>tol:out.append(((Ko,Lo,J,M),complex(z)))
    return tuple(out)

@functools.lru_cache(None)
def radial_local_matrix(v,c,spins_local,i,j,tol=1e-13):
    sp=tuple(spins_local); leg=R.LEGIDX[(v,c)]; s=sp[leg]; labels=all_local_labels(sp); idx={x:n for n,x in enumerate(labels)}
    M=np.zeros((len(labels),len(labels)),complex)
    kind1='cg' if c<v else 'adj'  # c->v, central endpoint v
    kind2='cg' if v<c else 'adj'  # v->c, central endpoint v
    for col,lab in enumerate(labels):
        # V then pair
        for lv,cv in _volume_local_outputs(sp,lab,tol):
            for so in (s-1,s+1):
                if so<0:continue
                sp1=list(sp);sp1[leg]=so;sp1=tuple(sp1);w=far_radial_weight(s,so,v,c)
                for l1,a1 in H.local_hit_outputs(v,sp,leg,so,kind1,i,lv,tol):
                    for l2,a2 in H.local_hit_outputs(v,sp1,leg,s,kind2,j,l1,tol):
                        M[idx[l2],col]+=w*cv*a1*a2
        # pair with V in the middle
        for so in (s-1,s+1):
            if so<0:continue
            sp1=list(sp);sp1[leg]=so;sp1=tuple(sp1);w=far_radial_weight(s,so,v,c)
            for l1,a1 in H.local_hit_outputs(v,sp,leg,so,kind1,i,lab,tol):
                for lv,cv in _volume_local_outputs(sp1,l1,tol):
                    for l2,a2 in H.local_hit_outputs(v,sp1,leg,s,kind2,j,lv,tol):
                        M[idx[l2],col]-=w*a1*cv*a2
    M[np.abs(M)<tol]=0
    return M

def apply_local_matrix_block(state,v,matfun,tol=1e-12):
    out={}
    for (spins,locs),amp in state.items():
        ls=R.local_spins(spins,v); labels=all_local_labels(tuple(ls)); li={x:n for n,x in enumerate(labels)}
        lab=locs[v]
        try: col=li[lab]
        except KeyError:continue
        M=matfun(tuple(ls),tol); vec=M[:,col]
        for r,z in enumerate(vec):
            if abs(z)>tol:
                ll=list(locs);ll[v]=labels[r];add_vec(out,(spins,tuple(ll)),amp*z,tol)
    return out

def apply_radial_collapsed_block(state,v,c,i,j,adj=False,tol=1e-12):
    def mf(ls,tol):
        A=radial_local_matrix(v,c,ls,i,j,tol)
        return A.conj().T if adj else A
    return apply_local_matrix_block(state,v,mf,tol)

def oriented_T_collapsed(init,v,a,b,c,Jmax2,tol=1e-12):
    out={}
    for i in range(2):
        for j in range(2):
            r=apply_radial_collapsed_block(init,v,c,i,j,False,tol)
            t=triangle_forward_block(r,v,a,b,i,j,Jmax2,tol)
            apply_scaled_merge(out,t,+0.5,tol)
            ld=triangle_adjoint_block(init,v,a,b,i,j,Jmax2,tol)
            rd=apply_radial_collapsed_block(ld,v,c,i,j,True,tol)
            apply_scaled_merge(out,rd,+0.5,tol)
    return out

def H_node_apply_fast(spins,v,B,tol=1e-10,Jmax2=None):
    spins=tuple(spins);B=np.asarray(B,complex)
    if Jmax2 is None:Jmax2=max(spins)+1
    init=initial_reduced_block(spins,B);acc={}
    for sign,spec in R.oriented_specs(v):
        vv,a,b,c=spec;rr=oriented_T_collapsed(init,vv,a,b,c,Jmax2,1e-12);apply_scaled_merge(acc,rr,sign,1e-11)
    return project_gauss_block(acc,tol)
# Use the exact radial-collapsed backend as the default block action.
H_node_apply = H_node_apply_fast

# ---- triangle-collapsed local transfer backend ----
def endpoint_kind_color(x,y,ci,cj,endpoint):
    if x<y:
        return ('cg', ci if endpoint==x else cj)
    else:
        return ('adj', ci if endpoint==x else cj)

@functools.lru_cache(None)
def local_two_hit_outputs(v,spins_local,leg1,so1,kind1,color1,leg2,so2,kind2,color2,label,tol=1e-13):
    sp=tuple(spins_local); out={}
    for l1,a1 in H.local_hit_outputs(v,sp,leg1,so1,kind1,color1,label,tol):
        sp1=list(sp); sp1[leg1]=so1; sp1=tuple(sp1)
        for l2,a2 in H.local_hit_outputs(v,sp1,leg2,so2,kind2,color2,l1,tol):
            z=out.get(l2,0j)+a1*a2
            if abs(z)>tol:out[l2]=z
            elif l2 in out:del out[l2]
    return tuple(out.items())

@functools.lru_cache(None)
def local_two_hit_shared_outputs(v,spins_local,leg1,so1,kind1,leg2,so2,kind2,label,tol=1e-13):
    out={}
    for m in range(2):
        for l,z in local_two_hit_outputs(v,spins_local,leg1,so1,kind1,m,leg2,so2,kind2,m,label,tol):
            out[l]=out.get(l,0j)+z
    return tuple((l,z) for l,z in out.items() if abs(z)>tol)

def edge_norm(s,so,x,y):
    return R.hit_mats(s,so,x,y,0,0)[2]

def apply_triangle_path_collapsed(state,path,iout,iin,Jmax2,tol=1e-12):
    # path=(v,a,b,v), three distinct edges, exact factorization over the two internal colors.
    v,a,b,v2=path; assert v==v2 and len({v,a,b})==3
    e1=R.EIDX[tuple(sorted((v,a)))]; e2=R.EIDX[tuple(sorted((a,b)))]; e3=R.EIDX[tuple(sorted((b,v)))]
    out={}
    for (spins,locs),amp in state.items():
        s1,s2,s3=spins[e1],spins[e2],spins[e3]
        for o1 in (s1-1,s1+1):
            if o1<0 or o1>Jmax2:continue
            for o2 in (s2-1,s2+1):
                if o2<0 or o2>Jmax2:continue
                for o3 in (s3-1,s3+1):
                    if o3<0 or o3>Jmax2:continue
                    nrm=edge_norm(s1,o1,v,a)*edge_norm(s2,o2,a,b)*edge_norm(s3,o3,b,v)
                    # local v: edge v-a first (color iout), edge b-v last (color iin)
                    lsv=R.local_spins(spins,v); lv1=R.LEGIDX[(v,a)]; lv3=R.LEGIDX[(v,b)]
                    k1,c1=endpoint_kind_color(v,a,iout,0,v); k3,c3=endpoint_kind_color(b,v,0,iin,v)
                    ov=local_two_hit_outputs(v,lsv,lv1,o1,k1,c1,lv3,o3,k3,c3,locs[v],tol)
                    if not ov:continue
                    # local a: shared m1 on v-a and a-b
                    lsa=R.local_spins(spins,a); la1=R.LEGIDX[(a,v)]; la2=R.LEGIDX[(a,b)]
                    ka1,_=endpoint_kind_color(v,a,0,0,a); ka2,_=endpoint_kind_color(a,b,0,0,a)
                    oa=local_two_hit_shared_outputs(a,lsa,la1,o1,ka1,la2,o2,ka2,locs[a],tol)
                    if not oa:continue
                    # local b: shared m2 on a-b and b-v
                    lsb=R.local_spins(spins,b); lb2=R.LEGIDX[(b,a)]; lb3=R.LEGIDX[(b,v)]
                    kb2,_=endpoint_kind_color(a,b,0,0,b); kb3,_=endpoint_kind_color(b,v,0,0,b)
                    ob=local_two_hit_shared_outputs(b,lsb,lb2,o2,kb2,lb3,o3,kb3,locs[b],tol)
                    if not ob:continue
                    spn=list(spins);spn[e1]=o1;spn[e2]=o2;spn[e3]=o3;spn=tuple(spn)
                    for lvo,av in ov:
                        for lao,aa in oa:
                            for lbo,ab in ob:
                                ll=list(locs);ll[v]=lvo;ll[a]=lao;ll[b]=lbo
                                add_vec(out,(spn,tuple(ll)),amp*(nrm*av*aa*ab),tol)
    return out

def triangle_forward_collapsed(state,v,a,b,i,j,Jmax2,tol=1e-12):
    out={}
    pf=apply_triangle_path_collapsed(state,(v,a,b,v),i,j,Jmax2,tol)
    pr=apply_triangle_path_collapsed(state,(v,b,a,v),i,j,Jmax2,tol)
    apply_scaled_merge(out,pf,+1,tol);apply_scaled_merge(out,pr,-1,tol);return out

def triangle_adjoint_collapsed(state,v,a,b,i,j,Jmax2,tol=1e-12):
    # P(vabv)_{ij}^dag = P(vbav)_{ji}
    out={}
    pfd=apply_triangle_path_collapsed(state,(v,b,a,v),j,i,Jmax2,tol)
    prd=apply_triangle_path_collapsed(state,(v,a,b,v),j,i,Jmax2,tol)
    apply_scaled_merge(out,pfd,+1,tol);apply_scaled_merge(out,prd,-1,tol);return out

def oriented_T_ultrafast(init,v,a,b,c,Jmax2,tol=1e-12):
    out={}
    for i in range(2):
        for j in range(2):
            r=apply_radial_collapsed_block(init,v,c,i,j,False,tol)
            t=triangle_forward_collapsed(r,v,a,b,i,j,Jmax2,tol)
            apply_scaled_merge(out,t,+0.5,tol)
            ld=triangle_adjoint_collapsed(init,v,a,b,i,j,Jmax2,tol)
            rd=apply_radial_collapsed_block(ld,v,c,i,j,True,tol)
            apply_scaled_merge(out,rd,+0.5,tol)
    return out

def H_node_apply_ultrafast(spins,v,B,tol=1e-10,Jmax2=None,no_prune=None):
    spins=tuple(spins);B=np.asarray(B,complex)
    if Jmax2 is None:Jmax2=max(spins)+1
    if no_prune is None:no_prune=(tol==0)
    prune_tol=0.0 if no_prune else 1e-12
    merge_tol=0.0 if no_prune else 1e-11
    init=initial_reduced_block(spins,B,prune_tol);acc={}
    for sign,spec in R.oriented_specs(v):
        vv,a,b,c=spec;rr=oriented_T_ultrafast(init,vv,a,b,c,Jmax2,prune_tol);apply_scaled_merge(acc,rr,sign,merge_tol)
    return project_gauss_block(acc,0.0 if no_prune else tol)
H_node_apply = H_node_apply_ultrafast
