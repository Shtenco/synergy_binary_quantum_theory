#!/usr/bin/env python3
"""Reference Peter-Weyl Euclidean operator on an arbitrary tetrahedral PL dual complex.

This is the graph-independent lift of the frozen K5 magnetic Peter-Weyl stack.
All SU(2), CG, edge-hit and absolute-volume primitives are reused unchanged
from k5_peter_weyl_safe_hda_column.py.  Only the graph regulator is generalized:
for local face slots (a,b), the curvature path is the complete dual 2-cell
around the corresponding primal edge, supplied by pl_dual_complex.DualComplex.

On the boundary of a 4-simplex the dual 2-cells are triangles, so the engine
must reproduce the historical K5 amplitudes up to the independently fixed
node-orientation sign.
"""
from __future__ import annotations
import functools,itertools,math,sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW


class PLPeterWeylEuclidean:
    def __init__(self,dual):
        self.dual=dual
        self.VERT=tuple(range(dual.n_tets))
        self.EDGES=tuple(dual.dual_edges())
        self.EIDX={e:i for i,e in enumerate(self.EDGES)}
        self.LEGIDX={}
        for v in self.VERT:
            for r in range(4):
                w=dual.neighbor[(v,r)]
                self.LEGIDX[(v,w)]=r

    def local_spins(self,spins,v):
        return tuple(spins[self.EIDX[tuple(sorted((v,self.dual.neighbor[(v,r)])))]]
                     for r in range(4))

    @functools.lru_cache(maxsize=None)
    def oriented_intertwiner(self,v,spins_local,K):
        T=PW.intertwiner_tensor_cached(tuple(spins_local),K).copy()
        for r in range(4):
            w=self.dual.neighbor[(v,r)]
            if w<v:
                T=PW.apply_axis_np(T,r,PW.epsilon_j(spins_local[r]))
        return T

    def initial_factorized_oriented(self,key):
        spins,Ks=key
        return (spins,tuple(self.oriented_intertwiner(v,self.local_spins(spins,v),Ks[v]).copy()
                           for v in self.VERT),1+0j)

    def apply_hit_branch(self,branch,x,y,i,j,Jmax2):
        spins,tensors,amp=branch
        ei=self.EIDX[tuple(sorted((x,y)))]
        s=spins[ei]; out=[]
        for so in (s-1,s+1):
            if so<0 or so>Jmax2:
                continue
            Ml,Mr,norm=PW.hit_mats(s,so,x,y,i,j)
            u,v=sorted((x,y)); tens=list(tensors)
            tens[u]=PW.apply_axis_np(tens[u],self.LEGIDX[(u,v)],Ml)
            tens[v]=PW.apply_axis_np(tens[v],self.LEGIDX[(v,u)],Mr)
            spn=list(spins); spn[ei]=so
            out.append((tuple(spn),tuple(tens),amp*norm))
        return out

    def apply_path_branch(self,branch,path,iout,iin,Jmax2):
        m=len(path)-1; out=[]
        for mids in itertools.product(range(2),repeat=max(0,m-1)):
            cols=(iout,)+tuple(mids)+(iin,)
            branches=[branch]
            for q,(x,y) in enumerate(zip(path[:-1],path[1:])):
                nxt=[]
                for br in branches:
                    nxt.extend(self.apply_hit_branch(br,x,y,cols[q],cols[q+1],Jmax2))
                branches=nxt
                if not branches:
                    break
            out.extend(branches)
        return out

    def apply_volume_tensor_oriented(self,T,spins_local,v):
        X=T
        for r in range(4):
            w=self.dual.neighbor[(v,r)]
            if w<v:
                X=PW.apply_axis_np(X,r,PW.epsilon_j(spins_local[r]).conj().T)
        X=PW.apply_volume_tensor(X,spins_local)
        for r in range(4):
            w=self.dual.neighbor[(v,r)]
            if w<v:
                X=PW.apply_axis_np(X,r,PW.epsilon_j(spins_local[r]))
        return X

    @staticmethod
    def add(dst,src,scale=1.0,tol=1e-13):
        for k,a in src.items():
            z=dst.get(k,0j)+scale*a
            if abs(z)>tol:
                dst[k]=z
            elif k in dst:
                del dst[k]

    @functools.lru_cache(maxsize=None)
    def primitive_items(self,input_key,seq,Jmax2):
        branches=[self.initial_factorized_oriented(input_key)]
        for op in seq:
            if op[0]=='V':
                v=op[1]; nxt=[]
                for spins,tensors,amp in branches:
                    t=list(tensors)
                    t[v]=self.apply_volume_tensor_oriented(t[v],self.local_spins(spins,v),v)
                    nxt.append((spins,tuple(t),amp))
                branches=nxt
            else:
                nxt=[]
                for br in branches:
                    nxt.extend(self.apply_path_branch(br,op[1],op[2],op[3],Jmax2))
                branches=nxt
            if not branches:
                break
        out={}
        for spins,tensors,amp in branches:
            local_opts=[]; ok=True
            for v in self.VERT:
                ls=self.local_spins(spins,v); opts=[]
                for K in PW.allowed_k2_t(*ls):
                    c=np.vdot(self.oriented_intertwiner(v,ls,K),tensors[v])
                    if abs(c)>1e-12:
                        opts.append((K,c))
                if not opts:
                    ok=False; break
                local_opts.append(opts)
            if not ok:
                continue
            for choice in itertools.product(*local_opts):
                val=amp
                for _,c in choice:
                    val*=c
                if abs(val)>1e-12:
                    key=(spins,tuple(k for k,_ in choice))
                    out[key]=out.get(key,0j)+val
        return tuple((k,a) for k,a in out.items() if abs(a)>1e-11)

    @staticmethod
    def adjoint_sequence(seq):
        out=[]
        for op in reversed(seq):
            if op[0]=='V':
                out.append(op)
            else:
                out.append(('P',tuple(reversed(op[1])),op[3],op[2]))
        return tuple(out)

    @functools.lru_cache(maxsize=None)
    def T_sequences(self,v,ra,rb,rc):
        cnode=self.dual.neighbor[(v,rc)]
        pf=self.dual.plaquette_path(v,ra,rb)
        pr=self.dual.plaquette_path(v,rb,ra)
        out=[]
        for i,j,k in itertools.product(range(2),repeat=3):
            base1=(('V',v),('P',(cnode,v),k,i),('P',(v,cnode),j,k))
            base2=(('P',(cnode,v),k,i),('V',v),('P',(v,cnode),j,k))
            f=('P',pf,i,j); r=('P',pr,i,j)
            out += [(+1,base1+(f,)),(-1,base1+(r,)),
                    (-1,base2+(f,)),(+1,base2+(r,))]
        return tuple(out)

    @functools.lru_cache(maxsize=None)
    def T_items(self,key,v,ra,rb,rc,Jmax2,adj):
        out={}
        for coef,seq0 in self.T_sequences(v,ra,rb,rc):
            seq=self.adjoint_sequence(seq0) if adj else seq0
            self.add(out,dict(self.primitive_items(key,seq,Jmax2)),coef)
        return tuple((k,a) for k,a in out.items() if abs(a)>1e-10)

    def oriented_specs(self,v):
        out=[]
        for omitted in range(4):
            tri=tuple(r for r in range(4) if r!=omitted)
            sign=self.dual.local_sign(v,omitted)
            a,b,c=tri
            out += [(sign,(v,a,b,c)),(sign,(v,b,c,a)),(sign,(v,c,a,b))]
        return tuple(out)

    def H_sine_basis(self,key,v,Jmax2,tol=1e-10):
        out={}
        for sign,spec in self.oriented_specs(v):
            rr=dict(self.T_items(key,*spec,Jmax2,False))
            aa=dict(self.T_items(key,*spec,Jmax2,True))
            self.add(out,rr,-0.5j*sign)
            self.add(out,aa,+0.5j*sign)
        return {k:a for k,a in out.items() if abs(a)>tol}

    def H_sine_state(self,state,v,Jmax2,tol=1e-10):
        out={}
        for key,a0 in state.items():
            self.add(out,self.H_sine_basis(key,v,Jmax2,tol),a0,tol)
        return {k:a for k,a in out.items() if abs(a)>tol}

    @staticmethod
    def norm(state):
        return math.sqrt(sum(abs(a)**2 for a in state.values()))
