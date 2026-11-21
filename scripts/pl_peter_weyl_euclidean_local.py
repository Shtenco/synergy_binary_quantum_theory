#!/usr/bin/env python3
"""Exact active-cone implementation of the PL Peter-Weyl Euclidean operator.

The reference graph-independent engine stores one intertwiner tensor for every
node in every primitive branch.  That is exact but unnecessarily expensive on
refined complexes: a local Euclidean term touches only the source node, its
source link, and the nodes on one dual plaquette.  Every untouched node remains
in its input Gauss intertwiner and contributes the exact overlap 1.

This class preserves the full global spin-network key but stores tensors only
for nodes actually touched by V or a holonomy hit.  Final Gauss projection is
performed only on that active cone; all untouched K labels are copied exactly.
No support, amplitude, cutoff, sign, ordering, or physical approximation is
changed.
"""
from __future__ import annotations

import functools
import itertools
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0,str(HERE))

import pl_peter_weyl_euclidean as BASE


class LocalPLPeterWeylEuclidean(BASE.PLPeterWeylEuclidean):
    """Drop-in exact local-cone backend for ``PLPeterWeylEuclidean``."""

    def _spin(self,base_spins,delta,ei):
        return delta.get(ei,base_spins[ei])

    def _local_spins(self,base_spins,delta,v):
        return tuple(
            self._spin(
                base_spins,
                delta,
                self.EIDX[tuple(sorted((v,self.dual.neighbor[(v,r)])))],
            )
            for r in range(4)
        )

    @functools.lru_cache(maxsize=None)
    def primitive_items(self,input_key,seq,Jmax2):
        base_spins,base_Ks=input_key

        # branch = (changed_spin_dict, active_tensor_dict, amplitude)
        branches=[({}, {}, 1+0j)]

        def ensure_tensor(tensors,delta,v):
            if v in tensors:
                return tensors
            ls=self._local_spins(base_spins,delta,v)
            K=base_Ks[v]
            if K not in BASE.PW.allowed_k2_t(*ls):
                raise RuntimeError(('untouched input K incompatible with current local spins',v,ls,K))
            out=tensors.copy()
            out[v]=self.oriented_intertwiner(v,ls,K).copy()
            return out

        for op in seq:
            if op[0]=='V':
                v=op[1]
                nxt=[]
                for delta,tensors,amp in branches:
                    tensors=ensure_tensor(tensors,delta,v)
                    t=tensors.copy()
                    ls=self._local_spins(base_spins,delta,v)
                    t[v]=self.apply_volume_tensor_oriented(t[v],ls,v)
                    nxt.append((delta,t,amp))
                branches=nxt
            else:
                path,iout,iin=op[1],op[2],op[3]
                m=len(path)-1
                out=[]
                for mids in itertools.product(range(2),repeat=max(0,m-1)):
                    cols=(iout,)+tuple(mids)+(iin,)
                    work=branches
                    for q,(x,y) in enumerate(zip(path[:-1],path[1:])):
                        hit=[]
                        for delta,tensors,amp in work:
                            ei=self.EIDX[tuple(sorted((x,y)))]
                            s=self._spin(base_spins,delta,ei)
                            tensors=ensure_tensor(tensors,delta,x)
                            tensors=ensure_tensor(tensors,delta,y)
                            for so in (s-1,s+1):
                                if so<0 or so>Jmax2:
                                    continue
                                Ml,Mr,norm=BASE.PW.hit_mats(s,so,x,y,cols[q],cols[q+1])
                                u,v=sorted((x,y))
                                t=tensors.copy()
                                t[u]=BASE.PW.apply_axis_np(t[u],self.LEGIDX[(u,v)],Ml)
                                t[v]=BASE.PW.apply_axis_np(t[v],self.LEGIDX[(v,u)],Mr)
                                d=delta.copy()
                                if so==base_spins[ei]:
                                    d.pop(ei,None)
                                else:
                                    d[ei]=so
                                hit.append((d,t,amp*norm))
                        work=hit
                        if not work:
                            break
                    out.extend(work)
                branches=out
            if not branches:
                break

        result={}
        for delta,tensors,amp in branches:
            if abs(amp)<=1e-14:
                continue
            spins=list(base_spins)
            for ei,s in delta.items():
                spins[ei]=s
            spins=tuple(spins)

            active=sorted(tensors)
            local_opts=[]
            ok=True
            for v in active:
                ls=self.local_spins(spins,v)
                opts=[]
                for K in BASE.PW.allowed_k2_t(*ls):
                    c=BASE.np.vdot(self.oriented_intertwiner(v,ls,K),tensors[v])
                    if abs(c)>1e-12:
                        opts.append((K,c))
                if not opts:
                    ok=False
                    break
                local_opts.append(opts)
            if not ok:
                continue

            for choice in itertools.product(*local_opts):
                val=amp
                Ks=list(base_Ks)
                for v,(K,c) in zip(active,choice):
                    Ks[v]=K
                    val*=c
                if abs(val)>1e-12:
                    key=(spins,tuple(Ks))
                    result[key]=result.get(key,0j)+val

        return tuple((k,a) for k,a in result.items() if abs(a)>1e-11)
