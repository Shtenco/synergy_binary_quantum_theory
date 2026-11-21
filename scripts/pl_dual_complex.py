#!/usr/bin/env python3
"""Reusable oriented dual-complex geometry for BCQG PL 3-manifolds.

The microscopic K5 implementation uses triangular regulator loops because the
boundary of a 4-simplex has primal-edge valence three.  On a general
triangulated PL 3-manifold the correct regulator loop is the complete dual
2-cell around the corresponding primal edge.  This module constructs that data
without assuming K5.
"""
from __future__ import annotations
import itertools
from collections import defaultdict,deque


def boundary_4simplex():
    return [tuple(v for v in range(5) if v!=omit) for omit in range(5)]


def seed_16cell_boundary():
    # one sign choice from every antipodal +/-e_i pair
    return sorted(tuple(2*i+b for i,b in enumerate(bits))
                  for bits in itertools.product((0,1),repeat=4))


def faces_by_dim(tets):
    out=defaultdict(set)
    for tet in tets:
        for size in range(1,5):
            for f in itertools.combinations(tet,size):
                out[size-1].add(tuple(sorted(f)))
    return out


def barycentric_subdivision(tets):
    F=faces_by_dim(tets)
    all_faces=sorted(set().union(*F.values()),key=lambda x:(len(x),x))
    fid={f:i for i,f in enumerate(all_faces)}
    out=set()
    for tet in tets:
        for p in itertools.permutations(tet):
            cur=[]; chain=[]
            for v in p:
                cur.append(v); chain.append(fid[tuple(sorted(cur))])
            out.add(tuple(chain))
    return sorted(out)


def permutation_sign(seq):
    inv=sum(seq[i]>seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))
    return -1 if inv%2 else 1


def face_map(tets):
    fm=defaultdict(list)
    for ti,t in enumerate(tets):
        for r in range(4):
            f=tuple(sorted(t[:r]+t[r+1:]))
            fm[f].append(ti)
    return fm


def face_orientation_factor(tet,r):
    ordered=list(tet[:r]+tet[r+1:])
    canon=sorted(ordered); pos={v:i for i,v in enumerate(canon)}
    return ((-1)**r)*permutation_sign([pos[v] for v in ordered])


def orientation_signs(tets):
    """Return +/- orientation multiplier for every tetrahedron.

    Shared triangular faces are required to carry opposite induced
    orientations.  The overall global sign is fixed by sign[0]=+1.
    """
    fm=face_map(tets); adj=[[] for _ in tets]
    for f,ts in fm.items():
        if len(ts)!=2:
            raise ValueError(f'non-closed/non-manifold triangle {f}: incidence={len(ts)}')
        a,b=ts
        ra=next(r for r in range(4) if tuple(sorted(tets[a][:r]+tets[a][r+1:]))==f)
        rb=next(r for r in range(4) if tuple(sorted(tets[b][:r]+tets[b][r+1:]))==f)
        fa=face_orientation_factor(tets[a],ra)
        fb=face_orientation_factor(tets[b],rb)
        rel=-fa*fb
        adj[a].append((b,rel)); adj[b].append((a,rel))
    s=[None]*len(tets); s[0]=1; q=deque([0])
    while q:
        a=q.popleft()
        for b,rel in adj[a]:
            want=s[a]*rel
            if s[b] is None:
                s[b]=want; q.append(b)
            elif s[b]!=want:
                raise ValueError('triangulation is not consistently orientable under the supplied tetrahedron orderings')
    if any(x is None for x in s):
        raise ValueError('disconnected tetrahedral complex')
    return tuple(int(x) for x in s)


class DualComplex:
    def __init__(self,tets):
        self.tets=tuple(tuple(t) for t in tets)
        self.n_tets=len(self.tets)
        self.face_incidence=face_map(self.tets)
        if any(len(ts)!=2 for ts in self.face_incidence.values()):
            raise ValueError('closed PL dual adapter requires every triangle to have incidence two')
        self.orientation=orientation_signs(self.tets)
        self.neighbor={}
        for ti,t in enumerate(self.tets):
            for r in range(4):
                f=tuple(sorted(t[:r]+t[r+1:]))
                a,b=self.face_incidence[f]
                self.neighbor[(ti,r)]=b if a==ti else a
        em=defaultdict(list)
        for ti,t in enumerate(self.tets):
            for e in itertools.combinations(t,2):
                em[tuple(sorted(e))].append(ti)
        self.edge_incidence={e:tuple(ts) for e,ts in em.items()}

    def local_sign(self,tet_id,omitted_face_index):
        """Oriented local epsilon sign for a face-frame slot."""
        return self.orientation[tet_id]*((-1)**omitted_face_index)

    def primal_edge_for_face_pair(self,tet_id,r,s):
        if r==s:
            raise ValueError('face pair must be distinct')
        t=self.tets[tet_id]
        return tuple(sorted(v for k,v in enumerate(t) if k not in (r,s)))

    def neighbors_around_primal_edge(self,tet_id,edge):
        out=[]; t=self.tets[tet_id]
        E=set(edge)
        for r in range(4):
            f=set(t[:r]+t[r+1:])
            if E.issubset(f):
                out.append(self.neighbor[(tet_id,r)])
        # Exactly two faces of a tetrahedron contain a fixed tetrahedral edge.
        if len(set(out))!=2:
            raise RuntimeError((tet_id,edge,out))
        return tuple(out)

    def plaquette_path(self,tet_id,r,s):
        """Closed dual 2-cell path, starting across face r and returning across s.

        The number of dual edges in the path equals the primal-edge valence.
        Swapping r,s returns the reversed path.
        """
        edge=self.primal_edge_for_face_pair(tet_id,r,s)
        support=set(self.edge_incidence[edge])
        nr=self.neighbor[(tet_id,r)]; ns=self.neighbor[(tet_id,s)]
        if nr not in support or ns not in support:
            raise RuntimeError('local face neighbors missing from primal-edge dual cell')
        path=[tet_id,nr]; prev=tet_id; cur=nr
        while cur!=tet_id:
            cand=list(self.neighbors_around_primal_edge(cur,edge))
            nxt=cand[0] if cand[1]==prev else cand[1]
            prev,cur=cur,nxt
            if cur!=tet_id:
                path.append(cur)
            if len(path)>len(support)+1:
                raise RuntimeError('dual-face walk did not close')
        path.append(tet_id)
        if path[-2]!=ns:
            raise RuntimeError(('wrong return face',tet_id,r,s,edge,path,ns))
        if len(path)-1!=len(support):
            raise RuntimeError(('dual-face length mismatch',edge,path,support))
        return tuple(path)

    def dual_edges(self):
        out=set()
        for v in range(self.n_tets):
            for r in range(4):
                out.add(tuple(sorted((v,self.neighbor[(v,r)]))))
        return tuple(sorted(out))
