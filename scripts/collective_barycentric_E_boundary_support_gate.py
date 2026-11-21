#!/usr/bin/env python3
"""Exact first-E boundary support of a canonical barycentric tetra block.

Work on the full L1 barycentric subdivision of the 16-cell S3 so regulator
plaquettes are never truncated at a block boundary.  For every coarse tetra
block and every local physical Euclidean oriented spec, count the microscopic
fundamental hits crossing the block boundary.  Starting from fine j=1/2 links,
SU(2) coupling then determines the coarse-face irreps that must be retained by
any target-independent collective producer.

This gate uses only operator support/representation theory, not GR target data
and not Euclidean amplitudes.
"""
from __future__ import annotations
import argparse,itertools,json,sys
from collections import Counter,defaultdict
from pathlib import Path
HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
from pl_dual_complex import DualComplex,seed_16cell_boundary,faces_by_dim


def barycentric_with_parent(tets):
    F=faces_by_dim(tets)
    all_faces=sorted(set().union(*F.values()),key=lambda x:(len(x),x))
    fid={f:i for i,f in enumerate(all_faces)}
    fine=[];parent=[]
    for pi,tet in enumerate(tets):
        for p in itertools.permutations(tet):
            cur=[];chain=[]
            for v in p:
                cur.append(v);chain.append(fid[tuple(sorted(cur))])
            fine.append(tuple(chain));parent.append(pi)
    return fine,parent


def spin_after_hits(n):
    # doubled spins reachable from s=1 under n fundamental +/-1 CG hits,
    # with nonnegative intermediate/final spins.  This gate only needs n<=2.
    S={1}
    for _ in range(n):
        S={y for x in S for y in (x-1,x+1) if y>=0}
    return S


def couple_support(spins):
    total={0}
    for s in spins:
        nxt=set()
        for a in total:
            nxt.update(range(abs(a-s),a+s+1,2))
        total=nxt
    return total


def face_support(edge_hits,ordered_edges):
    opts=[sorted(spin_after_hits(edge_hits.get(e,0))) for e in ordered_edges]
    out=set()
    for spins in itertools.product(*opts):
        out|=couple_support(spins)
    return out


def block_summary(D,parent,parent_id):
    inside={v for v,p in enumerate(parent) if p==parent_id}
    boundary=[e for e in D.dual_edges() if (e[0] in inside)^(e[1] in inside)]
    group=defaultdict(list);edge_group={}
    for e in boundary:
        a,b=e;outside=b if a in inside else a
        g=parent[outside];group[g].append(e);edge_group[e]=g
    if sorted(len(es) for es in group.values())!=[6,6,6,6]:
        raise RuntimeError(('bad boundary partition',parent_id,{g:len(es) for g,es in group.items()}))
    for g in group: group[g]=sorted(group[g])

    stats=Counter();group_hits=Counter();union={g:set() for g in group}
    same_face_support=set();split_face_support=set();baseline=None
    for v in sorted(inside):
        for omitted in range(4):
            tri=tuple(r for r in range(4) if r!=omitted)
            for a,b,c in ((tri[0],tri[1],tri[2]),(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])):
                hits=Counter()
                w=D.neighbor[(v,c)];hits[tuple(sorted((v,w)))]+=2
                p=D.plaquette_path(v,a,b)
                for x,y in zip(p[:-1],p[1:]):hits[tuple(sorted((x,y)))]+=1
                bh={e:n for e,n in hits.items() if e in edge_group}
                gs={edge_group[e] for e in bh}
                stats[(sum(bh.values()),len(bh),len(gs))]+=1
                by=Counter()
                for e,n in bh.items():by[edge_group[e]]+=n
                group_hits[tuple(sorted(by.values()))]+=1
                for g in group:
                    S=face_support(bh,group[g]);union[g]|=S
                    if g in by and tuple(sorted(by.values()))==(2,): same_face_support|=S
                    if g in by and tuple(sorted(by.values()))==(1,1): split_face_support|=S
                if not bh and baseline is None:
                    baseline=face_support({},next(iter(group.values())))
    return {
      'parent':parent_id,'fine_tetrahedra':len(inside),'boundary_links':len(boundary),
      'boundary_groups':len(group),'boundary_links_per_face':sorted(len(es) for es in group.values()),
      'pattern_stats':{str(k):v for k,v in sorted(stats.items())},
      'coarse_face_group_hit_totals':{str(k):v for k,v in sorted(group_hits.items())},
      'baseline_face_total_spin2_support':sorted(baseline),
      'same_face_two_hit_spin2_support':sorted(same_face_support),
      'split_faces_one_hit_each_spin2_support':sorted(split_face_support),
      'union_face_total_spin2_support':sorted(set().union(*union.values()))
    }


def run():
    coarse=seed_16cell_boundary();fine,parent=barycentric_with_parent(coarse);D=DualComplex(fine)
    rows=[block_summary(D,parent,p) for p in range(len(coarse))]
    canonical=rows[0]
    same=all({k:v for k,v in r.items() if k!='parent'}=={k:v for k,v in canonical.items() if k!='parent'} for r in rows)
    expect={
      'fine_tetrahedra':24,'boundary_links':24,'boundary_groups':4,
      'boundary_links_per_face':[6,6,6,6],
      'baseline_face_total_spin2_support':[0,2,4,6],
      'same_face_two_hit_spin2_support':[0,2,4,6,8],
      'split_faces_one_hit_each_spin2_support':[1,3,5,7],
      'union_face_total_spin2_support':list(range(9))
    }
    pattern_expected={str((0,0,0)):72,str((2,1,1)):72,str((2,2,1)):96,str((2,2,2)):48}
    group_expected={str(()):72,str((1,1)):48,str((2,)):168}
    checks={k:canonical[k]==v for k,v in expect.items()}
    checks['all_16_blocks_identical']=same
    checks['288_specs']=sum(canonical['pattern_stats'].values())==288
    checks['pattern_distribution']=canonical['pattern_stats']==pattern_expected
    checks['group_hit_distribution']=canonical['coarse_face_group_hit_totals']==group_expected
    return {
      'status':'target-independent one-E support of canonical barycentric block',
      'passed':all(checks.values()),'checks':checks,'canonical_block':canonical,
      'all_blocks_identical':same,
      'required_coarse_face_spins':[x/2 for x in range(9)],
      'interpretation':'A static maximal j=3 face is not closed under one production Euclidean action. The direct collective block basis must retain coarse-face irreps j=0,1/2,...,4 (with their multiplicity channels) before target comparison.',
      'amplitude_note':'This is exact support geometry/CG reachability; zero amplitudes from later dynamical cancellations may reduce realized support, but no sector may be removed before that amplitude calculation.',
      'scope_note':'One-E support only. Hermitian Lorentzian S and repeated collective actions require additional support closure gates.'
    }

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path)
    a=ap.parse_args();o=run();t=json.dumps(o,indent=2);print(t)
    if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(t+'\n',encoding='utf-8')
    return 0 if o['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
