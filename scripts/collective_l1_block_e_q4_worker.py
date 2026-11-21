#!/usr/bin/env python3
"""Exact q=4 projection of one L1 block Euclidean amplitude column.

For one fine chamber u in parent coarse tetrahedron 0, compute only oriented
Euclidean specs whose dual plaquette has length q=4.  The result is then
projected to final basis states with exactly four changed microscopic doubled
spins relative to the homogeneous all-j=1/2 seed.

This projection is exact, not approximate: every plaquette edge receives one
fundamental hit and therefore changes 2j=1 to 0 or 2.  Hence q=6 and q=8 terms
cannot contribute to the exactly-four-changed-edge sector.  The separate source
link receives two hits and is disjoint from the plaquette; exactly four changed
edges selects the branch where that source link returns to its seed spin.
"""
from __future__ import annotations
import argparse,json,math,sys,traceback
from pathlib import Path

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0,str(HERE))

import peter_weyl_zeroaware_volume_migration_experiment as ZVM
from pl_dual_complex import DualComplex,seed_16cell_boundary
from collective_barycentric_E_boundary_support_gate import barycentric_with_parent
from pl_peter_weyl_euclidean_local import LocalPLPeterWeylEuclidean

JMAX2=3
TOL=1e-10


def add(dst,src,scale):
    for k,a in src.items():
        z=dst.get(k,0j)+scale*a
        if abs(z)>TOL:
            dst[k]=z
        elif k in dst:
            del dst[k]


def signature(key):
    spins,Ks=key
    return {
        'spin_changes':[[i,int(s)] for i,s in enumerate(spins) if s!=1],
        'K_changes':[[i,int(K)] for i,K in enumerate(Ks) if K!=0],
    }


def run(local_index,parent_id=0):
    ZVM.patch_and_clear()
    coarse=seed_16cell_boundary()
    fine,parent=barycentric_with_parent(coarse)
    D=DualComplex(fine)
    G=LocalPLPeterWeylEuclidean(D)
    inside=sorted(v for v,p in enumerate(parent) if p==parent_id)
    if len(inside)!=24:
        raise RuntimeError(('parent block size',len(inside)))
    u=inside[local_index]
    seed=((1,)*len(G.EDGES),(0,)*D.n_tets)

    q4=[]
    source_disjoint=[]
    for sign,spec in G.oriented_specs(u):
        v,ra,rb,rc=spec
        path=D.plaquette_path(v,ra,rb)
        q=len(path)-1
        if q!=4:
            continue
        pe={tuple(sorted(e)) for e in zip(path[:-1],path[1:])}
        se=tuple(sorted((v,D.neighbor[(v,rc)])))
        source_disjoint.append(se not in pe)
        q4.append((sign,spec))

    col={}
    for sign,spec in q4:
        rr=dict(G.T_items(seed,*spec,JMAX2,False))
        aa=dict(G.T_items(seed,*spec,JMAX2,True))
        add(col,rr,-0.5j*sign)
        add(col,aa,+0.5j*sign)

    projected={
        k:a for k,a in col.items()
        if abs(a)>TOL and sum(s!=1 for s in k[0])==4
    }
    rows=[]
    for k,a in sorted(projected.items(),key=lambda kv:repr(kv[0])):
        sig=signature(k)
        rows.append({**sig,'re':float(a.real),'im':float(a.imag)})

    n=math.sqrt(sum(abs(a)**2 for a in projected.values()))
    finite=all(math.isfinite(a.real) and math.isfinite(a.imag) for a in projected.values())
    hard={
        'L1_nodes_384':D.n_tets==384,
        'L1_dual_links_768':len(G.EDGES)==768,
        'parent_has_24_fine_nodes':len(inside)==24,
        'six_q4_oriented_specs':len(q4)==6,
        'q4_source_link_disjoint_from_plaquette':all(source_disjoint),
        'projected_column_nonzero':len(projected)>0 and n>TOL,
        'all_projected_states_have_exactly_four_spin_changes':all(len(r['spin_changes'])==4 for r in rows),
        'finite_amplitudes':finite,
    }
    return {
        'status':'exact L1 q4-projected Euclidean amplitude column',
        'passed':bool(all(hard.values())),
        'science_status':'L1_BLOCK_E_Q4_EXACT_PROJECTION',
        'local_fine_index':local_index,
        'global_fine_node':u,
        'parent_coarse_tetra':parent_id,
        'L1_nodes':D.n_tets,
        'L1_dual_links':len(G.EDGES),
        'q4_oriented_specs':len(q4),
        'support':len(projected),
        'norm':n,
        'hard_checks':hard,
        'states':rows,
        'isolation_theorem':'Exactly-four-changed-edge projection receives q=4 curvature terms only; q>=6 terms necessarily change at least q plaquette edges from the all-j=1/2 seed.',
        'scope_note':'This proves a lower bound on the full fine-Hilbert Euclidean column span. It is not yet the coarse boundary isometry or the GR constraint rank.'
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--local-index',type=int,choices=range(24),required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    try:
        o=run(a.local_index);code=0 if o['passed'] else 1
    except Exception as exc:
        o={'status':'worker exception','passed':False,'science_status':'INFRASTRUCTURE_DIAGNOSTIC','local_fine_index':a.local_index,'error_type':type(exc).__name__,'error':str(exc),'traceback':traceback.format_exc()};code=1
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(o,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in o.items() if k!='states'},indent=2))
    return code

if __name__=='__main__':
    raise SystemExit(main())
