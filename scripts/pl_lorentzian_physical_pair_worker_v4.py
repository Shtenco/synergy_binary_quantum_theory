#!/usr/bin/env python3
"""Direct physical Hermitian PL Lorentzian pair worker (v4 execution path).

For one omitted local slot and one cyclic order (a,b,c), compute the four V2
primitive columns needed for

  S_pair = -i/2 * eta * [ F_abc - F_bac - Fdag_abc + Fdag_bac ]

in one process with the same Peter-Weyl caches and tetrahedral charged-volume
backend.  This is exactly

  -i/2 eta Tr_aux({[C_a(K),C_b(K)], C_c(V_tet)}) |Omega>.

Every primitive retains the frozen V2 acceptance checks.  No tetrahedral orbit
reconstruction or GR target is used.
"""
from __future__ import annotations
import argparse,itertools,json,math,traceback,time,hashlib
from pathlib import Path
import numpy as np
import peter_weyl_lorentzian_logical_projection_gate as LP
import pl_lorentzian_triple_worker as W
import pl_lorentzian_group_worker as GW
from pl_dual_complex import DualComplex,seed_16cell_boundary
from pl_peter_weyl_euclidean import PLPeterWeylEuclidean
from pl_covariant_backend import install_pl_graph
from tetrahedral_volume_backend import install_tetrahedral_volume_backend

JMAX2=W.JMAX2;TOL=W.TOL;VERSION='direct-hermitian-commutator-v4'

def add(dst,src,scale=1.0,tol=1e-10):
    for k,a in src.items():
        z=dst.get(k,0j)+scale*a
        if abs(z)>tol:dst[k]=z
        elif k in dst:del dst[k]
def norm(s):return math.sqrt(sum(abs(z)**2 for z in s.values()))
def max_spin(s):return max((max(k[0]) for k in s),default=0)/2.0

def pair_indices(D,source,omit,cycle):
    base=tuple(r for r in range(4) if r!=omit)
    cyc=(base,(base[1],base[2],base[0]),(base[2],base[0],base[1]))
    p=cyc[cycle];q=(p[1],p[0],p[2])
    lookup={W.ordered_spec(D,source,i)[2]:i for i in range(24)}
    return p,q,lookup[p],lookup[q]

def run(omit,cycle,source=0):
    t0=time.time();D=DualComplex(seed_16cell_boundary());G=PLPeterWeylEuclidean(D)
    seed=((1,)*len(G.EDGES),(0,)*D.n_tets);p,q,ip,iq=pair_indices(D,source,omit,cycle)
    oldJ=LP.JMAX2;states={};meta={}
    with install_tetrahedral_volume_backend():
      with install_pl_graph(G):
        LP.JMAX2=JMAX2;restore,caches=LP.install_sine_cached_stack()
        try:
          for mode,idx,label in [('forward',ip,'fp'),('forward',iq,'fq'),('adjoint',ip,'ap'),('adjoint',iq,'aq')]:
            st,m=GW.evaluate(D,G,seed,source,idx,mode,caches);states[label]=st;meta[label]=m
        finally:
          restore();LP.JMAX2=oldJ
    cp=int(meta['fp']['PL_epsilon_coefficient']);cq=int(meta['fq']['PL_epsilon_coefficient'])
    S={}
    add(S,states['fp'],-0.5j*cp);add(S,states['fq'],-0.5j*cq)
    add(S,states['ap'],+0.5j*cp);add(S,states['aq'],+0.5j*cq)
    primitive=[meta[x] for x in ('fp','fq','ap','aq')]
    finite=all(np.isfinite([z.real,z.imag]).all() for z in S.values())
    maxleak=max(float(x['physical_acceptance_max_leakage']) for x in primitive)
    maxrej=max(float(x['nonscalar_rejected_norm']) for x in primitive)
    checks={
      'all_four_v2_primitives_pass':all(x['passed'] for x in primitive),
      'paired_indices_distinct':ip!=iq,
      'epsilon_coefficients_opposite':cp==-cq,
      'same_omitted_slot':all(int(x['omitted_local_slot'])==omit for x in primitive),
      'finite_physical_pair_amplitudes':finite,
      'primitive_physical_leakage_below_1e-8':maxleak<1e-8,
      'primitive_nonscalar_rejection_below_1e-8':maxrej<1e-8,
      'single_L_spin_wall':max_spin(S)<=JMAX2/2+1e-12,
    }
    out={'status':'direct physical Hermitian PL Lorentzian pair','passed':bool(all(checks.values())),
      'science_status':'AMPLITUDE_PRECURSOR_S_PAIR_V4','operator_version':VERSION,
      'volume_definition':'V_tet=sqrt(abs((1/4) sum_r (-1)^r q_hat_r)) with production zero-aware spectrum',
      'source_node':source,'omitted_local_slot':omit,'cyclic_index':cycle,
      'cyclic_order':list(p),'first_two_swap':list(q),'forward_indices':[ip,iq],'adjoint_indices':[ip,iq],
      'epsilon_coefficients':[cp,cq],'Jmax':JMAX2/2,'checks':checks,
      'S_pair_support':len(S),'S_pair_norm':norm(S),'S_pair_max_spin':max_spin(S),
      'max_primitive_physical_leakage':maxleak,'max_primitive_nonscalar_rejected_norm':maxrej,
      'primitive_metadata':primitive,'runtime_seconds':time.time()-t0,
      'definition':'S_pair=-i/2*(coef_p F_p + coef_q F_q - coef_p Fdag_p - coef_q Fdag_q)',
      'identity':'=-i/2 eta Tr_aux({[C_a(K),C_b(K)],C_c(V_tet)})|Omega>',
      'scope_note':'One of 12 exact physical contributions. No slot-orbit covariance or reconstructed heavy term is used.'}
    return S,out

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--omit',type=int,choices=range(4),required=True);p.add_argument('--cycle',type=int,choices=range(3),required=True);p.add_argument('--source',type=int,default=0);p.add_argument('--json-output',type=Path,required=True);p.add_argument('--state-output',type=Path,required=True);a=p.parse_args()
    try:S,out=run(a.omit,a.cycle,a.source);code=0 if out['passed'] else 1
    except Exception as exc:S={};out={'status':'v4 pair worker exception','passed':False,'operator_version':VERSION,'omit':a.omit,'cycle':a.cycle,'error_type':type(exc).__name__,'error':str(exc),'traceback':traceback.format_exc()};code=1
    D=DualComplex(seed_16cell_boundary());W.save_state(a.state_output,S,len(D.dual_edges()),D.n_tets)
    if out.get('passed') is True:
        out['state_sha256']=hashlib.sha256(a.state_output.read_bytes()).hexdigest()
    a.json_output.parent.mkdir(parents=True,exist_ok=True);a.json_output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8');print(json.dumps(out,indent=2,sort_keys=True));return code
if __name__=='__main__':raise SystemExit(main())
