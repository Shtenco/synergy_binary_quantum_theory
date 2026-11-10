#!/usr/bin/env python3
"""Cache-sharing six-term shard for the preregistered 16-cell Lorentzian orbit.

All 24 terms remain individually evaluated, checked and serialized.  The only
change from 24 independent processes is resource scheduling:
- forward shards group by c, the first C(V) acting on the ket;
- adjoint shards group by a, the first C(K) acting on the ket.
This reuses exact frozen-state C(V)/C(K) caches without changing any amplitude,
orientation coefficient, truncation or acceptance criterion.
"""
from __future__ import annotations
import argparse,json,math,time,traceback
from pathlib import Path
import numpy as np
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_lorentzian_logical_projection_gate as LP
import peter_weyl_lorentzian_gauss_action_gate as LGA
from pl_dual_complex import DualComplex,seed_16cell_boundary
from pl_peter_weyl_euclidean import PLPeterWeylEuclidean
from pl_covariant_backend import install_pl_graph
import pl_lorentzian_triple_worker as W

JMAX2=W.JMAX2;TOL=W.TOL

from v2_term_guards import validate_term_diagnostics


def term_indices(D,source,mode,slot):
    out=[]
    for idx in range(24):
        _,_,perm,_,_=W.ordered_spec(D,source,idx)
        first=perm[-1] if mode=='forward' else perm[0]
        if first==slot:out.append(idx)
    if len(out)!=6:raise RuntimeError((mode,slot,out))
    return out

def evaluate(D,G,seed,source,idx,mode,caches):
    omit,base,perm,targets,coef=W.ordered_spec(D,source,idx);t0=time.time()
    if mode=='forward':cov,diag=LP.ordered_triple_state(seed,source,*targets)
    else:cov,diag=W.ordered_dagger_state(seed,source,*targets)
    gauss,accepted2,rejected2=LGA.project_scalar_gauss(cov,source,TOL)
    validated=validate_term_diagnostics(diag,accepted2,rejected2)
    accepted2=validated['accepted2'];rejected2=validated['rejected2']
    total=accepted2+rejected2;frac=1.0 if total<1e-30 else accepted2/total
    physical=max(validated['CV_complete_basis_leakage'],
                 validated['CK_outer_complete_basis_leakage'],
                 validated['CK_internal_volume_sector_leakage'])
    charge=validated['CK_complete_charge_basis_leakage']
    cn=math.sqrt(W.norm2(cov));gn=math.sqrt(W.norm2(gauss))
    checks={
      'finite_covariant_norm':math.isfinite(cn),'finite_gauss_norm':math.isfinite(gn),
      'finite_covariant_amplitudes':all(np.isfinite([z.real,z.imag]).all() for z in cov.values()),
      'finite_gauss_amplitudes':all(np.isfinite([z.real,z.imag]).all() for z in gauss.values()),
      'physical_basis_volume_leakage':physical<1e-8,
      'scalar_closure_fraction_or_exact_zero':frac>1-1e-10,
      'nonscalar_rejected_norm':math.sqrt(max(rejected2,0.0))<1e-8,
      'single_L_spin_wall':W.max_spin(gauss)<=JMAX2/2+1e-12,
      'PL_orientation_coefficient':coef in (-1,1)}
    cache_info={name:{'hits':fn.cache_info().hits,'misses':fn.cache_info().misses,'currsize':fn.cache_info().currsize} for name,fn in caches.items()}
    meta={'status':'exact cache-sharing PL-S3 Lorentzian ordered term','passed':bool(all(checks.values())),
          'mode':mode,'index':idx,'source_node':source,'omitted_local_slot':omit,'base_local_slots':list(base),
          'permuted_local_slots':list(perm),'ordered_target_nodes':list(targets),'PL_epsilon_coefficient':coef,
          'Jmax':JMAX2/2,'input_key':repr(seed),'covariant_support':len(cov),'covariant_norm':cn,
          'gauss_support':len(gauss),'gauss_norm':gn,'gauss_max_spin':W.max_spin(gauss),
          'exact_zero_ordered_term':len(cov)==0 and len(gauss)==0,'scalar_closure_fraction':frac,
          'nonscalar_rejected_norm':math.sqrt(max(rejected2,0.0)),'physical_acceptance_max_leakage':physical,
          'complete_charge_diagnostic':charge,'cache_info_after_term':cache_info,'runtime_seconds':time.time()-t0,
          'checks':checks,'weighted_here':False,
          'scope_note':'One individually serialized term from a six-term cache-sharing process shard; mathematical result is identical to the independent worker definition.'}
    return gauss,meta

def run(mode,slot,outdir,source=0):
    D=DualComplex(seed_16cell_boundary());G=PLPeterWeylEuclidean(D);indices=term_indices(D,source,mode,slot)
    seed=((1,)*len(G.EDGES),(0,)*D.n_tets);oldJ=LP.JMAX2;summary=[];allpass=True
    with install_pl_graph(G):
        LP.JMAX2=JMAX2;restore,caches=LP.install_sine_cached_stack()
        try:
            for idx in indices:
                state,meta=evaluate(D,G,seed,source,idx,mode,caches);allpass&=meta['passed'];summary.append(meta)
                W.save_state(Path(outdir)/f'term_{mode}_{idx}.npz',state,len(G.EDGES),D.n_tets)
                (Path(outdir)/f'term_{mode}_{idx}.json').write_text(json.dumps(meta,indent=2,sort_keys=True)+'\n',encoding='utf-8')
        finally:
            restore();LP.JMAX2=oldJ
    return {'status':'six-term PL Lorentzian cache-sharing shard','passed':bool(allpass),'mode':mode,'first_actual_local_slot':slot,
            'indices':indices,'terms_passed':sum(x['passed'] for x in summary),'terms':summary,
            'execution_only_note':'Cache sharing changes no operator, coefficient, cutoff, pruning rule or hard threshold.'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mode',choices=('forward','adjoint'),required=True);p.add_argument('--first-slot',type=int,choices=range(4),required=True);p.add_argument('--source',type=int,default=0);p.add_argument('--out-dir',type=Path,required=True);a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    try:o=run(a.mode,a.first_slot,a.out_dir,a.source);code=0 if o['passed'] else 1
    except Exception as exc:o={'status':'shard exception','passed':False,'mode':a.mode,'first_actual_local_slot':a.first_slot,'error_type':type(exc).__name__,'error':str(exc),'traceback':traceback.format_exc()};code=1
    (a.out_dir/f'shard_{a.mode}_{a.first_slot}.json').write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='utf-8');print(json.dumps(o,indent=2,sort_keys=True));return code
if __name__=='__main__':raise SystemExit(main())
