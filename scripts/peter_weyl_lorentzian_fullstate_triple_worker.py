#!/usr/bin/env python3
"""Exact full-state worker for one ordered Lorentzian triple.

Computes the genuine covariant state

    Tr_aux[C_a(K_sine) C_b(K_sine) C_c(V)] |g>

at the frozen single-H_L wall Jmax=7/2, then performs only the exact scalar
covariant->Gauss closure. No logical projection is used. The 24 epsilon-
oriented terms are independent and can therefore be evaluated in parallel and
assembled linearly by a separate collector.

An individual ordered triple is allowed to vanish exactly. Nonzero support is
not a physical acceptance criterion; only finite arithmetic, scalar closure,
leakage control and cutoff safety are required before the 24-term sum.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0,str(HERE))

import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_lorentzian_logical_projection_gate as LP
import peter_weyl_lorentzian_gauss_action_gate as LGA

JMAX2=7
TOL=1e-10


def norm2(state):
    return float(sum(abs(a)**2 for a in state.values()))


def max_spin(state):
    return max((max(k[0]) for k in state),default=0)/2.0


def save_state(path,state):
    rows=sorted(state.items(),key=lambda kv:repr(kv[0]))
    if rows:
        spins=np.asarray([k[0] for k,_ in rows],dtype=np.int16)
        Ks=np.asarray([k[1] for k,_ in rows],dtype=np.int16)
        amp=np.asarray([a for _,a in rows],dtype=np.complex128)
    else:
        spins=np.zeros((0,len(PW.EDGES)),dtype=np.int16)
        Ks=np.zeros((0,len(PW.VERT)),dtype=np.int16)
        amp=np.zeros((0,),dtype=np.complex128)
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,spins=spins,Ks=Ks,amp=amp)


def run(a,b,c,coefficient,source=0,column_index=0):
    neigh=PW.NEIG[source]
    if len({a,b,c})!=3 or any(x not in neigh for x in (a,b,c)):
        raise ValueError(f'(a,b,c) must be distinct neighbors of source {source}: {neigh}')
    if coefficient not in (-1,1):
        raise ValueError('coefficient must be +/-1')

    basis=PW.basis_full_jhalf()
    if not (0<=column_index<len(basis)):
        raise ValueError(f'column_index must be in [0,{len(basis)-1}], got {column_index}')
    LP.JMAX2=JMAX2
    initial=basis[column_index]
    restore,caches=LP.install_sine_cached_stack()
    try:
        cov,diag=LP.ordered_triple_state(initial,source,a,b,c)
        gauss,accepted2,rejected2=LGA.project_scalar_gauss(cov,source,TOL)
        physical=max(
            float(diag.get('CV_complete_basis_leakage',0.0)),
            float(diag.get('CK_outer_complete_basis_leakage',0.0)),
            float(diag.get('CK_internal_volume_sector_leakage',0.0)),
        )
        charge=float(diag.get('CK_complete_charge_basis_leakage',0.0))
        total_scalar_test=accepted2+rejected2
        scalar_fraction=1.0 if total_scalar_test<1e-30 else accepted2/total_scalar_test
        cache_info={
            name:{'hits':fn.cache_info().hits,'misses':fn.cache_info().misses,'currsize':fn.cache_info().currsize}
            for name,fn in caches.items()
        }
        covnorm=math.sqrt(norm2(cov)); gnorm=math.sqrt(norm2(gauss))
        checks={
            'finite_covariant_norm':math.isfinite(covnorm),
            'finite_gauss_norm':math.isfinite(gnorm),
            'finite_covariant_amplitudes':all(np.isfinite([z.real,z.imag]).all() for z in cov.values()),
            'finite_gauss_amplitudes':all(np.isfinite([z.real,z.imag]).all() for z in gauss.values()),
            'physical_basis_volume_leakage':physical<1e-8,
            'scalar_closure_fraction_or_exact_zero':scalar_fraction>1-1e-10,
            'nonscalar_rejected_norm':math.sqrt(max(rejected2,0.0))<1e-8,
            'single_HL_spin_wall':max_spin(gauss)<=JMAX2/2+1e-12,
        }
        return gauss,{
            'status':'exact full-state ordered Lorentzian triple worker',
            'passed':all(checks.values()),
            'source_node':source,
            'ordered_edges':[a,b,c],
            'epsilon_coefficient':coefficient,
            'Jmax':JMAX2/2,
            'input_key':repr(initial),
            'column_index':int(column_index),
            'input_Ks':list(initial[1]),
            'source_pair_K01':list(initial[1][:2]),
            'covariant_support':len(cov),
            'covariant_norm':covnorm,
            'gauss_support':len(gauss),
            'gauss_norm':gnorm,
            'gauss_max_spin':max_spin(gauss),
            'exact_zero_ordered_term':len(cov)==0 and len(gauss)==0,
            'scalar_closure_fraction':scalar_fraction,
            'nonscalar_rejected_norm':math.sqrt(max(rejected2,0.0)),
            'physical_acceptance_max_leakage':physical,
            'historical_charge_diagnostic':charge,
            'cache_info':cache_info,
            'checks':checks,
            'weighted_here':False,
            'scope':'One of 24 exact ordered full-state terms; exact zero is allowed and the final epsilon sum is produced only by the collector.',
        }
    finally:
        restore()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--a',type=int,required=True); p.add_argument('--b',type=int,required=True); p.add_argument('--c',type=int,required=True)
    p.add_argument('--coefficient',type=int,required=True); p.add_argument('--source',type=int,default=0)
    p.add_argument('--json-output',type=Path,required=True); p.add_argument('--state-output',type=Path,required=True)
    p.add_argument('--column-index',type=int,default=0,help='Index in PW.basis_full_jhalf() (0,8,16,24 for K0/K1 logical pair, others 0)')
    x=p.parse_args(); state,out=run(x.a,x.b,x.c,x.coefficient,x.source,x.column_index)
    x.json_output.parent.mkdir(parents=True,exist_ok=True)
    x.json_output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    save_state(x.state_output,state)
    print(json.dumps(out,indent=2,sort_keys=True))
    return 0 if out['passed'] else 1

if __name__=='__main__': raise SystemExit(main())
