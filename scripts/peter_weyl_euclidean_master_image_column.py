#!/usr/bin/env python3
"""Produce one raw-reference Euclidean master image Y_i=M_E|b_i>.

The frozen safe Peter-Weyl Euclidean node Hamiltonian is implemented by
explicit primitive/adjoint symmetrization. At the declared finite reference
convention H_v^dagger=H_v, hence

    M_E |b_i> = sum_v H_v^2 |b_i>.

This producer intentionally does not call compose_on_sparse, whose final 1e-8
prune would turn the history moment into a retained-map composition.  Instead
both hits use apply_H_cached_state directly and the five node contributions are
summed without any additional tolerance prune.
"""
from __future__ import annotations

import argparse, gc, hashlib, json, math, sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW

JMAX2=5


def add_exact(dst,src):
    for k,z in src.items(): dst[k]=dst.get(k,0j)+complex(z)
    for k in [k for k,z in dst.items() if z==0j]: del dst[k]


def encode(state):
    return [
        {'spins':[int(x) for x in spins],'K_labels':[int(x) for x in Ks],
         'amp':[float(complex(z).real),float(complex(z).imag)]}
        for (spins,Ks),z in sorted(state.items(),key=lambda kv:repr(kv[0]))
    ]


def norm(s): return math.sqrt(sum(abs(z)**2 for z in s.values()))
def max_spin(s): return max((max(k[0])/2.0 for k in s),default=0.0)

def canonical_hash(rows):
    return hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def run(input_index:int):
    basis=PW.basis_full_jhalf()
    if not (0<=input_index<len(basis)): raise ValueError('input_index must be 0..31')
    initial=basis[input_index]
    total={};nodes=[];finite=True;maxsp=0.0
    operator_source_sha=hashlib.sha256(Path(PW.__file__).read_bytes()).hexdigest()

    PW.T_cached.cache_clear()
    for v in range(5):
        first=PW.apply_H_cached_state({initial:1+0j},v,JMAX2)
        second=PW.apply_H_cached_state(first,v,JMAX2)
        add_exact(total,second)
        finite &= all(np.isfinite(complex(z).real) and np.isfinite(complex(z).imag) for z in first.values())
        finite &= all(np.isfinite(complex(z).real) and np.isfinite(complex(z).imag) for z in second.values())
        maxsp=max(maxsp,max_spin(first),max_spin(second))
        nodes.append({
            'node':v,
            'first_raw_support':len(first),'first_raw_norm':norm(first),'first_raw_max_spin':max_spin(first),
            'second_raw_support':len(second),'second_raw_norm':norm(second),'second_raw_max_spin':max_spin(second),
        })
        # Cache clearing is performance-only and occurs after the complete H_v^2
        # contribution for this node has been formed.
        PW.T_cached.cache_clear();gc.collect()

    boundary_return=[]
    for j,key in enumerate(basis):
        z=total.get(key,0j)
        if z!=0j: boundary_return.append({'boundary_index':j,'amp':[float(z.real),float(z.imag)]})
    rows=encode(total)
    hard={
        'five_node_square_contributions':len(nodes)==5,
        'all_sparse_amplitudes_finite':bool(finite),
        'spin_cutoff_respected':maxsp<=2.5+1e-12,
        'no_compose_on_sparse_final_1e-8_prune':True,
        'no_additional_post_second_hit_tolerance_prune':True,
    }
    return {
        'schema':'BQG_EUCLIDEAN_MASTER_IMAGE_COLUMN_V1','passed':bool(all(hard.values())),
        'status':'raw-reference Euclidean master image M_E|b_i>',
        'input_index':input_index,'input_K_labels':list(initial[1]),'Jmax':2.5,
        'operator_identity':'M_E b_i = sum_v H_v^dagger H_v b_i = sum_v H_v^2 b_i in the frozen explicitly Hermitian Euclidean reference convention',
        'operator_source_sha256':operator_source_sha,
        'node_contributions':nodes,'support':len(total),'norm':norm(total),'max_spin_reached':maxsp,
        'boundary_return':boundary_return,'state_sha256':canonical_hash(rows),'hard_integrity_checks':hard,'state':rows,
        'claim_boundary':'One raw-reference Euclidean master-image column only. Full mu2 requires all 32 columns plus the pruning-bound consistency check against the certified one-hit packet. No Lorentzian/HDA/P_phys claim follows.'
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input-index',type=int,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    out=run(a.input_index);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,separators=(',',':'))+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='state'},indent=2));return 0 if out['passed'] else 1
if __name__=='__main__': raise SystemExit(main())
