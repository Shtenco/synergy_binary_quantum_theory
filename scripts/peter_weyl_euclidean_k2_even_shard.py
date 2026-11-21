#!/usr/bin/env python3
"""Produce 32 reusable two-hit Euclidean states for one (target,source) shard.

Loads the certified one-hit packet g_(v,i)=H_v b_i and computes only
q_(w,v,i)=H_w g_(v,i), i=0..31.  No first-layer Hamiltonian action is repeated.
The complete sparse second-action states are serialized for later global K2-even
rank reveal and physical-history/block-Krylov work.
"""
from __future__ import annotations

import argparse, gc, hashlib, json, math, sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW

JMAX2=5
PRUNE=1.0e-8


def decode(rows):
    out={}
    for r in rows:
        key=(tuple(int(x) for x in r['spins']),tuple(int(x) for x in r['K_labels']))
        z=complex(float(r['amp'][0]),float(r['amp'][1]));out[key]=out.get(key,0j)+z
    return out


def encode(state):
    rows=[]
    for (spins,K),z in sorted(state.items()):
        rows.append({'spins':list(spins),'K_labels':list(K),'amp':[float(z.real),float(z.imag)]})
    return rows


def norm(s): return math.sqrt(sum(abs(z)**2 for z in s.values()))
def max_spin(s): return max((max(k[0])/2.0 for k in s),default=0.0)

def state_hash(rows):
    return hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def run(packet_dir:Path,target_node:int,source_node:int,clear_every:int=4):
    if target_node not in range(5) or source_node not in range(5): raise ValueError('target/source nodes must be 0..4')
    packet_dir=Path(packet_dir);manifest=json.loads((packet_dir/'euclidean_packet_manifest.json').read_text(encoding='utf-8'))
    if manifest.get('schema')!='BQG_MICROSCOPIC_CONSTRAINT_PACKET_V2' or manifest.get('family')!='E' or not manifest.get('passed',False):
        raise RuntimeError('invalid reusable Euclidean packet')
    if float(manifest.get('Jmax',-1))!=2.5: raise RuntimeError('unexpected Euclidean packet Jmax')

    basis=PW.basis_full_jhalf()
    if len(basis)!=32: raise RuntimeError('frozen 32D boundary basis unavailable')
    outputs=[];supports=[];norms=[];maxsp=[];finite=True;first_gate=True
    PW.T_cached.cache_clear()
    for i in range(32):
        p=packet_dir/'columns'/f'E_node{source_node}_input{i:02d}.json'
        d=json.loads(p.read_text(encoding='utf-8'))
        first_gate &= bool(d.get('passed',False))
        first=decode(d['complete_gauss_outgoing_column']['state'])
        second=PW.compose_on_sparse(first,target_node,JMAX2)
        rows=encode(second)
        finite &= all(np.isfinite(r['amp'][0]) and np.isfinite(r['amp'][1]) for r in rows)
        supports.append(len(second));norms.append(norm(second));maxsp.append(max_spin(second))
        boundary_return=[]
        for bi,key in enumerate(basis):
            z=second.get(key,0j)
            if z!=0j: boundary_return.append({'boundary_index':bi,'amp':[float(z.real),float(z.imag)]})
        outputs.append({
            'input_index':i,'input_K_labels':list(basis[i][1]),
            'first_source_node':source_node,'second_target_node':target_node,
            'support':len(second),'norm':norm(second),'max_spin':max_spin(second),
            'boundary_return':boundary_return,'state_sha256':state_hash(rows),'state':rows,
        })
        if clear_every>0 and (i+1)%clear_every==0:
            PW.T_cached.cache_clear();gc.collect()

    hard={
        'source_packet_passed':manifest.get('passed') is True,
        'all_32_source_columns_passed':bool(first_gate),
        'exact_32_second_actions':len(outputs)==32,
        'all_sparse_amplitudes_finite':bool(finite),
        'spin_cutoff_respected':max(maxsp,default=0.0)<=2.5+1e-12,
    }
    return {
        'schema':'BQG_EUCLIDEAN_K2_EVEN_SHARD_V1','passed':bool(all(hard.values())),
        'status':'reusable 32-column Euclidean two-hit even shard',
        'target_node':target_node,'source_node':source_node,'Jmax':2.5,'prune_threshold':PRUNE,
        'source_packet_sha256':manifest.get('packet_sha256'),'source_domain_label':manifest.get('domain_label'),
        'columns':32,'support_min':min(supports),'support_max':max(supports),'support_mean':float(np.mean(supports)),
        'norm_min':min(norms),'norm_max':max(norms),'norm_mean':float(np.mean(norms)),'max_spin_reached':max(maxsp),
        'hard_integrity_checks':hard,'states':outputs,
        'claim_boundary':'Reusable forward two-hit states H_w H_v b_i only. They form K2-even carrier data but are not identified with M_E V0, do not give mu2 by themselves, and do not emit P_phys.'
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--packet-dir',type=Path,required=True);ap.add_argument('--target-node',type=int,required=True);ap.add_argument('--source-node',type=int,required=True);ap.add_argument('--clear-cache-every',type=int,default=4);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    out=run(a.packet_dir,a.target_node,a.source_node,a.clear_cache_every);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,separators=(',',':'))+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='states'},indent=2));return 0 if out['passed'] else 1
if __name__=='__main__': raise SystemExit(main())
