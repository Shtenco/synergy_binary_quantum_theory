#!/usr/bin/env python3
"""Collect the 12 direct-Hermitian physical Lorentzian pair contributions."""
from __future__ import annotations
import argparse,json,math,traceback,hashlib
from pathlib import Path
import numpy as np
import pl_lorentzian_48_collect as BASE

VERSION='direct-hermitian-commutator-v4';TOL=1e-10

def load_state(path):return BASE.load_state(path)
def add(dst,src):
    for k,a in src.items():dst[k]=dst.get(k,0j)+a
def norm(s):return math.sqrt(sum(abs(a)**2 for a in s.values()))
def max_spin(s):return max((max(k[0]) for k in s),default=0)/2.0

def run(root):
    rows=[];seen=set();S={}
    for p in Path(root).rglob('pair_*.json'):
        d=json.loads(p.read_text(encoding='utf-8'))
        if d.get('operator_version')!=VERSION:continue
        if d.get('passed') is not True:raise RuntimeError(f'v4 pair must explicitly pass with boolean true: {p}')
        for field in ('max_primitive_physical_leakage','max_primitive_nonscalar_rejected_norm'):
            value=d.get(field)
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(float(value)) or float(value)<0.0:
                raise RuntimeError(f'invalid {field} in {p}')
        key=(int(d['omitted_local_slot']),int(d['cyclic_index']))
        if key in seen:raise RuntimeError(f'duplicate v4 pair {key}')
        seen.add(key)
        cand=list(Path(root).rglob(f'pair_{key[0]}_{key[1]}.npz'))
        if len(cand)!=1:raise RuntimeError(f'need one NPZ for pair {key}, got {len(cand)}')
        expected_hash=d.get('state_sha256')
        if not isinstance(expected_hash,str) or len(expected_hash)!=64:
            raise RuntimeError(f'missing state_sha256 for pair {key}')
        actual_hash=hashlib.sha256(cand[0].read_bytes()).hexdigest()
        if actual_hash!=expected_hash:raise RuntimeError(f'pair NPZ/JSON SHA256 mismatch for {key}')
        st=load_state(cand[0]);add(S,st)
        rows.append({'omit':key[0],'cycle':key[1],'cyclic_order':d['cyclic_order'],'support':len(st),'norm':norm(st),
                     'max_spin':max_spin(st),'max_primitive_leakage':d['max_primitive_physical_leakage'],
                     'max_primitive_rejected':d['max_primitive_nonscalar_rejected_norm'],'runtime_seconds':d.get('runtime_seconds')})
    expect={(o,c) for o in range(4) for c in range(3)}
    if seen!=expect:raise RuntimeError(f'v4 pair orbit mismatch missing={sorted(expect-seen)} extra={sorted(seen-expect)}')
    S={k:a for k,a in S.items() if abs(a)>TOL}
    seed=((1,)*32,(0,)*16);seedpar=sum(seed[0])%2;wrong=sum(1 for k in S if sum(k[0])%2!=seedpar)
    finite=all(np.isfinite([z.real,z.imag]).all() for z in S.values())
    maxleak=max(r['max_primitive_leakage'] for r in rows);maxrej=max(r['max_primitive_rejected'] for r in rows)
    checks={'all_12_physical_pairs_loaded':len(rows)==12,'orbit_complete_unique':seen==expect,
            'finite_combined_amplitudes':finite,'max_primitive_leakage_below_1e-8':maxleak<1e-8,
            'max_primitive_nonscalar_rejection_below_1e-8':maxrej<1e-8,
            'S_preserves_even_valence_seed_parity':wrong==0,'single_L_spin_wall':max_spin(S)<=3.5+1e-12}
    out={'status':'exact direct-Hermitian 16-cell PL-S3 Lorentzian node column','passed':bool(all(checks.values())),
      'science_status':'AMPLITUDE_PRECURSOR_S_NODE0_V4','operator_version':VERSION,'source_node':0,'Jmax':3.5,
      'checks':checks,'pair_count':12,'primitive_orderings_covered':48,
      'S_support':len(S),'S_norm':norm(S),'S_max_spin':max_spin(S),'wrong_seed_parity_outputs':wrong,
      'max_primitive_physical_leakage':maxleak,'max_primitive_nonscalar_rejected_norm':maxrej,'pairs':sorted(rows,key=lambda r:(r['omit'],r['cycle'])),
      'definition':'sum of 12 exact -i/2 eta Tr_aux({[C_a(K),C_b(K)],C_c(V_tet)}) physical pair columns',
      'provenance_note':'No slot-orbit reconstruction. Each pair directly evaluates the same four V2 primitive orderings with shared caches and combines them before final collection.',
      'science_guard':'Nonzero S is reported, not required. Promotion to canonical C1 additionally requires the identity gate and at least one independently serialized V2 primitive-pair comparison.'}
    return S,out

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--json-output',type=Path,required=True);p.add_argument('--state-output',type=Path,required=True);a=p.parse_args()
    try:S,out=run(a.root);code=0 if out['passed'] else 1
    except Exception as exc:S={};out={'status':'v4 collector exception','passed':False,'operator_version':VERSION,'error_type':type(exc).__name__,'error':str(exc),'traceback':traceback.format_exc()};code=1
    a.json_output.parent.mkdir(parents=True,exist_ok=True);a.json_output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8');BASE.save_bundle(a.state_output,{}, {},S);print(json.dumps(out,indent=2,sort_keys=True));return code
if __name__=='__main__':raise SystemExit(main())
