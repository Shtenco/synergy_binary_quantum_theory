#!/usr/bin/env python3
"""Compute the globally factorized direct-K part of the first raw H_L column.

For the diagonal auxiliary channels of

    C_a(K)_ii = K - (h_a K h_a^-1)_ii,

the direct K operator is independent of the outer edge a.  Hence

    sum_a K S_{a,ii} = K (sum_a S_{a,ii}).

This worker constructs the signed all-a middle sums for i=0 and i=1 and applies
K only twice.  It does not assume that the pre-outer sums vanish; their actual
residuals and the resulting direct-K state are serialized.
"""
from __future__ import annotations

import argparse, hashlib, json, math
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_lorentzian_middle_prefix_gate as MID
import peter_weyl_lorentzian_pruned_prefix_worker as PLAN
import peter_weyl_zeroaware_volume_migration_experiment as ZVM


def decode(rows):
    out={}
    for r in rows:
        key=(tuple(int(x) for x in r['spins']),tuple(int(x) for x in r['Kother']),int(r['J2']),int(r['M2']),int(r['K12']),int(r['K34']))
        out[key]=out.get(key,0j)+complex(float(r['amp'][0]),float(r['amp'][1]))
    return out


def add_exact(dst,src,scale=1):
    for k,z in src.items(): dst[k]=dst.get(k,0j)+complex(scale)*z
    for k in [k for k,z in dst.items() if z==0j]: del dst[k]


def scalar(state): return {k:z for k,z in state.items() if int(k[2]) in (0,2)}
def norm(s): return math.sqrt(sum(abs(z)**2 for z in s.values()))
def sha_json(obj): return hashlib.sha256(json.dumps(obj,separators=(',',':'),sort_keys=True).encode()).hexdigest()


def run(prefix_dir:Path,middle_run_id=None,middle_head_sha=None):
    ZVM.patch_and_clear();summary=json.loads((prefix_dir/'middle_prefix_summary.json').read_text(encoding='utf-8'))
    if summary.get('schema')!='BQG_LORENTZIAN_ALL_MIDDLE_PREFIXES_V1' or summary.get('passed') is not True: raise RuntimeError('invalid all-middle summary')
    source=int(summary['source_node']);input_index=int(summary['input_logical_basis_index']);neighbors=tuple(FULL.RAW.PW.NEIG[source]);pairs=MID.ordered_pairs(source)
    if len(neighbors)!=4 or len(pairs)!=12: raise RuntimeError('bad source valence/prefix family')
    prefixes=[]
    for idx in range(12):
        p=json.loads((prefix_dir/f'prefix_{idx}.json').read_text(encoding='utf-8'))
        if p.get('passed') is not True or int(p['ordered_pair']['pair_index'])!=idx: raise RuntimeError(f'invalid prefix {idx}')
        prefixes.append(p)
    basis=FULL.RAW.PW.basis_full_jhalf();initial=basis[input_index]

    diagonal_groups={0:{},1:{}}
    for i in range(2):
        g=diagonal_groups[i]
        for a in neighbors:
            for p in prefixes:
                b=int(p['ordered_pair']['b']);c=int(p['ordered_pair']['c'])
                rows=[r for r in PLAN.epsilon_outer_terms(source,b,c) if int(r['a'])==a]
                if not rows: continue
                row=rows[0];sign=int(row['sign']);pm={tuple(int(x) for x in q['indices']):q for q in p['paths']}
                for k in range(2): add_exact(g,scalar(decode(pm[(i,i,k)]['middle_state'])),scale=sign)

    total={};rows=[];max_v=0.0;max_b=0.0;calls=0;max_spin=0.0
    old,caches=FULL.install_sine_ordering()
    try:
        for i,g in diagonal_groups.items():
            pre_norm=norm(g);pre_max=max((abs(z) for z in g.values()),default=0.0);max_spin=max(max_spin,FULL.max_spin(g))
            if g:
                out,vleak,bleak=FULL.RAW.KCOMP.direct_K_covariant(g,source,FULL.JMAX2);calls+=1;max_v=max(max_v,float(vleak));max_b=max(max_b,float(bleak))
            else: out={}
            add_exact(total,out);max_spin=max(max_spin,FULL.max_spin(out));rows.append({'diagonal_i':i,'global_preouter_support':len(g),'global_preouter_norm':pre_norm,'global_preouter_max_abs':pre_max,'direct_K_support':len(out),'direct_K_norm':FULL.norm(out)})
        cache_info={name:{'hits':f.cache_info().hits,'misses':f.cache_info().misses,'currsize':f.cache_info().currsize} for name,f in caches.items()}
    finally: FULL.restore_ordering(old)

    sd=FULL.scalar_diagnostics(total);conv=PLAN.convention_descriptor(source);hab=PLAN.habitat_descriptor(source)
    hard={'two_diagonal_channels':len(rows)==2,'direct_K_calls_at_most_two':calls<=2,'no_pre_direct_tolerance_pruning':True,'volume_sector_leakage_below_1e-9':max_v<1e-9,'direct_output_scalar_within_frozen_threshold':FULL.scalar_ok(sd),'spin_cutoff_respected':max_spin<=FULL.JMAX2/2+1e-12}
    return {'schema':'BQG_LORENTZIAN_GLOBAL_DIRECT_K_V1','passed':bool(all(hard.values())),'science_status':'GLOBAL_DIRECT_K_ZERO' if not total else 'GLOBAL_DIRECT_K_NONZERO','execution_mode':'factor_direct_K_across_all_outer_edges_v1',
            'source_node':source,'input_logical_basis_index':input_index,'input_K_labels':list(initial[1]),'Jmax':FULL.JMAX2/2,'direct_K_call_count':calls,'diagonal_channels':rows,
            'direct_support':len(total),'direct_norm':FULL.norm(total),'direct_scalar_diagnostics':sd,'max_spin_reached':max_spin,'max_volume_sector_leakage':max_v,
            'runtime_exact_cache':cache_info,'all_middle_provenance':{'run_id':middle_run_id,'head_sha':middle_head_sha,'summary_sha256':sha_json(summary)},'habitat_descriptor':hab,'habitat_hash':PLAN.canonical_hash(hab),
            'boundary_domain_hash':PLAN.boundary_domain_hash(basis),'convention_descriptor':conv,'convention_hash':PLAN.canonical_hash(conv),'hard_integrity_checks':hard,'state':PLAN.encode_state(total),
            'claim_boundary':'Globally factorized outer-edge-independent direct K contribution for the two diagonal auxiliary channels of one raw Lorentzian boundary column. It must be combined with all four matching -h_a K h_a^-1 packets; alone it is not H_L, HDA, P_phys or cosmology.'}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--middle-dir',type=Path,required=True);ap.add_argument('--middle-run-id',type=int);ap.add_argument('--middle-head-sha');ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    out=run(a.middle_dir,a.middle_run_id,a.middle_head_sha);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in out.items() if k!='state'},indent=2));return 0 if out['passed'] else 1
if __name__=='__main__': raise SystemExit(main())
