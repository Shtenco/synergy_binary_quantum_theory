#!/usr/bin/env python3
"""Compute only the edge-dependent -h_a K h_a^-1 part of one outer-a shard.

The globally grouped Lorentzian first-column calculation separates

    C_a(K)_ij = delta_ij K - (h_a K h_a^-1)_ij.

This worker evaluates the second term for all four auxiliary (i,j) channels
after the exact signed prefix/k grouping for one fixed outer edge a.  The
outer-edge-independent direct K contribution is computed by a separate global
worker, allowing the eight repeated diagonal direct-K evaluations to be
replaced by only two global diagonal evaluations.

No symmetry assumption and no amplitude pruning is used before the hKh action.
"""
from __future__ import annotations

import argparse, hashlib, json
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
def sha_json(obj): return hashlib.sha256(json.dumps(obj,separators=(',',':'),sort_keys=True).encode()).hexdigest()


def load(prefix_dir:Path):
    s=json.loads((prefix_dir/'middle_prefix_summary.json').read_text(encoding='utf-8'))
    if s.get('schema')!='BQG_LORENTZIAN_ALL_MIDDLE_PREFIXES_V1' or s.get('passed') is not True: raise RuntimeError('invalid all-middle summary')
    ps=[]
    for idx in range(12):
        p=json.loads((prefix_dir/f'prefix_{idx}.json').read_text(encoding='utf-8'))
        if p.get('passed') is not True or int(p['ordered_pair']['pair_index'])!=idx: raise RuntimeError(f'invalid prefix {idx}')
        ps.append(p)
    return s,ps


def grouped_input(prefixes,source,a,i,j):
    g={};triples=[]
    for p in prefixes:
        b=int(p['ordered_pair']['b']);c=int(p['ordered_pair']['c'])
        rows=[r for r in PLAN.epsilon_outer_terms(source,b,c) if int(r['a'])==a]
        if not rows: continue
        if len(rows)!=1: raise RuntimeError('duplicate fixed-a epsilon row')
        row=rows[0];triples.append(tuple(int(x) for x in row['ordered_edges']));sign=int(row['sign'])
        pm={tuple(int(x) for x in q['indices']):q for q in p['paths']}
        for k in range(2): add_exact(g,scalar(decode(pm[(i,j,k)]['middle_state'])),scale=sign)
    return g,triples


def hkh_only(state,source,target,i,j,Jmax2):
    KCOMP=FULL.RAW.KCOMP;COMP=FULL.RAW.COMP
    hKh={};max_outer=0.0;max_v=0.0;max_b=0.0
    for k in range(2):
        inv,oleak=COMP.inverse_complete(state,source,target,k,j,Jmax2);max_outer=max(max_outer,float(oleak))
        Kinv,vleak,bleak=KCOMP.apply_K_complete_custom(inv,source,Jmax2,(source,target));max_v=max(max_v,float(vleak));max_b=max(max_b,float(bleak))
        COMP.add(hKh,COMP.close_complete(Kinv,source,target,i,k,Jmax2))
    out={};COMP.add(out,hKh,-1)
    return out,{'outer_complete_basis_leakage':max_outer,'internal_volume_sector_leakage':max_v,'complete_charge_basis_leakage':max_b}


def run(prefix_dir:Path,outer_a:int,middle_run_id=None,middle_head_sha=None):
    ZVM.patch_and_clear();summary,prefixes=load(prefix_dir);source=int(summary['source_node']);input_index=int(summary['input_logical_basis_index']);outer_a=int(outer_a)
    neighbors=tuple(FULL.RAW.PW.NEIG[source]);pairs=MID.ordered_pairs(source)
    if outer_a not in neighbors or len(pairs)!=12: raise RuntimeError('bad source/outer edge')
    if any(int(p['source_node'])!=source or int(p['input_logical_basis_index'])!=input_index for p in prefixes): raise RuntimeError('prefix provenance mismatch')
    basis=FULL.RAW.PW.basis_full_jhalf();initial=basis[input_index]

    groups={};triple_union=set()
    for i in range(2):
        for j in range(2):
            g,t=grouped_input(prefixes,source,outer_a,i,j);groups[(i,j)]=g;triple_union.update(t)
    if len(triple_union)!=6: raise RuntimeError(f'outer a={outer_a}: expected six ordered triples')

    total={};channels=[];diag={'outer_complete_basis_leakage':0.0,'internal_volume_sector_leakage':0.0,'complete_charge_basis_leakage':0.0};max_spin=0.0;calls=0
    old,caches=FULL.install_sine_ordering()
    try:
        for (i,j),g in groups.items():
            max_spin=max(max_spin,FULL.max_spin(g))
            if g:
                s3,d3=hkh_only(g,source,outer_a,i,j,FULL.JMAX2);calls+=1;FULL.update_diag(diag,d3)
            else: s3={}
            add_exact(total,s3);max_spin=max(max_spin,FULL.max_spin(s3))
            channels.append({'indices_ij':[i,j],'grouped_input_support':len(g),'grouped_input_norm':FULL.norm(g),'minus_hKh_support':len(s3),'minus_hKh_norm':FULL.norm(s3)})
        cache_info={name:{'hits':f.cache_info().hits,'misses':f.cache_info().misses,'currsize':f.cache_info().currsize} for name,f in caches.items()}
    finally: FULL.restore_ordering(old)

    sd=FULL.scalar_diagnostics(total);conv=PLAN.convention_descriptor(source);hab=PLAN.habitat_descriptor(source)
    hard={'six_outer_a_triples_accounted':len(triple_union)==6,'four_aux_channels':len(channels)==4,'hKh_calls_at_most_four':calls<=4,'no_pre_outer_tolerance_pruning':True,
          'outer_complete_basis_leakage_below_1e-9':diag['outer_complete_basis_leakage']<1e-9,'internal_volume_sector_leakage_below_1e-9':diag['internal_volume_sector_leakage']<1e-9,
          'minus_hKh_output_scalar_within_frozen_threshold':FULL.scalar_ok(sd),'spin_cutoff_respected':max_spin<=FULL.JMAX2/2+1e-12}
    return {'schema':'BQG_LORENTZIAN_GLOBAL_HKH_OUTER_A_V1','passed':bool(all(hard.values())),'science_status':'MINUS_HKH_OUTER_A_ZERO' if not total else 'MINUS_HKH_OUTER_A_NONZERO',
            'execution_mode':'global_signed_prefix_k_grouping_hKh_only_v1','direct_K_part_included':False,'source_node':source,'input_logical_basis_index':input_index,'input_K_labels':list(initial[1]),
            'outer_a':outer_a,'Jmax':FULL.JMAX2/2,'all_middle_provenance':{'run_id':middle_run_id,'head_sha':middle_head_sha,'summary_sha256':sha_json(summary)},
            'accounted_ordered_triples':[list(x) for x in sorted(triple_union)],'hKh_call_count':calls,'channels':channels,'partial_support':len(total),'partial_norm':FULL.norm(total),
            'partial_scalar_diagnostics':sd,'max_spin_reached':max_spin,'max_diagnostics':diag,'runtime_exact_cache':cache_info,'habitat_descriptor':hab,'habitat_hash':PLAN.canonical_hash(hab),
            'boundary_domain_hash':PLAN.boundary_domain_hash(basis),'convention_descriptor':conv,'convention_hash':PLAN.canonical_hash(conv),'hard_integrity_checks':hard,'state':PLAN.encode_state(total),
            'claim_boundary':'Only the signed edge-dependent -h_a K h_a^-1 contribution for one outer edge of one raw Lorentzian boundary column. The global direct-K term is deliberately absent and must be supplied by the matching global-direct packet before reconstructing the candidate H_L column.'}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--middle-dir',type=Path,required=True);ap.add_argument('--outer-a',type=int,required=True);ap.add_argument('--middle-run-id',type=int);ap.add_argument('--middle-head-sha');ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    out=run(a.middle_dir,a.outer_a,a.middle_run_id,a.middle_head_sha);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in out.items() if k not in ('state','channels')},indent=2));return 0 if out['passed'] else 1
if __name__=='__main__': raise SystemExit(main())
