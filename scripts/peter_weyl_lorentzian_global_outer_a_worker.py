#!/usr/bin/env python3
"""Exact global outer-edge continuation from the certified all-middle packet.

For one frozen outer edge a, regroup all surviving Lorentzian middle-prefix
states before the expensive final C_a(K) hit:

  sum_{b,c,k} eps_{abc} C_a(K)_{ij} Xi_bc^{ijk}
    = C_a(K)_{ij} [sum_{b,c,k} eps_{abc} Xi_bc^{ijk}].

Thus each outer-a worker performs at most four final C(K) component actions,
one for each (i,j), rather than repeating the same linear operator for every
prefix and k.  For off-diagonal auxiliary components the exact definition

  C_ij(K) = delta_ij K - (h K h^-1)_ij

has delta_ij=0, so the expensive direct-K branch is not evaluated for i!=j.
This changes no operator or amplitude; it only removes a computation whose
coefficient is exactly zero.  No tolerance pruning or symmetry reconstruction
is used.
"""
from __future__ import annotations

import argparse, hashlib, json
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_lorentzian_pruned_prefix_worker as PLAN
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

EXPECTED_ZERO=[2,5,8,9,10,11]
EXPECTED_LIVE=[0,1,3,4,6,7]


def decode(rows):
    out={}
    for r in rows:
        key=(tuple(int(x) for x in r['spins']),tuple(int(x) for x in r['Kother']),int(r['J2']),int(r['M2']),int(r['K12']),int(r['K34']))
        z=complex(float(r['amp'][0]),float(r['amp'][1]));out[key]=out.get(key,0j)+z
    return {k:z for k,z in out.items() if z!=0j}


def add_exact(dst,src,scale=1):
    for k,z in src.items(): dst[k]=dst.get(k,0j)+scale*z
    for k in [k for k,z in dst.items() if z==0j]: del dst[k]


def scalar_channel(state): return {k:z for k,z in state.items() if int(k[2]) in (0,2)}


def path_map(prefix):
    rows={tuple(int(x) for x in p['indices']):p for p in prefix.get('paths',[])}
    expected={(i,j,k) for i in range(2) for j in range(2) for k in range(2)}
    if set(rows)!=expected: raise RuntimeError('prefix does not contain the exact eight auxiliary paths')
    return rows


def sha_json(obj): return hashlib.sha256(json.dumps(obj,separators=(',',':'),sort_keys=True).encode()).hexdigest()


def C_K_component_trace_optimized(state_cov,source_v,target_v,i,j,Jmax2):
    """Frozen C_K_component with an exact-zero branch removed for i!=j.

    The reference implementation evaluates direct_K_covariant unconditionally
    and multiplies it by delta_ij only afterwards.  For i!=j the direct term is
    identically absent, so skipping its evaluation is algebraically exact.
    """
    KCOMP=FULL.RAW.KCOMP; COMP=FULL.RAW.COMP
    if i==j:
        out,diag=KCOMP.C_K_component(state_cov,source_v,target_v,i,j,Jmax2)
        return out,diag,False
    hKh={};max_outer=0.0;max_v=0.0;max_b=0.0
    for k in range(2):
        inv,oleak=COMP.inverse_complete(state_cov,source_v,target_v,k,j,Jmax2);max_outer=max(max_outer,oleak)
        Kinv,vleak,bleak=KCOMP.apply_K_complete_custom(inv,source_v,Jmax2,(source_v,target_v));max_v=max(max_v,vleak);max_b=max(max_b,bleak)
        COMP.add(hKh,COMP.close_complete(Kinv,source_v,target_v,i,k,Jmax2))
    out={};COMP.add(out,hKh,-1)
    return out,{'outer_complete_basis_leakage':max_outer,'internal_volume_sector_leakage':max_v,'complete_charge_basis_leakage':max_b},True


def run(middle_dir:Path,outer_a:int,source_v=0,input_index=0,middle_run_id=None,middle_head_sha=None):
    if outer_a not in FULL.RAW.PW.NEIG[source_v]: raise ValueError('outer edge is not a source neighbor')
    summary=json.loads((middle_dir/'middle_prefix_summary.json').read_text(encoding='utf-8'))
    if summary.get('schema')!='BQG_LORENTZIAN_ALL_MIDDLE_PREFIXES_V1' or summary.get('passed') is not True: raise RuntimeError('invalid all-middle summary')
    if int(summary.get('source_node',-1))!=source_v or int(summary.get('input_logical_basis_index',-1))!=input_index: raise RuntimeError('all-middle source/input mismatch')
    if summary.get('zero_prefix_indices')!=EXPECTED_ZERO or summary.get('nonzero_prefix_indices')!=EXPECTED_LIVE: raise RuntimeError('frozen zero/live partition mismatch')
    if abs(float(summary.get('Jmax',-1))-FULL.JMAX2/2)>1e-12: raise RuntimeError('Jmax mismatch')

    basis=FULL.RAW.PW.basis_full_jhalf();initial=basis[input_index];groups={(i,j):{} for i in range(2) for j in range(2)}
    live_triples=[];zero_triples=[];prefix_audit=[]
    for idx in range(12):
        p=json.loads((middle_dir/f'prefix_{idx}.json').read_text(encoding='utf-8'));pair=p.get('ordered_pair',{});b=int(pair.get('b',-1));c=int(pair.get('c',-1))
        if int(pair.get('pair_index',-1))!=idx: raise RuntimeError(f'prefix_{idx}: pair index mismatch')
        plan=[r for r in PLAN.epsilon_outer_terms(source_v,b,c) if int(r['a'])==outer_a]
        if len(plan)>1: raise RuntimeError('a prefix cannot have two identical outer edges')
        if not plan: continue
        row=plan[0];triple=tuple(int(x) for x in row['ordered_edges']);sign=int(row['sign']);is_zero=bool(p.get('prefix_zero_pathwise')) or p.get('science_status')=='MIDDLE_PREFIX_ZERO_PATHWISE'
        if idx in EXPECTED_ZERO:
            if not is_zero: raise RuntimeError(f'prefix {idx} must be pathwise zero')
            zero_triples.append(triple);prefix_audit.append({'pair_index':idx,'b':b,'c':c,'triple':list(triple),'sign':sign,'status':'ZERO_BEFORE_OUTER'});continue
        if idx not in EXPECTED_LIVE or is_zero or p.get('science_status')!='MIDDLE_PREFIX_NONZERO': raise RuntimeError(f'prefix {idx} live classification mismatch')
        pm=path_map(p)
        for i in range(2):
            for j in range(2):
                for k in range(2): add_exact(groups[(i,j)],scalar_channel(decode(pm[(i,j,k)].get('middle_state',[]))),scale=sign)
        live_triples.append(triple);prefix_audit.append({'pair_index':idx,'b':b,'c':c,'triple':list(triple),'sign':sign,'status':'GROUPED_BEFORE_OUTER'})
    if len(live_triples)+len(zero_triples)!=6: raise RuntimeError(f'outer a={outer_a}: expected exactly six epsilon triples')

    ZVM.patch_and_clear();total={};rows=[];diag={'CK_outer_complete_basis_leakage':0.0,'CK_internal_volume_sector_leakage':0.0,'CK_complete_charge_basis_leakage':0.0};max_spin=0.0;calls=0;offdiag_skips=0
    old,caches=FULL.install_sine_ordering()
    try:
        for (i,j),grouped in groups.items():
            max_spin=max(max_spin,FULL.max_spin(grouped))
            if grouped:
                s3,d3,skipped=C_K_component_trace_optimized(grouped,source_v,outer_a,i,j,FULL.JMAX2);calls+=1;offdiag_skips+=int(skipped);FULL.update_diag(diag,d3)
            else: s3={};skipped=False
            add_exact(total,s3);max_spin=max(max_spin,FULL.max_spin(s3));rows.append({'indices_ij':[i,j],'grouped_input_support':len(grouped),'grouped_input_norm':FULL.norm(grouped),'outer_output_support':len(s3),'outer_output_norm':FULL.norm(s3),'offdiagonal_direct_K_skipped':bool(skipped)})
        cache_info={name:{'hits':f.cache_info().hits,'misses':f.cache_info().misses,'currsize':f.cache_info().currsize} for name,f in caches.items()}
    finally: FULL.restore_ordering(old)

    scalar=FULL.scalar_diagnostics(total);hard={'six_outer_a_triples_accounted':len(live_triples)+len(zero_triples)==6,'live_zero_triples_disjoint':not(set(live_triples)&set(zero_triples)),'at_most_four_outer_CK_calls':calls<=4,'offdiagonal_direct_K_removed_only_by_exact_delta_ij':offdiag_skips<=2,'no_pre_outer_tolerance_pruning':True,'outer_complete_basis_leakage_below_1e-9':diag['CK_outer_complete_basis_leakage']<1e-9,'outer_internal_volume_sector_leakage_below_1e-9':diag['CK_internal_volume_sector_leakage']<1e-9,'partial_output_scalar_within_frozen_threshold':FULL.scalar_ok(scalar),'spin_cutoff_respected':max_spin<=FULL.JMAX2/2+1e-12}
    conv=PLAN.convention_descriptor(source_v);hab=PLAN.habitat_descriptor(source_v)
    return {'schema':'BQG_LORENTZIAN_GLOBAL_OUTER_A_V1','passed':bool(all(hard.values())),'science_status':'OUTER_A_ZERO' if not total else 'OUTER_A_NONZERO','execution_mode':'exact_global_prefix_k_grouping_by_outer_a_v2_skip_deltaij_zero_directK','linearity_identity':'sum_bc,k eps*C_a(K)_ij Xi_bc^ijk = C_a(K)_ij sum_bc,k eps*Xi_bc^ijk','source_node':source_v,'input_logical_basis_index':input_index,'input_K_labels':list(initial[1]),'outer_a':outer_a,'Jmax':FULL.JMAX2/2,'all_middle_provenance':{'run_id':middle_run_id,'head_sha':middle_head_sha,'summary_sha256':sha_json(summary)},'zero_prefix_indices':EXPECTED_ZERO,'nonzero_prefix_indices':EXPECTED_LIVE,'live_triples':[list(x) for x in live_triples],'zero_before_outer_triples':[list(x) for x in zero_triples],'prefix_audit':prefix_audit,'outer_channels':rows,'outer_CK_call_count':calls,'offdiagonal_direct_K_evaluations_skipped':offdiag_skips,'partial_support':len(total),'partial_norm':FULL.norm(total),'partial_scalar_diagnostics':scalar,'max_spin_reached':max_spin,'max_diagnostics':diag,'runtime_exact_cache':cache_info,'habitat_descriptor':hab,'habitat_hash':PLAN.canonical_hash(hab),'boundary_domain_hash':PLAN.boundary_domain_hash(basis),'convention_descriptor':conv,'convention_hash':PLAN.canonical_hash(conv),'hard_integrity_checks':hard,'state':PLAN.encode_state(total),'claim_boundary':'Exact signed contribution of one outer edge to the first raw Lorentzian boundary column. Global prefix/k regrouping and removal of the exactly zero off-diagonal direct-K branch change no operator. No symmetry assumption, physical projector or cosmology is inferred.'}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--middle-dir',type=Path,required=True);ap.add_argument('--outer-a',type=int,required=True);ap.add_argument('--source-node',type=int,default=0);ap.add_argument('--input-index',type=int,default=0);ap.add_argument('--middle-run-id',type=int);ap.add_argument('--middle-head-sha');ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();out=run(a.middle_dir,a.outer_a,a.source_node,a.input_index,a.middle_run_id,a.middle_head_sha);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in out.items() if k not in ('state','prefix_audit','outer_channels')},indent=2));return 0 if out['passed'] else 1
if __name__=='__main__': raise SystemExit(main())
