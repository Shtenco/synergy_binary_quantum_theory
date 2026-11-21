#!/usr/bin/env python3
"""Aggregate one global direct-K packet and four -h_a K h_a^-1 packets.

This reconstructs the same globally regrouped candidate raw Lorentzian first
column while factoring the outer-edge-independent diagonal K term out of the
four edge workers.  The result remains fail-closed for production until an
independent frozen-reference sparse-equivalence certificate passes.
"""
from __future__ import annotations

import argparse, glob, json
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_lorentzian_pruned_prefix_worker as PLAN


def decode(rows):
    out={}
    for r in rows:
        key=(tuple(int(x) for x in r['spins']),tuple(int(x) for x in r['Kother']),int(r['J2']),int(r['M2']),int(r['K12']),int(r['K34']))
        out[key]=out.get(key,0j)+complex(float(r['amp'][0]),float(r['amp'][1]))
    return out


def add_exact(dst,src):
    for k,z in src.items(): dst[k]=dst.get(k,0j)+z
    for k in [k for k,z in dst.items() if z==0j]: del dst[k]


def run(direct_path:Path,hkh_paths):
    direct=json.loads(Path(direct_path).read_text(encoding='utf-8'))
    if direct.get('schema')!='BQG_LORENTZIAN_GLOBAL_DIRECT_K_V1' or direct.get('passed') is not True: raise RuntimeError('invalid global direct-K packet')
    hkh=[json.loads(Path(p).read_text(encoding='utf-8')) for p in hkh_paths]
    if len(hkh)!=4 or any(p.get('schema')!='BQG_LORENTZIAN_GLOBAL_HKH_OUTER_A_V1' or p.get('passed') is not True for p in hkh): raise RuntimeError('expected four valid hKh packets')
    hkh.sort(key=lambda p:int(p['outer_a']))
    if [int(p['outer_a']) for p in hkh]!=[1,2,3,4]: raise RuntimeError('hKh outer-edge coverage must be 1..4')

    packets=[direct]+hkh
    for field in ('source_node','input_logical_basis_index','Jmax','habitat_hash','boundary_domain_hash','convention_hash'):
        vals={json.dumps(p[field],sort_keys=True) for p in packets}
        if len(vals)!=1: raise RuntimeError(f'factorized packets disagree on {field}')
    prov={(p['all_middle_provenance'].get('run_id'),p['all_middle_provenance'].get('head_sha'),p['all_middle_provenance'].get('summary_sha256')) for p in packets}
    if len(prov)!=1: raise RuntimeError('factorized packets disagree on all-middle provenance')

    total={};add_exact(total,decode(direct['state']));triples=[];hrows=[]
    for p in hkh:
        s=decode(p['state']);add_exact(total,s);tt=[tuple(int(x) for x in q) for q in p['accounted_ordered_triples']];triples.extend(tt)
        hrows.append({'outer_a':int(p['outer_a']),'support':len(s),'norm':FULL.norm(s),'hKh_call_count':int(p['hKh_call_count'])})
    if len(triples)!=24 or len(set(triples))!=24: raise RuntimeError('hKh packets do not account for 24 unique ordered triples')

    source=int(direct['source_node']);input_index=int(direct['input_logical_basis_index']);basis=FULL.RAW.PW.basis_full_jhalf();initial=basis[input_index]
    scalar=FULL.scalar_diagnostics(total);gauss,mapdiag=FULL.project_covariant_J0_to_gauss(total,source);logical=FULL.logical_projection(gauss);logical_norm=FULL.norm(logical);full_norm=FULL.norm(total);initial_amp=logical.get(initial,0j)
    logical_rows=[]
    for idx,key in enumerate(basis):
        z=logical.get(key,0j)
        if abs(z)>FULL.TOL: logical_rows.append({'logical_basis_index':idx,'K_labels':list(key[1]),'amp':[float(z.real),float(z.imag)],'abs':abs(z)})
    cov_gauss_rel=abs(full_norm-FULL.norm(gauss))/max(full_norm,FULL.norm(gauss),1e-300)
    hard={'four_hKh_outer_edges_once':[int(p['outer_a']) for p in hkh]==[1,2,3,4],'all_24_ordered_triples_unique':len(triples)==24 and len(set(triples))==24,
          'global_direct_K_calls_at_most_two':int(direct['direct_K_call_count'])<=2,'total_hKh_calls_at_most_16':sum(int(p['hKh_call_count']) for p in hkh)<=16,
          'full_signed_output_scalar_within_frozen_threshold':FULL.scalar_ok(scalar),'J0_reverse_projection_has_no_invalid_keys':not mapdiag.get('invalid_J0_covariant_keys'),
          'J0_reverse_projection_has_no_collisions':int(mapdiag.get('mapping_collisions',0))==0,'covariant_to_gauss_norm_isometry':cov_gauss_rel<2e-10}
    if not total: science='GLOBAL_FACTORIZED_RAW_HL_COLUMN_ZERO'
    elif logical_norm>FULL.NONZERO_TOL: science='GLOBAL_FACTORIZED_RAW_HL_COLUMN_NONZERO_WITH_LOGICAL_RETURN'
    else: science='GLOBAL_FACTORIZED_RAW_HL_COLUMN_NONZERO_LOGICAL_RETURN_ZERO'
    return {'schema':'BQG_LORENTZIAN_GLOBAL_FACTORIZED_FIRST_COLUMN_CANDIDATE_V1','passed':bool(all(hard.values())),'science_status':science,
            'execution_mode':'global_direct_K_plus_four_outer_hKh_v1','production_equivalence_certified':False,'physical_projector_eligible':False,
            'source_node':source,'input_logical_basis_index':input_index,'input_K_labels':list(initial[1]),'Jmax':float(direct['Jmax']),
            'direct_K':{'support':int(direct['direct_support']),'norm':float(direct['direct_norm']),'call_count':int(direct['direct_K_call_count']),'diagonal_channels':direct['diagonal_channels']},
            'hKh_outer_edges':hrows,'total_operator_calls':int(direct['direct_K_call_count'])+sum(int(p['hKh_call_count']) for p in hkh),
            'full_outgoing_support':len(total),'full_outgoing_norm':full_norm,'full_scalar_diagnostics':scalar,
            'gauss_reverse_projection':{'support':len(gauss),'norm':FULL.norm(gauss),'relative_norm_error':cov_gauss_rel,'diagnostics':mapdiag},
            'logical_return':{'support':len(logical),'norm':logical_norm,'fraction_of_full_norm':logical_norm/max(full_norm,1e-300),'initial_return_amplitude':[float(initial_amp.real),float(initial_amp.imag)],'nonzero_amplitudes':logical_rows},
            'all_middle_provenance':{'run_id':next(iter(prov))[0],'head_sha':next(iter(prov))[1],'summary_sha256':next(iter(prov))[2]},
            'habitat_hash':direct['habitat_hash'],'boundary_domain_hash':direct['boundary_domain_hash'],'convention_hash':direct['convention_hash'],'hard_integrity_checks':hard,'state':PLAN.encode_state(total),
            'promotion_requirement':'Independent full-state sparse-equivalence against the frozen non-factorized reference evaluator at identical hashes.',
            'claim_boundary':'Complete factorized candidate first raw Lorentzian boundary column. The direct K term is actually evaluated globally (not assumed zero), and all four edge-dependent -h_a K h_a^-1 pieces are included. Production equivalence is still open, so this cannot enter the 32-column Gram, HDA master or P_phys yet.'}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--direct',type=Path,required=True);ap.add_argument('--hkh',type=Path,action='append');ap.add_argument('--hkh-glob');ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();paths=list(a.hkh or [])
    if a.hkh_glob: paths.extend(Path(x) for x in glob.glob(a.hkh_glob))
    out=run(a.direct,paths);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in out.items() if k!='state'},indent=2));return 0 if out['passed'] else 1
if __name__=='__main__': raise SystemExit(main())
