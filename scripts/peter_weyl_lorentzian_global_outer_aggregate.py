#!/usr/bin/env python3
"""Aggregate four globally grouped outer-edge packets into one candidate raw H_L column.

The four packets partition all 24 frozen epsilon ordered triples by outer edge.
Each packet absorbs the epsilon signs before its final C_a(K) actions.  This
aggregate performs only a sparse sum, then the existing J=0 -> Gauss reverse
map and logical-return diagnostics.

Important fail-closed boundary: global pre-outer regrouping is algebraically
linear, but the sparse implementation contains finite numerical pruning inside
low-level operator routines.  Therefore this accelerated execution is NOT
emitted under the production BQG_LORENTZIAN_FULL_COLUMN_DAG_V1 schema until an
independent frozen-reference sparse-equivalence certificate passes.  The
candidate microscopic result is still recorded in full, but cannot flow into
the 32-column Gram/master/projector producer yet.
"""
from __future__ import annotations

import argparse, glob, json, math
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_lorentzian_pruned_prefix_worker as PLAN


def decode(rows):
    out={}
    for r in rows:
        key=(tuple(int(x) for x in r['spins']),tuple(int(x) for x in r['Kother']),
             int(r['J2']),int(r['M2']),int(r['K12']),int(r['K34']))
        z=complex(float(r['amp'][0]),float(r['amp'][1]));out[key]=out.get(key,0j)+z
    return {k:z for k,z in out.items() if z!=0j}


def add_exact(dst,src):
    for k,z in src.items(): dst[k]=dst.get(k,0j)+z
    for k in [k for k,z in dst.items() if z==0j]: del dst[k]


def run(paths):
    packets=[json.loads(Path(p).read_text(encoding='utf-8')) for p in paths]
    if len(packets)!=4: raise RuntimeError(f'expected four outer-edge packets, got {len(packets)}')
    if any(p.get('schema')!='BQG_LORENTZIAN_GLOBAL_OUTER_A_V1' or p.get('passed') is not True for p in packets):
        raise RuntimeError('invalid outer-edge packet')
    packets.sort(key=lambda p:int(p['outer_a']))
    if [int(p['outer_a']) for p in packets]!=[1,2,3,4]: raise RuntimeError('outer edge coverage must be exactly 1..4')

    for field in ('source_node','input_logical_basis_index','Jmax','habitat_hash','boundary_domain_hash','convention_hash'):
        vals={json.dumps(p[field],sort_keys=True) for p in packets}
        if len(vals)!=1: raise RuntimeError(f'outer packets disagree on {field}')
    prov={(p['all_middle_provenance'].get('run_id'),p['all_middle_provenance'].get('head_sha'),p['all_middle_provenance'].get('summary_sha256')) for p in packets}
    if len(prov)!=1: raise RuntimeError('outer packets disagree on all-middle provenance')

    total={};all_live=[];all_zero=[];outer_rows=[];calls=0;max_spin=0.0
    for p in packets:
        s=decode(p.get('state',[]));add_exact(total,s)
        live=[tuple(x) for x in p.get('live_triples',[])];zero=[tuple(x) for x in p.get('zero_before_outer_triples',[])]
        all_live.extend(live);all_zero.extend(zero);calls+=int(p.get('outer_CK_call_count',0));max_spin=max(max_spin,float(p.get('max_spin_reached',0.0)))
        outer_rows.append({'outer_a':int(p['outer_a']),'science_status':p['science_status'],'support':len(s),'norm':FULL.norm(s),
                           'live_triple_count':len(live),'zero_before_outer_triple_count':len(zero),'outer_CK_call_count':int(p['outer_CK_call_count'])})

    accounted=all_live+all_zero
    if len(accounted)!=24 or len(set(accounted))!=24: raise RuntimeError('outer partition does not account for 24 unique epsilon triples')
    if len(all_live)!=12 or len(all_zero)!=12: raise RuntimeError('expected measured 12 live + 12 zero triple partition')

    source=int(packets[0]['source_node']);input_index=int(packets[0]['input_logical_basis_index']);basis=FULL.RAW.PW.basis_full_jhalf();initial=basis[input_index]
    scalar=FULL.scalar_diagnostics(total);gauss,mapdiag=FULL.project_covariant_J0_to_gauss(total,source);logical=FULL.logical_projection(gauss)
    logical_rows=[]
    for idx,key in enumerate(basis):
        z=logical.get(key,0j)
        if abs(z)>FULL.TOL: logical_rows.append({'logical_basis_index':idx,'K_labels':list(key[1]),'amp':[float(z.real),float(z.imag)],'abs':abs(z)})
    full_norm=FULL.norm(total);logical_norm=FULL.norm(logical);initial_amp=logical.get(initial,0j)
    cov_gauss_rel=abs(full_norm-FULL.norm(gauss))/max(full_norm,FULL.norm(gauss),1e-300)

    hard={
        'four_outer_edges_once':[int(p['outer_a']) for p in packets]==[1,2,3,4],
        'all_24_ordered_triples_unique':len(accounted)==24 and len(set(accounted))==24,
        'measured_12_live_12_zero_partition':len(all_live)==12 and len(all_zero)==12 and not(set(all_live)&set(all_zero)),
        'global_outer_CK_calls_at_most_16':calls<=16,
        'full_signed_output_scalar_within_frozen_threshold':FULL.scalar_ok(scalar),
        'spin_cutoff_respected':max_spin<=FULL.JMAX2/2+1e-12,
        'J0_reverse_projection_has_no_invalid_keys':not mapdiag.get('invalid_J0_covariant_keys'),
        'J0_reverse_projection_has_no_collisions':int(mapdiag.get('mapping_collisions',0))==0,
        'covariant_to_gauss_norm_isometry':cov_gauss_rel<2e-10,
    }
    if not total: science='GLOBAL_GROUPED_RAW_HL_COLUMN_ZERO'
    elif logical_norm>FULL.NONZERO_TOL: science='GLOBAL_GROUPED_RAW_HL_COLUMN_NONZERO_WITH_LOGICAL_RETURN'
    else: science='GLOBAL_GROUPED_RAW_HL_COLUMN_NONZERO_LOGICAL_RETURN_ZERO'
    return {
        'schema':'BQG_LORENTZIAN_GLOBAL_FIRST_COLUMN_CANDIDATE_V1',
        'passed':bool(all(hard.values())),
        'science_status':science,
        'execution_mode':'certified_all_middle_global_outer_a_linearity_v2',
        'production_equivalence_certified':False,
        'physical_projector_eligible':False,
        'reference_schema_after_equivalence':'BQG_LORENTZIAN_FULL_COLUMN_DAG_V1',
        'promotion_requirement':'Independent full-state sparse-equivalence against the frozen reference evaluator at the same source/input/Jmax/habitat/domain/convention hashes.',
        'source_node':source,'input_logical_basis_index':input_index,'input_K_labels':list(initial[1]),'Jmax':FULL.JMAX2/2,
        'zero_prefix_indices':[2,5,8,9,10,11],'nonzero_prefix_indices':[0,1,3,4,6,7],'zero_prefix_count':6,'nonzero_prefix_count':6,
        'unique_CV_state_count':16,'middle_CK_call_count':72,'outer_CK_call_count':calls,
        'outer_edge_partials':outer_rows,'live_ordered_triples':[list(x) for x in all_live],'zero_before_outer_ordered_triples':[list(x) for x in all_zero],
        'all_middle_provenance':{'run_id':next(iter(prov))[0],'head_sha':next(iter(prov))[1],'summary_sha256':next(iter(prov))[2]},
        'full_outgoing_support':len(total),'full_outgoing_norm':full_norm,'full_scalar_diagnostics':scalar,'max_spin_reached':max_spin,
        'gauss_reverse_projection':{'support':len(gauss),'norm':FULL.norm(gauss),'relative_norm_error':cov_gauss_rel,'diagnostics':mapdiag},
        'logical_return':{'support':len(logical),'norm':logical_norm,'fraction_of_full_norm':logical_norm/max(full_norm,1e-300),
                          'initial_return_amplitude':[float(initial_amp.real),float(initial_amp.imag)],'nonzero_amplitudes':logical_rows},
        'hard_integrity_checks':hard,
        'habitat_descriptor':packets[0]['habitat_descriptor'],'habitat_hash':packets[0]['habitat_hash'],
        'boundary_domain_hash':packets[0]['boundary_domain_hash'],'convention_descriptor':packets[0]['convention_descriptor'],'convention_hash':packets[0]['convention_hash'],
        'state':PLAN.encode_state(total),
        'claim_boundary':'Complete candidate first raw Lorentzian boundary column from global linear regrouping of the certified all-middle packet. Microscopic integrity is audited, but numerical equivalence to the frozen reference execution is not yet certified; therefore this packet cannot enter the 32-column Gram, HDA master, P_phys or cosmology producer.'
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--packet',type=Path,action='append');ap.add_argument('--packet-glob');ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    paths=list(a.packet or [])
    if a.packet_glob: paths.extend(Path(x) for x in glob.glob(a.packet_glob))
    out=run(paths);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='state'},indent=2));return 0 if out['passed'] else 1
if __name__=='__main__': raise SystemExit(main())
