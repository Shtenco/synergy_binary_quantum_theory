#!/usr/bin/env python3
"""Compute one complete raw Lorentzian H_L boundary column in one exact DAG session.

For one K5 source node and one of the 32 frozen all-j=1/2 logical boundary
states this gate evaluates the preregistered

    H_L^raw ~ epsilon Tr_aux[C(K) C(K) C(V)]

without assuming the first-column zero/live prefix pattern.

Execution DAG:
  1. compute the 16 unique C_c(V)_{ki}|psi> states once;
  2. evaluate all 12 ordered middle prefixes C_b(K) C_c(V), sharing caches;
  3. classify each prefix pathwise zero/nonzero as a measured outcome;
  4. for every nonzero prefix, group k=0,1 exactly for fixed (i,j) before the
     final linear C_a(K)_{ij} hit, with no pre-outer tolerance pruning;
  5. assemble the two epsilon-signed triples per prefix and verify exact
     12-prefix / 24-triple coverage;
  6. emit the complete sparse covariant outgoing column plus J=0 Gauss/logical
     return diagnostics.

The operator ordering, scalar-channel rule, Jmax=7/2 and numerical conventions
are identical to the frozen first-column calculation.  This is a single-hit
boundary column only; it is not the Hermitian physical H_L choice, HH-safe
habitat, master projector or cosmology.
"""
from __future__ import annotations

import argparse, json
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_lorentzian_middle_prefix_gate as MID
import peter_weyl_lorentzian_pruned_prefix_worker as PLAN
import peter_weyl_zeroaware_volume_migration_experiment as ZVM


def sum_no_prune(*states):
    out={}
    for state in states:
        for key,z in state.items(): out[key]=out.get(key,0j)+complex(z)
    return {k:z for k,z in out.items() if z!=0j}


def scalar_channel(state):
    return {key:z for key,z in state.items() if int(key[2]) in (0,2)}


def run(source_v=0,input_index=0):
    ZVM.patch_and_clear()
    basis=FULL.RAW.PW.basis_full_jhalf()
    if len(basis)!=32: raise RuntimeError(f"expected 32 logical inputs, found {len(basis)}")
    if not (0<=input_index<len(basis)): raise ValueError("input_index outside logical basis")
    initial=basis[input_index]
    psi=FULL.RAW.CV.gauss_to_covariant({initial:1+0j},source_v)
    neighbors=tuple(FULL.RAW.PW.NEIG[source_v]); pairs=MID.ordered_pairs(source_v)
    if len(neighbors)!=4 or len(pairs)!=12: raise RuntimeError("expected four-valent source and 12 ordered prefixes")

    cv_states={};cv_leak={};prefix_rows=[];total={};all_triples=[];max_spin=0.0
    global_diag={
        "CV_complete_basis_leakage":0.0,
        "CK_outer_complete_basis_leakage":0.0,
        "CK_internal_volume_sector_leakage":0.0,
        "CK_complete_charge_basis_leakage":0.0,
    }
    max_term_nonscalar_fraction=0.0
    max_term_nonscalar_norm_if_near_zero=0.0
    middle_ck_calls=0;outer_ck_calls=0

    old,caches=FULL.install_sine_ordering()
    try:
        # 16 unique right-most volume commutator states.
        for c in neighbors:
            for i in range(2):
                for k in range(2):
                    s1,leak=FULL.RAW.COMP.C_volume_component(psi,source_v,c,k,i,FULL.JMAX2)
                    cv_states[(c,i,k)]=s1;cv_leak[(c,i,k)]=float(leak)
                    global_diag["CV_complete_basis_leakage"]=max(global_diag["CV_complete_basis_leakage"],float(leak))
                    max_spin=max(max_spin,FULL.max_spin(s1))

        for pair_index,(b,c) in enumerate(pairs):
            middle={};middle_paths=[];prefix_diag={
                "CV_complete_basis_leakage":max(cv_leak[(c,i,k)] for i in range(2) for k in range(2)),
                "CK_outer_complete_basis_leakage":0.0,
                "CK_internal_volume_sector_leakage":0.0,
                "CK_complete_charge_basis_leakage":0.0,
            }
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        s1=cv_states[(c,i,k)]
                        if s1:
                            s2,d2=FULL.RAW.KCOMP.C_K_component(s1,source_v,b,j,k,FULL.JMAX2)
                            middle_ck_calls+=1
                            FULL.update_diag(prefix_diag,d2);FULL.update_diag(global_diag,d2)
                        else:
                            s2={}
                        middle[(i,j,k)]=s2
                        max_spin=max(max_spin,FULL.max_spin(s2))
                        middle_paths.append({
                            "indices":[i,j,k],
                            "after_CV_support":len(s1),"after_CV_norm":FULL.norm(s1),
                            "after_middle_CK_support":len(s2),"after_middle_CK_norm":FULL.norm(s2),
                        })
            zero=all(not s for s in middle.values())
            plan=PLAN.epsilon_outer_terms(source_v,b,c)
            if len(plan)!=2: raise RuntimeError(f"prefix {(b,c)} does not map to exactly two outer triples")
            all_triples.extend(tuple(int(x) for x in row["ordered_edges"]) for row in plan)
            term_states=[{} for _ in plan]
            grouped_rows=[]

            if not zero:
                for i in range(2):
                    for j in range(2):
                        q0=scalar_channel(middle[(i,j,0)])
                        q1=scalar_channel(middle[(i,j,1)])
                        grouped=sum_no_prune(q0,q1)
                        max_spin=max(max_spin,FULL.max_spin(q0),FULL.max_spin(q1),FULL.max_spin(grouped))
                        outer=[]
                        for ti,row in enumerate(plan):
                            if grouped:
                                s3,d3=FULL.RAW.KCOMP.C_K_component(grouped,source_v,row["a"],i,j,FULL.JMAX2)
                                outer_ck_calls+=1
                                FULL.update_diag(prefix_diag,d3);FULL.update_diag(global_diag,d3)
                            else:
                                s3={}
                            FULL.add(term_states[ti],s3)
                            max_spin=max(max_spin,FULL.max_spin(s3))
                            outer.append({"a":int(row["a"]),"support":len(s3),"norm":FULL.norm(s3)})
                        grouped_rows.append({
                            "indices_ij":[i,j],"contributing_k":[0,1],
                            "k0_scalar_support":len(q0),"k0_scalar_norm":FULL.norm(q0),
                            "k1_scalar_support":len(q1),"k1_scalar_norm":FULL.norm(q1),
                            "grouped_scalar_support":len(grouped),"grouped_scalar_norm":FULL.norm(grouped),
                            "outer_actions":outer,
                        })

            prefix_partial={};term_rows=[]
            for row,state in zip(plan,term_states):
                sd=FULL.scalar_diagnostics(state)
                if sd["norm"]>FULL.NONZERO_TOL:
                    max_term_nonscalar_fraction=max(max_term_nonscalar_fraction,sd["nonscalar_weight_fraction"])
                else:
                    max_term_nonscalar_norm_if_near_zero=max(max_term_nonscalar_norm_if_near_zero,sd["nonscalar_norm"])
                FULL.add(prefix_partial,state,scale=row["sign"])
                term_rows.append({
                    **row,"support":len(state),"norm":FULL.norm(state),
                    "scalar_diagnostics":sd,"max_spin_reached":FULL.max_spin(state),
                })
            FULL.add(total,prefix_partial)
            prefix_scalar=FULL.scalar_diagnostics(prefix_partial)
            prefix_rows.append({
                "pair_index":int(pair_index),"b":int(b),"c":int(c),
                "middle_science_status":"MIDDLE_PREFIX_ZERO_PATHWISE" if zero else "MIDDLE_PREFIX_NONZERO",
                "prefix_zero_pathwise":bool(zero),
                "middle_paths":middle_paths,"grouped_outer_channels":grouped_rows,
                "outer_terms":term_rows,
                "signed_prefix_support":len(prefix_partial),"signed_prefix_norm":FULL.norm(prefix_partial),
                "signed_prefix_scalar_diagnostics":prefix_scalar,"max_diagnostics":prefix_diag,
            })
        cache_info={name:{"hits":fun.cache_info().hits,"misses":fun.cache_info().misses,"currsize":fun.cache_info().currsize} for name,fun in caches.items()}
    finally:
        FULL.restore_ordering(old)

    zero_indices=[r["pair_index"] for r in prefix_rows if r["prefix_zero_pathwise"]]
    live_indices=[r["pair_index"] for r in prefix_rows if not r["prefix_zero_pathwise"]]
    scalar=FULL.scalar_diagnostics(total)
    gauss,mapdiag=FULL.project_covariant_J0_to_gauss(total,source_v)
    logical=FULL.logical_projection(gauss);logical_norm=FULL.norm(logical);initial_amp=logical.get(initial,0j)
    logical_rows=[]
    for idx,key in enumerate(basis):
        amp=logical.get(key,0j)
        if abs(amp)>FULL.TOL:
            logical_rows.append({"logical_basis_index":idx,"K_labels":list(key[1]),"amp":[float(amp.real),float(amp.imag)],"abs":abs(amp)})

    expected_triples=set()
    for b,c in pairs:
        expected_triples.update(tuple(int(x) for x in row["ordered_edges"]) for row in PLAN.epsilon_outer_terms(source_v,b,c))
    hard={
        "sixteen_unique_CV_states":len(cv_states)==16,
        "twelve_prefixes_once":len(prefix_rows)==12 and {r["pair_index"] for r in prefix_rows}==set(range(12)),
        "zero_live_partition_complete":len(zero_indices)+len(live_indices)==12 and not(set(zero_indices)&set(live_indices)),
        "all_24_ordered_triples_unique":len(all_triples)==24 and len(set(all_triples))==24 and set(all_triples)==expected_triples,
        "outer_call_count_matches_live_grouped_DAG":outer_ck_calls==8*len(live_indices),
        "CV_complete_basis_leakage_below_1e-9":global_diag["CV_complete_basis_leakage"]<1e-9,
        "CK_outer_complete_basis_leakage_below_1e-9":global_diag["CK_outer_complete_basis_leakage"]<1e-9,
        "CK_internal_volume_sector_leakage_below_1e-9":global_diag["CK_internal_volume_sector_leakage"]<1e-9,
        "outer_terms_scalar_within_frozen_threshold":max_term_nonscalar_fraction<1e-8 and max_term_nonscalar_norm_if_near_zero<FULL.NONZERO_TOL,
        "full_signed_output_scalar_within_frozen_threshold":FULL.scalar_ok(scalar),
        "spin_cutoff_respected":max_spin<=FULL.JMAX2/2+1e-12,
        "J0_reverse_projection_has_no_invalid_keys":not mapdiag["invalid_J0_covariant_keys"],
        "J0_reverse_projection_has_no_collisions":int(mapdiag["mapping_collisions"])==0,
    }
    if not total: science="FULL_RAW_HL_COLUMN_ZERO"
    elif logical_norm>FULL.NONZERO_TOL: science="FULL_RAW_HL_COLUMN_NONZERO_WITH_LOGICAL_RETURN"
    else: science="FULL_RAW_HL_COLUMN_NONZERO_LOGICAL_RETURN_ZERO"
    conv=PLAN.convention_descriptor(source_v);hab=PLAN.habitat_descriptor(source_v)
    return {
        "schema":"BQG_LORENTZIAN_FULL_COLUMN_DAG_V1",
        "passed":bool(all(hard.values())),"science_status":science,
        "execution_mode":"single_session_16CV_12prefix_grouped_k_outer_v1",
        "source_node":int(source_v),"input_logical_basis_index":int(input_index),"input_K_labels":list(initial[1]),
        "Jmax":FULL.JMAX2/2,
        "zero_prefix_indices":zero_indices,"nonzero_prefix_indices":live_indices,
        "zero_prefix_count":len(zero_indices),"nonzero_prefix_count":len(live_indices),
        "unique_CV_state_count":len(cv_states),"middle_CK_call_count":middle_ck_calls,"outer_CK_call_count":outer_ck_calls,
        "prefixes":prefix_rows,
        "full_outgoing_support":len(total),"full_outgoing_norm":FULL.norm(total),"full_scalar_diagnostics":scalar,"max_spin_reached":max_spin,
        "gauss_reverse_projection":{"support":len(gauss),"norm":FULL.norm(gauss),"diagnostics":mapdiag},
        "logical_return":{"support":len(logical),"norm":logical_norm,"fraction_of_full_norm":logical_norm/max(FULL.norm(total),1e-300),
                          "initial_return_amplitude":[float(initial_amp.real),float(initial_amp.imag)],"nonzero_amplitudes":logical_rows},
        "max_diagnostics":global_diag,"hard_integrity_checks":hard,"runtime_exact_cache":cache_info,
        "habitat_descriptor":hab,"habitat_hash":PLAN.canonical_hash(hab),
        "boundary_domain_hash":PLAN.boundary_domain_hash(basis),"convention_descriptor":conv,"convention_hash":PLAN.canonical_hash(conv),
        "state":PLAN.encode_state(total),
        "claim_boundary":"Complete raw single-hit Lorentzian boundary column at Jmax=7/2 from one exact cache-preserving DAG. It is not yet the Hermitian physical H_L convention, not HH-safe Jmax=13/2 closure data, and does not emit P_phys or cosmology."
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-node',type=int,default=0);ap.add_argument('--input-index',type=int,default=0);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();out=run(a.source_node,a.input_index);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('state','prefixes')},indent=2));return 0 if out['passed'] else 1

if __name__=='__main__': raise SystemExit(main())
