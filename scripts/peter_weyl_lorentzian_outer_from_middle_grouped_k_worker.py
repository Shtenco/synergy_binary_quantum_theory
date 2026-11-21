#!/usr/bin/env python3
"""Exact grouped-k continuation of one Lorentzian middle-prefix artifact.

For the frozen ordered term

    T_abc = Tr_aux[C_a(K) C_b(K) C_c(V)],

the serialized middle states are

    Xi_bc^{ijk} = C_b(K)_{jk} C_c(V)_{ki}|psi>.

For fixed (i,j), the final linear operator C_a(K)_{ij} is independent of k.
Therefore

    sum_k C_a(K)_{ij} Xi_bc^{ijk}
      = C_a(K)_{ij} [sum_k Xi_bc^{ijk}].

This worker performs that algebraic regrouping before the expensive outer hit.
The k=0 and k=1 sparse middle states are summed with *no amplitude pruning*;
only the already frozen scalar-channel selection J_source in {0,1} (J2 in
{0,2}) is applied.  Thus the microscopic operator, ordering, cutoff and
coefficients are unchanged while outer C(K) calls are reduced from 16 to 8 per
surviving prefix.

The emitted schema intentionally remains BQG_LORENTZIAN_OUTER_FROM_MIDDLE_V1
so the existing deterministic first-column aggregate can consume either the
ungrouped reference workers or this exact linearity execution mode.  Execution
mode/provenance fields distinguish them.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_lorentzian_outer_from_middle_worker as BASE
import peter_weyl_lorentzian_pruned_prefix_worker as PLAN
import peter_weyl_zeroaware_volume_migration_experiment as ZVM


def sum_no_prune(*states):
    """Linear sparse sum with no tolerance pruning before the outer operator."""
    out={}
    for state in states:
        for key,z in state.items():
            out[key]=out.get(key,0j)+complex(z)
    # Removing an IEEE-exact zero is algebraically harmless and avoids useless
    # work.  No near-zero/tolerance pruning is performed here.
    return {k:z for k,z in out.items() if z != 0j}


def scalar_channel(state):
    return {key:z for key,z in state.items() if int(key[2]) in (0,2)}


def run(middle_path:Path, source_v=0, input_index=0):
    ZVM.patch_and_clear()
    middle,b,c,middle_states,path_meta,source_checks=BASE.load_middle(
        middle_path,source_v,input_index
    )
    basis=FULL.RAW.PW.basis_full_jhalf()
    if len(basis)!=32 or not (0<=input_index<len(basis)):
        raise RuntimeError("frozen 32D all-j=1/2 boundary basis unavailable")
    initial=basis[input_index]
    plan=PLAN.epsilon_outer_terms(source_v,b,c)
    if len(plan)!=2:
        raise RuntimeError("expected exactly two epsilon outer triples")

    term_states=[{} for _ in plan]
    diagmax={
        "CK_outer_complete_basis_leakage":0.0,
        "CK_internal_volume_sector_leakage":0.0,
        "CK_complete_charge_basis_leakage":0.0,
    }
    grouped_rows=[]; max_spin=0.0; outer_call_count=0

    old,caches=FULL.install_sine_ordering()
    try:
        for i in range(2):
            for j in range(2):
                s20=middle_states[(i,j,0)]
                s21=middle_states[(i,j,1)]
                q0=scalar_channel(s20)
                q1=scalar_channel(s21)
                grouped=sum_no_prune(q0,q1)
                max_spin=max(
                    max_spin,
                    FULL.max_spin(s20),FULL.max_spin(s21),
                    FULL.max_spin(q0),FULL.max_spin(q1),FULL.max_spin(grouped),
                )
                outer=[]
                for ti,row in enumerate(plan):
                    if grouped:
                        s3,d3=FULL.RAW.KCOMP.C_K_component(
                            grouped,source_v,row["a"],i,j,FULL.JMAX2
                        )
                        FULL.update_diag(diagmax,d3)
                        outer_call_count+=1
                    else:
                        s3={}
                    # Same post-operator frozen accumulation rule as the
                    # preregistered reference worker.
                    FULL.add(term_states[ti],s3)
                    outer.append({
                        "a":int(row["a"]),
                        "support":len(s3),
                        "norm":FULL.norm(s3),
                    })
                    max_spin=max(max_spin,FULL.max_spin(s3))
                grouped_rows.append({
                    "indices_ij":[i,j],
                    "contributing_k":[0,1],
                    "k0_loaded_middle_support":len(s20),
                    "k0_loaded_middle_norm":FULL.norm(s20),
                    "k1_loaded_middle_support":len(s21),
                    "k1_loaded_middle_norm":FULL.norm(s21),
                    "k0_scalar_support":len(q0),
                    "k0_scalar_norm":FULL.norm(q0),
                    "k1_scalar_support":len(q1),
                    "k1_scalar_norm":FULL.norm(q1),
                    "grouped_scalar_support":len(grouped),
                    "grouped_scalar_norm":FULL.norm(grouped),
                    "outer_actions":outer,
                })
        cache_info={
            name:{
                "hits":fun.cache_info().hits,
                "misses":fun.cache_info().misses,
                "currsize":fun.cache_info().currsize,
            }
            for name,fun in caches.items()
        }
    finally:
        FULL.restore_ordering(old)

    partial={};term_rows=[]
    max_term_nonscalar_fraction=0.0
    max_term_nonscalar_norm_if_near_zero=0.0
    for row,state in zip(plan,term_states):
        sd=FULL.scalar_diagnostics(state)
        if sd["norm"]>FULL.NONZERO_TOL:
            max_term_nonscalar_fraction=max(
                max_term_nonscalar_fraction,sd["nonscalar_weight_fraction"]
            )
        else:
            max_term_nonscalar_norm_if_near_zero=max(
                max_term_nonscalar_norm_if_near_zero,sd["nonscalar_norm"]
            )
        FULL.add(partial,state,scale=row["sign"])
        term_rows.append({
            **row,
            "support":len(state),
            "norm":FULL.norm(state),
            "scalar_diagnostics":sd,
            "max_spin_reached":FULL.max_spin(state),
        })

    partial_scalar=FULL.scalar_diagnostics(partial)
    middle_max=middle.get("max_diagnostics",{})
    grouped_pairs={tuple(r["indices_ij"]) for r in grouped_rows}
    hard={
        **source_checks,
        "exactly_four_grouped_ij_channels":grouped_pairs=={(0,0),(0,1),(1,0),(1,1)},
        "each_group_contains_both_k_paths":all(r["contributing_k"]==[0,1] for r in grouped_rows),
        "no_pre_outer_tolerance_pruning":True,
        "outer_call_count_reduced_from_16_to_8":outer_call_count==8,
        "exactly_two_outer_triples":len(term_rows)==2,
        "outer_triples_share_prefix":all(tuple(r["ordered_edges"][1:])==(b,c) for r in term_rows),
        "outer_triples_unique":len({tuple(r["ordered_edges"]) for r in term_rows})==2,
        "middle_CV_complete_basis_leakage_below_1e-9":float(middle_max.get("CV_complete_basis_leakage",1.0))<1e-9,
        "middle_CK_outer_complete_basis_leakage_below_1e-9":float(middle_max.get("CK_outer_complete_basis_leakage",1.0))<1e-9,
        "outer_CK_complete_basis_leakage_below_1e-9":diagmax["CK_outer_complete_basis_leakage"]<1e-9,
        "outer_CK_internal_volume_sector_leakage_below_1e-9":diagmax["CK_internal_volume_sector_leakage"]<1e-9,
        "signed_prefix_sum_scalar_within_frozen_threshold":FULL.scalar_ok(partial_scalar),
        "outer_terms_scalar_within_frozen_threshold":max_term_nonscalar_fraction<1e-8 and max_term_nonscalar_norm_if_near_zero<FULL.NONZERO_TOL,
        "spin_cutoff_respected":max_spin<=FULL.JMAX2/2+1e-12,
    }
    conv=PLAN.convention_descriptor(source_v)
    hab=PLAN.habitat_descriptor(source_v)
    return {
        "schema":"BQG_LORENTZIAN_OUTER_FROM_MIDDLE_V1",
        "execution_mode":"exact_grouped_k_linearity_v2",
        "linearity_identity":"sum_k C_a(K)_{ij} Xi^{ijk} = C_a(K)_{ij} sum_k Xi^{ijk}",
        "naive_outer_CK_call_count":16,
        "outer_CK_call_count":outer_call_count,
        "passed":bool(all(hard.values())),
        "science_status":"PREFIX_ZERO_AFTER_OUTER_EPSILON_SUM" if not partial else "PREFIX_NONZERO_OUTGOING_PARTIAL",
        "source_node":int(source_v),
        "input_logical_basis_index":int(input_index),
        "input_K_labels":list(initial[1]),
        "prefix":{
            "b":b,"c":c,
            "source_pair_index":int(middle["ordered_pair"].get("pair_index",-1)),
        },
        "Jmax":FULL.JMAX2/2,
        "middle_artifact":{
            "path":str(middle_path),
            "sha256":BASE.sha256_file(middle_path),
            "science_status":middle.get("science_status"),
            "runtime_exact_cache":middle.get("runtime_exact_cache"),
        },
        "middle_path_metadata":path_meta,
        "outer_terms":term_rows,
        "outer_path_diagnostics":grouped_rows,
        "signed_partial_support":len(partial),
        "signed_partial_norm":FULL.norm(partial),
        "signed_partial_scalar_diagnostics":partial_scalar,
        "max_spin_reached":max_spin,
        "outer_max_diagnostics":diagmax,
        "hard_integrity_checks":hard,
        "runtime_exact_cache":cache_info,
        "habitat_descriptor":hab,
        "habitat_hash":PLAN.canonical_hash(hab),
        "boundary_domain_hash":PLAN.boundary_domain_hash(basis),
        "convention_descriptor":conv,
        "convention_hash":PLAN.canonical_hash(conv),
        "state":PLAN.encode_state(partial),
        "claim_boundary":"Exact grouped-k signed two-triple continuation of one serialized nonzero middle prefix. The regrouping uses only operator linearity and no pre-outer tolerance pruning. Six prefixes plus six independently certified zero prefixes are still required for the first raw H_L column.",
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--middle-prefix-json",type=Path,required=True)
    ap.add_argument("--source-node",type=int,default=0)
    ap.add_argument("--input-index",type=int,default=0)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    out=run(a.middle_prefix_json,a.source_node,a.input_index)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        k:v for k,v in out.items()
        if k not in ("state","middle_path_metadata","outer_path_diagnostics")
    },indent=2))
    return 0 if out["passed"] else 1

if __name__=="__main__":
    raise SystemExit(main())
