#!/usr/bin/env python3
"""Continue one exact Lorentzian prefix from a serialized middle-state artifact.

The expensive states

    Xi_bc^{ijk} = C_b(K)_{jk} C_c(V)_{ki} |psi>

are consumed from a completed exact middle-prefix run.  This worker evaluates
only the two epsilon-ordered outer C_a(K) continuations sharing the same (b,c).
No operator order, cutoff, scalar-channel rule, or coefficient is changed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_lorentzian_pruned_prefix_worker as PLAN
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

EXPECTED_DEFINITION="all eight auxiliary states C_b(K) C_c(V)|psi> before the outer C_a(K)"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decode_state(rows):
    out={}
    for r in rows:
        key=(
            tuple(int(x) for x in r["spins"]),
            tuple(int(x) for x in r["Kother"]),
            int(r["J2"]),int(r["M2"]),int(r["K12"]),int(r["K34"]),
        )
        z=complex(float(r["amp"][0]),float(r["amp"][1]))
        out[key]=out.get(key,0j)+z
    return out


def load_middle(path: Path, source_v: int, input_index: int):
    data=json.loads(path.read_text(encoding="utf-8"))
    pair=data.get("ordered_pair",{})
    b=int(pair.get("b",-1)); c=int(pair.get("c",-1))
    if (b,c) not in PLAN.CANDIDATE_PREFIXES:
        raise ValueError(f"middle artifact prefix {(b,c)} is not in frozen surviving set")
    checks={
        "middle_artifact_passed":bool(data.get("passed",False)),
        "middle_artifact_reports_nonzero":data.get("science_status")=="MIDDLE_PREFIX_NONZERO",
        "source_node_matches":int(data.get("source_node",-1))==int(source_v),
        "input_index_matches":int(data.get("input_logical_basis_index",-1))==int(input_index),
        "Jmax_matches_frozen_gate":abs(float(data.get("Jmax",-1))-FULL.JMAX2/2)<1e-12,
        "definition_matches":data.get("definition")==EXPECTED_DEFINITION,
        "middle_not_pathwise_zero":data.get("prefix_zero_pathwise") is False,
        "middle_integrity_checks_pass":all(bool(x) for x in data.get("hard_integrity_checks",{}).values()),
    }
    paths=data.get("paths",[])
    expected={(i,j,k) for i in range(2) for j in range(2) for k in range(2)}
    seen=set(); states={}; path_meta=[]
    for row in paths:
        idx=tuple(int(x) for x in row.get("indices",[]))
        if len(idx)!=3 or idx in seen:
            raise ValueError("duplicate or malformed auxiliary path in middle artifact")
        seen.add(idx)
        state=decode_state(row.get("middle_state",[]))
        states[idx]=state
        path_meta.append({
            "indices":list(idx),
            "after_CV_support":int(row.get("after_CV_support",0)),
            "after_CV_norm":float(row.get("after_CV_norm",0.0)),
            "after_middle_CK_support":int(row.get("after_middle_CK_support",0)),
            "after_middle_CK_norm":float(row.get("after_middle_CK_norm",0.0)),
            "serialized_middle_support":len(state),
            "serialized_middle_norm":FULL.norm(state),
        })
    checks["exactly_eight_unique_aux_paths"]=(seen==expected)
    checks["every_serialized_middle_path_nonzero"]=all(bool(states[x]) for x in expected)
    # Verify serialized states reproduce the recorded path norms/supports.
    tol=5e-12
    checks["serialized_support_matches_metadata"]=all(
        r["after_middle_CK_support"]==r["serialized_middle_support"] for r in path_meta
    )
    checks["serialized_norm_matches_metadata"]=all(
        abs(r["after_middle_CK_norm"]-r["serialized_middle_norm"])<=tol*max(1.0,r["after_middle_CK_norm"])
        for r in path_meta
    )
    if not all(checks.values()):
        raise RuntimeError(f"middle artifact failed provenance/integrity checks: {checks}")
    return data,b,c,states,path_meta,checks


def run(middle_path: Path, source_v=0, input_index=0):
    ZVM.patch_and_clear()
    middle,b,c,middle_states,path_meta,source_checks=load_middle(middle_path,source_v,input_index)
    basis=FULL.RAW.PW.basis_full_jhalf()
    if len(basis)!=32 or not (0<=input_index<len(basis)):
        raise RuntimeError("frozen 32D all-j=1/2 boundary basis unavailable")
    initial=basis[input_index]
    plan=PLAN.epsilon_outer_terms(source_v,b,c)
    term_states=[{} for _ in plan]
    diagmax={
        "CK_outer_complete_basis_leakage":0.0,
        "CK_internal_volume_sector_leakage":0.0,
        "CK_complete_charge_basis_leakage":0.0,
    }
    path_rows=[]; max_spin=0.0

    old,caches=FULL.install_sine_ordering()
    try:
        for idx in sorted(middle_states):
            i,j,k=idx; s2=middle_states[idx]
            # Frozen exact scalar-channel selection from the preregistered full gate.
            s2_scalar={key:amp for key,amp in s2.items() if key[2] in (0,2)}
            outer=[]
            for ti,row in enumerate(plan):
                if s2_scalar:
                    s3,d3=FULL.RAW.KCOMP.C_K_component(s2_scalar,source_v,row["a"],i,j,FULL.JMAX2)
                    FULL.update_diag(diagmax,d3)
                else:
                    s3={}
                FULL.add(term_states[ti],s3)
                outer.append({"a":row["a"],"support":len(s3),"norm":FULL.norm(s3)})
                max_spin=max(max_spin,FULL.max_spin(s3))
            max_spin=max(max_spin,FULL.max_spin(s2),FULL.max_spin(s2_scalar))
            path_rows.append({
                "indices":list(idx),
                "loaded_middle_support":len(s2),"loaded_middle_norm":FULL.norm(s2),
                "scalar_relevant_middle_support":len(s2_scalar),
                "outer_actions":outer,
            })
        cache_info={name:{"hits":fun.cache_info().hits,"misses":fun.cache_info().misses,"currsize":fun.cache_info().currsize} for name,fun in caches.items()}
    finally:
        FULL.restore_ordering(old)

    partial={}; term_rows=[]
    max_term_nonscalar_fraction=0.0; max_term_nonscalar_norm_if_near_zero=0.0
    for row,state in zip(plan,term_states):
        sd=FULL.scalar_diagnostics(state)
        if sd["norm"]>FULL.NONZERO_TOL:
            max_term_nonscalar_fraction=max(max_term_nonscalar_fraction,sd["nonscalar_weight_fraction"])
        else:
            max_term_nonscalar_norm_if_near_zero=max(max_term_nonscalar_norm_if_near_zero,sd["nonscalar_norm"])
        FULL.add(partial,state,scale=row["sign"])
        term_rows.append({**row,"support":len(state),"norm":FULL.norm(state),"scalar_diagnostics":sd,"max_spin_reached":FULL.max_spin(state)})

    partial_scalar=FULL.scalar_diagnostics(partial)
    middle_max=middle.get("max_diagnostics",{})
    hard={
        **source_checks,
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
    conv=PLAN.convention_descriptor(source_v); hab=PLAN.habitat_descriptor(source_v)
    return {
        "schema":"BQG_LORENTZIAN_OUTER_FROM_MIDDLE_V1",
        "passed":bool(all(hard.values())),
        "science_status":"PREFIX_ZERO_AFTER_OUTER_EPSILON_SUM" if not partial else "PREFIX_NONZERO_OUTGOING_PARTIAL",
        "source_node":int(source_v),"input_logical_basis_index":int(input_index),"input_K_labels":list(initial[1]),
        "prefix":{"b":b,"c":c,"source_pair_index":int(middle["ordered_pair"].get("pair_index",-1))},
        "Jmax":FULL.JMAX2/2,
        "middle_artifact":{"path":str(middle_path),"sha256":sha256_file(middle_path),"science_status":middle.get("science_status"),"runtime_exact_cache":middle.get("runtime_exact_cache")},
        "middle_path_metadata":path_meta,
        "outer_terms":term_rows,"outer_path_diagnostics":path_rows,
        "signed_partial_support":len(partial),"signed_partial_norm":FULL.norm(partial),
        "signed_partial_scalar_diagnostics":partial_scalar,"max_spin_reached":max_spin,
        "outer_max_diagnostics":diagmax,"hard_integrity_checks":hard,"runtime_exact_cache":cache_info,
        "habitat_descriptor":hab,"habitat_hash":PLAN.canonical_hash(hab),
        "boundary_domain_hash":PLAN.boundary_domain_hash(basis),
        "convention_descriptor":conv,"convention_hash":PLAN.canonical_hash(conv),
        "state":PLAN.encode_state(partial),
        "claim_boundary":"Exact signed two-triple continuation of one previously serialized nonzero middle prefix. Six such workers plus the six independently certified pathwise-zero prefixes are required for the first raw H_L column.",
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--middle-prefix-json",type=Path,required=True)
    ap.add_argument("--source-node",type=int,default=0);ap.add_argument("--input-index",type=int,default=0)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args(); out=run(a.middle_prefix_json,a.source_node,a.input_index)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in out.items() if k not in ("state","middle_path_metadata","outer_path_diagnostics")},indent=2))
    return 0 if out["passed"] else 1

if __name__=="__main__": raise SystemExit(main())
