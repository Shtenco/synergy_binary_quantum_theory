#!/usr/bin/env python3
"""Exact DAG worker for one non-pruned Lorentzian epsilon middle prefix.

For a frozen ordered pair (b,c), compute the shared middle state

    Xi_bc(i,j,k) = C_b(K) C_c(V) |psi>

exactly once for each of the eight auxiliary paths.  The same Xi_bc is then fed
to the two and only two epsilon-ordered outer actions C_a(K) compatible with the
four-valent face assembly.  This is algebraically identical to the frozen
24-term operator, but avoids recomputing the expensive C(V)->C(K) prefix for the
two triples that share (b,c).

No approximation, amplitude pruning beyond the frozen base gate tolerance,
representation reduction, fitted beta, or altered operator ordering is used.
The emitted sparse `state` is the signed covariant outgoing partial column for
this prefix and can be summed by the deterministic aggregate gate.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

CANDIDATE_PREFIXES=((1,2),(2,1),(1,3),(3,1),(2,3),(3,2))


def canonical_hash(obj):
    raw=json.dumps(obj,sort_keys=True,separators=(",",":"),ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def encode_state(state):
    rows=[]
    for key,amp in sorted(state.items(),key=lambda kv:repr(kv[0])):
        spins,Kother,J2,M2,K12,K34=key
        rows.append({
            "spins":list(spins),"Kother":list(Kother),"J2":int(J2),"M2":int(M2),
            "K12":int(K12),"K34":int(K34),"amp":[float(amp.real),float(amp.imag)]
        })
    return rows


def boundary_domain_hash(basis):
    serial=[{"spins":list(key[0]),"K_labels":list(key[1])} for key in basis]
    return canonical_hash({"logical_basis":serial,"dimension":len(serial)})


def convention_descriptor(source_v):
    return {
        "graph":"K5",
        "source_node":int(source_v),
        "Jmax2":int(FULL.JMAX2),
        "euclidean_ordering":"H_E^sine=(T-T^dagger)/(2i)",
        "K_definition":"K=[V,H_E^sine]",
        "lorentzian_raw":"sum_r (-1)^r sum_perm sgn(perm) Tr_aux[C_a(K) C_b(K) C_c(V)]",
        "auxiliary_paths":"i,j,k in {0,1}",
        "scalar_channel_pruning":"before final C_a(K), keep source J2 in {0,2}; rank-(0+1) C(K) cannot return source J2=4 to J=0",
        "numerical_add_tolerance":float(FULL.TOL),
        "nonzero_reporting_tolerance":float(FULL.NONZERO_TOL),
    }


def habitat_descriptor(source_v):
    return {
        "graph":"K5",
        "representation":"finite Peter-Weyl graph-changing covariant habitat",
        "source_node":int(source_v),
        "boundary_sector":"all-j=1/2 gauge-invariant logical basis",
        "transient_Jmax2":int(FULL.JMAX2),
        "purpose":"single raw H_L hit; not HH-safe closure habitat",
    }


def epsilon_outer_terms(source_v,b,c):
    neighbors=tuple(FULL.RAW.PW.NEIG[source_v])
    if len(neighbors)!=4:
        raise RuntimeError("four-valent source required")
    if b==c or b not in neighbors or c not in neighbors:
        raise ValueError("(b,c) must be distinct source-node neighbors")
    rows=[]
    for r,omitted in enumerate(neighbors):
        base=tuple(x for x in neighbors if x!=omitted)
        face_sign=-1 if r%2 else 1
        for perm in itertools.permutations(base):
            a,pb,pc=perm
            if (pb,pc)!=(b,c):
                continue
            sign=face_sign*FULL.parity(base,perm)
            rows.append({
                "a":int(a),"b":int(b),"c":int(c),"ordered_edges":[int(a),int(b),int(c)],
                "omitted_neighbor":int(omitted),"face_index":int(r),"base_face":[int(x) for x in base],
                "sign":int(sign),
            })
    if len(rows)!=2:
        raise RuntimeError(f"expected exactly two outer triples for prefix {(b,c)}, found {len(rows)}")
    if len({tuple(x["ordered_edges"]) for x in rows})!=2:
        raise RuntimeError("duplicate outer triples")
    return rows


def run(source_v=0,input_index=0,b=1,c=2,allow_any_prefix=False):
    b=int(b);c=int(c)
    if not allow_any_prefix and (b,c) not in CANDIDATE_PREFIXES:
        raise ValueError(f"prefix {(b,c)} not in frozen six-prefix candidate set {CANDIDATE_PREFIXES}")

    ZVM.patch_and_clear()
    basis=FULL.RAW.PW.basis_full_jhalf()
    if len(basis)!=32:
        raise RuntimeError(f"expected 32 logical inputs, found {len(basis)}")
    if not (0<=input_index<len(basis)):
        raise ValueError("input_index outside logical basis")
    initial=basis[input_index]
    psi=FULL.RAW.CV.gauss_to_covariant({initial:1+0j},source_v)
    plan=epsilon_outer_terms(source_v,b,c)
    term_states=[{} for _ in plan]
    diagmax={
        "CV_complete_basis_leakage":0.0,
        "CK_outer_complete_basis_leakage":0.0,
        "CK_internal_volume_sector_leakage":0.0,
        "CK_complete_charge_basis_leakage":0.0,
    }
    paths=[];max_spin=0.0;middle_nonzero_paths=0

    old,caches=FULL.install_sine_ordering()
    try:
        for i,j,k in itertools.product(range(2),repeat=3):
            s1,leakV=FULL.RAW.COMP.C_volume_component(psi,source_v,c,k,i,FULL.JMAX2)
            diagmax["CV_complete_basis_leakage"]=max(diagmax["CV_complete_basis_leakage"],float(leakV))
            if s1:
                s2,d2=FULL.RAW.KCOMP.C_K_component(s1,source_v,b,j,k,FULL.JMAX2)
                FULL.update_diag(diagmax,d2)
            else:
                s2={}
            if s2:
                middle_nonzero_paths+=1
            # Frozen exact scalar-channel selection from the preregistered full gate.
            s2_scalar={key:amp for key,amp in s2.items() if key[2] in (0,2)}
            outer_support=[]
            for ti,row in enumerate(plan):
                if s2_scalar:
                    s3,d3=FULL.RAW.KCOMP.C_K_component(s2_scalar,source_v,row["a"],i,j,FULL.JMAX2)
                    FULL.update_diag(diagmax,d3)
                else:
                    s3={}
                FULL.add(term_states[ti],s3)
                outer_support.append({"a":row["a"],"support":len(s3),"norm":FULL.norm(s3)})
                max_spin=max(max_spin,FULL.max_spin(s3))
            max_spin=max(max_spin,FULL.max_spin(s1),FULL.max_spin(s2),FULL.max_spin(s2_scalar))
            paths.append({
                "indices":[int(i),int(j),int(k)],
                "after_CV_support":len(s1),"after_CV_norm":FULL.norm(s1),
                "after_middle_CK_support":len(s2),"after_middle_CK_norm":FULL.norm(s2),
                "after_middle_CK_scalar_relevant_support":len(s2_scalar),
                "outer_actions":outer_support,
            })
        cache_info={name:{"hits":fun.cache_info().hits,"misses":fun.cache_info().misses,"currsize":fun.cache_info().currsize} for name,fun in caches.items()}
    finally:
        FULL.restore_ordering(old)

    partial={};term_rows=[]
    max_term_nonscalar_fraction=0.0;max_term_nonscalar_norm_if_near_zero=0.0
    for row,state in zip(plan,term_states):
        sd=FULL.scalar_diagnostics(state)
        if sd["norm"]>FULL.NONZERO_TOL:
            max_term_nonscalar_fraction=max(max_term_nonscalar_fraction,sd["nonscalar_weight_fraction"])
        else:
            max_term_nonscalar_norm_if_near_zero=max(max_term_nonscalar_norm_if_near_zero,sd["nonscalar_norm"])
        FULL.add(partial,state,scale=row["sign"])
        term_rows.append({**row,"support":len(state),"norm":FULL.norm(state),"scalar_diagnostics":sd,"max_spin_reached":FULL.max_spin(state)})

    partial_scalar=FULL.scalar_diagnostics(partial)
    hard={
        "exactly_eight_aux_paths":len(paths)==8,
        "exactly_two_outer_triples":len(term_rows)==2,
        "outer_triples_share_requested_prefix":all(tuple(x["ordered_edges"][1:])==(b,c) for x in term_rows),
        "all_outer_terms_scalar_within_frozen_threshold":max_term_nonscalar_fraction<1e-8 and max_term_nonscalar_norm_if_near_zero<FULL.NONZERO_TOL,
        "signed_prefix_sum_scalar_within_frozen_threshold":FULL.scalar_ok(partial_scalar),
        "CV_complete_basis_leakage_below_1e-9":diagmax["CV_complete_basis_leakage"]<1e-9,
        "CK_outer_complete_basis_leakage_below_1e-9":diagmax["CK_outer_complete_basis_leakage"]<1e-9,
        "CK_internal_volume_sector_leakage_below_1e-9":diagmax["CK_internal_volume_sector_leakage"]<1e-9,
        "spin_cutoff_respected":max_spin<=FULL.JMAX2/2+1e-12,
    }
    conv=convention_descriptor(source_v);hab=habitat_descriptor(source_v)
    full_domain_hash=boundary_domain_hash(basis)
    status=("PREFIX_ZERO_AFTER_FULL_OUTER_SUM" if not partial else "PREFIX_NONZERO_OUTGOING_PARTIAL")
    return {
        "schema":"BQG_LORENTZIAN_PRUNED_PREFIX_WORKER_V1",
        "status":"exact pruned Lorentzian epsilon prefix worker",
        "passed":bool(all(hard.values())),"science_status":status,
        "source_node":int(source_v),"input_logical_basis_index":int(input_index),"input_K_labels":list(initial[1]),
        "prefix":{"b":b,"c":c},"candidate_prefix_set":[list(x) for x in CANDIDATE_PREFIXES],
        "Jmax":FULL.JMAX2/2,"middle_nonzero_aux_path_count":middle_nonzero_paths,
        "outer_term_count":len(term_rows),"outer_terms":term_rows,"paths":paths,
        "signed_partial_support":len(partial),"signed_partial_norm":FULL.norm(partial),
        "signed_partial_scalar_diagnostics":partial_scalar,"max_spin_reached":max_spin,
        "max_diagnostics":diagmax,"hard_integrity_checks":hard,"runtime_exact_cache":cache_info,
        "habitat_descriptor":hab,"habitat_hash":canonical_hash(hab),
        "boundary_domain_hash":full_domain_hash,"convention_descriptor":conv,"convention_hash":canonical_hash(conv),
        "state":encode_state(partial),
        "historical_primitive_charge_basis_diagnostic":{"value":diagmax["CK_complete_charge_basis_leakage"],"hard_acceptance":False},
        "claim_boundary":"One exact signed two-triple partial of raw H_L on one boundary column. It is not yet the full H_L column, not an HH-safe habitat certificate, and implies no physical projector or cosmology.",
    }


def plan_only(source_v=0):
    all_pairs=[]
    neighbors=tuple(FULL.RAW.PW.NEIG[source_v])
    for b in neighbors:
        for c in neighbors:
            if b==c: continue
            rows=epsilon_outer_terms(source_v,b,c)
            all_pairs.append({"prefix":[b,c],"triples":rows})
    triples=[tuple(t["ordered_edges"]) for row in all_pairs for t in row["triples"]]
    return {
        "passed":len(all_pairs)==12 and len(triples)==24 and len(set(triples))==24,
        "source_node":source_v,"neighbors":list(neighbors),"ordered_prefix_count":len(all_pairs),
        "ordered_triple_count":len(triples),"unique_ordered_triple_count":len(set(triples)),"pairs":all_pairs,
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-node",type=int,default=0);ap.add_argument("--input-index",type=int,default=0)
    ap.add_argument("--b",type=int);ap.add_argument("--c",type=int);ap.add_argument("--allow-any-prefix",action="store_true")
    ap.add_argument("--plan-only",action="store_true");ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    if a.plan_only:
        out=plan_only(a.source_node)
    else:
        if a.b is None or a.c is None: ap.error("--b and --c required unless --plan-only")
        out=run(a.source_node,a.input_index,a.b,a.c,a.allow_any_prefix)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    summary={k:v for k,v in out.items() if k not in ("state","paths","pairs")};print(json.dumps(summary,indent=2))
    return 0 if out.get("passed",False) else 1

if __name__=="__main__": raise SystemExit(main())
