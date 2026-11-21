#!/usr/bin/env python3
"""Deterministically assemble the first raw Lorentzian outgoing column.

Twelve ordered triples are evaluated through six exact outer-from-middle worker
artifacts.  The complementary twelve triples are admitted only because the
frozen middle-prefix result records pathwise-zero C_b(K)C_c(V) states for their
six prefixes.  The union must reproduce all 24 epsilon-oriented ordered terms.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_lorentzian_pruned_prefix_worker as PLAN

ZERO_PREFIXES=((1,4),(2,4),(3,4),(4,1),(4,2),(4,3))
SURVIVING_PREFIXES=PLAN.CANDIDATE_PREFIXES


def file_sha256(path: Path):
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


def verify_zero_prefix_record(path: Path):
    text=path.read_text(encoding="utf-8")
    checks={
        "records_source_node_v0":"source node `v=0`" in text,
        "records_input_index_0":"logical input index `0`" in text,
        "records_Jmax_7_over_2":"`Jmax=7/2`" in text,
        "records_pathwise_zero_equation":"C_b(K)C_c(V)|\\psi_0\\rangle=0" in text,
        "records_twelve_of_twentyfour_zero":"twelve of the twenty-four" in text,
    }
    for b,c in ZERO_PREFIXES:
        checks[f"zero_prefix_{b}_{c}_present"]=(f"`({b},{c})`" in text or f"`({b}, {c})`" in text)
    if not all(checks.values()):
        raise RuntimeError(f"zero-prefix provenance document failed frozen checks: {checks}")
    return {"path":str(path),"sha256":file_sha256(path),"checks":checks}


def run(worker_paths, zero_doc: Path):
    zero_provenance=verify_zero_prefix_record(zero_doc)
    workers=[]; seen=set(); common={}; total={}; worker_hashes=[]; explicit_triples=[]
    for path in sorted(worker_paths,key=lambda p:str(p)):
        data=json.loads(path.read_text(encoding="utf-8"))
        if data.get("schema")!="BQG_LORENTZIAN_OUTER_FROM_MIDDLE_V1" or not data.get("passed",False):
            raise RuntimeError(f"invalid outer worker artifact {path}")
        p=data.get("prefix",{}); prefix=(int(p["b"]),int(p["c"]))
        if prefix in seen: raise RuntimeError(f"duplicate surviving prefix {prefix}")
        seen.add(prefix)
        if prefix not in SURVIVING_PREFIXES: raise RuntimeError(f"unexpected surviving prefix {prefix}")
        for k in ("source_node","input_logical_basis_index","Jmax","habitat_hash","boundary_domain_hash","convention_hash"):
            if k not in common: common[k]=data[k]
            elif data[k]!=common[k]: raise RuntimeError(f"worker mismatch for {k}: {path}")
        state=decode_state(data.get("state",[])); FULL.add(total,state)
        terms=[tuple(int(x) for x in r["ordered_edges"]) for r in data.get("outer_terms",[])]
        if len(terms)!=2 or len(set(terms))!=2: raise RuntimeError(f"bad outer term coverage {path}")
        explicit_triples.extend(terms)
        worker_hashes.append({"prefix":list(prefix),"path":str(path),"sha256":file_sha256(path),"middle_artifact":data.get("middle_artifact")})
        workers.append({
            "prefix":list(prefix),"science_status":data.get("science_status"),
            "signed_partial_support":int(data.get("signed_partial_support",0)),
            "signed_partial_norm":float(data.get("signed_partial_norm",0.0)),
            "max_spin_reached":float(data.get("max_spin_reached",0.0)),
            "outer_terms":[list(x) for x in terms],
        })
    if seen!=set(SURVIVING_PREFIXES):
        raise RuntimeError(f"surviving prefix coverage incomplete: have {sorted(seen)}, need {sorted(SURVIVING_PREFIXES)}")

    zero_triples=[]
    for b,c in ZERO_PREFIXES:
        zero_triples.extend(tuple(r["ordered_edges"]) for r in PLAN.epsilon_outer_terms(int(common["source_node"]),b,c))
    all_triples=explicit_triples+zero_triples
    coverage_checks={
        "six_surviving_prefixes_once":len(seen)==6,
        "twelve_explicit_outer_triples":len(explicit_triples)==12 and len(set(explicit_triples))==12,
        "twelve_certified_zero_outer_triples":len(zero_triples)==12 and len(set(zero_triples))==12,
        "all_24_ordered_triples_unique":len(all_triples)==24 and len(set(all_triples))==24,
    }

    basis=FULL.RAW.PW.basis_full_jhalf()
    initial=basis[int(common["input_logical_basis_index"])]
    scalar=FULL.scalar_diagnostics(total)
    gauss,mapdiag=FULL.project_covariant_J0_to_gauss(total,int(common["source_node"]))
    logical=FULL.logical_projection(gauss)
    logical_norm=FULL.norm(logical)
    initial_amp=logical.get(initial,0j)
    logical_rows=[]
    for idx,key in enumerate(basis):
        amp=logical.get(key,0j)
        if abs(amp)>FULL.TOL:
            logical_rows.append({"logical_basis_index":idx,"K_labels":list(key[1]),"amp":[float(amp.real),float(amp.imag)],"abs":abs(amp)})

    max_spin=max([float(x["max_spin_reached"]) for x in workers]+[FULL.max_spin(total)])
    hard={
        **coverage_checks,
        "common_source_node_is_0":int(common["source_node"])==0,
        "common_input_index_is_0":int(common["input_logical_basis_index"])==0,
        "common_Jmax_is_7_over_2":abs(float(common["Jmax"])-3.5)<1e-12,
        "full_signed_output_scalar_within_frozen_threshold":FULL.scalar_ok(scalar),
        "spin_cutoff_respected":max_spin<=3.5+1e-12,
        "J0_reverse_projection_has_no_invalid_keys":not mapdiag["invalid_J0_covariant_keys"],
        "J0_reverse_projection_has_no_collisions":int(mapdiag["mapping_collisions"])==0,
    }
    if not total:
        science="FULL_FIRST_RAW_HL_COLUMN_ZERO"
    elif logical_norm>FULL.NONZERO_TOL:
        science="FULL_FIRST_RAW_HL_COLUMN_NONZERO_WITH_LOGICAL_RETURN"
    else:
        science="FULL_FIRST_RAW_HL_COLUMN_NONZERO_LOGICAL_RETURN_ZERO"

    return {
        "schema":"BQG_LORENTZIAN_PRUNED_FIRST_COLUMN_V1",
        "passed":bool(all(hard.values())),"science_status":science,
        **common,
        "operator":"H_L^raw ~ epsilon Tr_aux[C(K)C(K)C(V)] under frozen sine ordering",
        "coverage":{"surviving_prefixes":[list(x) for x in SURVIVING_PREFIXES],"pathwise_zero_prefixes":[list(x) for x in ZERO_PREFIXES],
                    "explicit_ordered_triples":[list(x) for x in sorted(explicit_triples)],"certified_zero_ordered_triples":[list(x) for x in sorted(zero_triples)]},
        "coverage_checks":coverage_checks,"zero_prefix_provenance":zero_provenance,"worker_artifacts":worker_hashes,"workers":workers,
        "full_outgoing_support":len(total),"full_outgoing_norm":FULL.norm(total),"full_scalar_diagnostics":scalar,"max_spin_reached":max_spin,
        "gauss_reverse_projection":{"support":len(gauss),"norm":FULL.norm(gauss),"diagnostics":mapdiag},
        "logical_return":{"support":len(logical),"norm":logical_norm,"fraction_of_full_norm":logical_norm/max(FULL.norm(total),1e-300),
                          "initial_return_amplitude":[float(initial_amp.real),float(initial_amp.imag)],"nonzero_amplitudes":logical_rows},
        "hard_integrity_checks":hard,"state":PLAN.encode_state(total),
        "claim_boundary":"Complete first raw Lorentzian outgoing column at the declared finite cutoff/order, assembled from 12 explicitly continued triples plus 12 independently certified pathwise-zero triples. This is not yet a Hermitian physical H_L convention, not the remaining 31 columns, not an HH-safe habitat certificate, and not P_phys or cosmology.",
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--worker",type=Path,action="append",required=True)
    ap.add_argument("--zero-prefix-result",type=Path,default=Path("PETER_WEYL_LORENTZIAN_MIDDLE_PREFIX_RESULT.md"))
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args(); out=run(a.worker,a.zero_prefix_result)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in out.items() if k not in ("state","worker_artifacts","workers","coverage")},indent=2))
    return 0 if out["passed"] else 1

if __name__=="__main__": raise SystemExit(main())
