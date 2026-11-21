#!/usr/bin/env python3
"""Aggregate six surviving Lorentzian outer-continuation packets.

Together with the six previously certified pathwise-zero middle-prefix classes,
the six surviving prefix packets reconstruct all 24 signed epsilon terms for
the preregistered first logical input.  The output serializes the resulting
Gauss state so it can be reused later in E/L master-pencil overlaps.

A zero/nonzero result is strictly a FIRST LOGICAL COLUMN result.  It is not a
global H_L statement and it is not a physical-projector certificate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL

SURVIVING_PAIR_INDICES = (0, 1, 3, 4, 6, 7)
ZERO_PREFIXES = ((3,4),(2,4),(1,4),(4,3),(4,2),(4,1))
ZERO_PREFIX_EVIDENCE_RUN = 33962645194


def encode_state(state):
    rows=[]
    for key,amp in sorted(state.items(),key=lambda kv:repr(kv[0])):
        spins,Kother,J2,M2,K12,K34=key
        rows.append({
            "spins":list(spins),"Kother":list(Kother),"J2":int(J2),"M2":int(M2),
            "K12":int(K12),"K34":int(K34),"amp":[float(amp.real),float(amp.imag)]
        })
    return rows


def decode_state(rows):
    out={}
    for row in rows:
        key=(tuple(int(x) for x in row["spins"]),tuple(int(x) for x in row["Kother"]),
             int(row["J2"]),int(row["M2"]),int(row["K12"]),int(row["K34"]))
        amp=complex(float(row["amp"][0]),float(row["amp"][1]))
        if abs(amp)>FULL.TOL:
            out[key]=amp
    return out


def state_hash(state):
    raw=json.dumps(encode_state(state),separators=(",",":"),sort_keys=True).encode()
    return hashlib.sha256(raw).hexdigest()


def load_packets(directory:Path):
    files=sorted(directory.glob("continuation_*.json"))
    return files,[json.loads(p.read_text(encoding="utf-8")) for p in files]


def aggregate(packets,files):
    errors=[]
    by_index={}
    for path,p in zip(files,packets):
        if p.get("passed") is not True:
            errors.append(f"{path.name}: continuation packet did not pass")
            continue
        pair=p.get("ordered_pair",{})
        idx=int(pair.get("pair_index",-1))
        if idx in by_index:
            errors.append(f"duplicate pair_index {idx}")
        by_index[idx]=p
    if set(by_index)!=set(SURVIVING_PAIR_INDICES):
        errors.append(f"expected pair indices {SURVIVING_PAIR_INDICES}, got {tuple(sorted(by_index))}")
    if errors:
        return {"passed":False,"science_status":"INCOMPLETE_FIRST_COLUMN_PACKET","errors":errors}

    source_nodes={int(p["source_node"]) for p in by_index.values()}
    inputs={int(p["input_logical_basis_index"]) for p in by_index.values()}
    jmax={float(p["Jmax"]) for p in by_index.values()}
    prefix_runs={p.get("upstream_prefix_provenance",{}).get("run_id") for p in by_index.values()}
    prefix_heads={p.get("upstream_prefix_provenance",{}).get("head_sha") for p in by_index.values()}
    if len(source_nodes)!=1: errors.append("continuation packets disagree on source node")
    if len(inputs)!=1: errors.append("continuation packets disagree on logical input")
    if jmax!={FULL.JMAX2/2}: errors.append("continuation packets disagree on frozen Jmax")
    if len(prefix_runs)!=1: errors.append("continuation packets disagree on upstream prefix run")
    if len(prefix_heads)!=1: errors.append("continuation packets disagree on upstream prefix head")

    continued_terms=sum(len(p.get("continued_full_terms",[])) for p in by_index.values())
    if continued_terms!=12: errors.append(f"expected 12 continued triples, found {continued_terms}")
    if len(ZERO_PREFIXES)*2!=12: errors.append("frozen zero-prefix ledger must cover 12 triples")
    if errors:
        return {"passed":False,"science_status":"INCONSISTENT_FIRST_COLUMN_PACKET","errors":errors}

    total={}
    pair_rows=[]
    for idx in sorted(by_index):
        p=by_index[idx]
        s=decode_state(p.get("signed_pair_state",[]))
        FULL.add(total,s)
        pair_rows.append({
            "pair_index":idx,
            "b":int(p["ordered_pair"]["b"]),"c":int(p["ordered_pair"]["c"]),
            "science_status":p["science_status"],
            "signed_pair_support":len(s),"signed_pair_norm":FULL.norm(s),
            "signed_pair_state_sha256":state_hash(s),
        })

    scalar=FULL.scalar_diagnostics(total)
    gauss,map_diag=FULL.project_covariant_J0_to_gauss(total,next(iter(source_nodes)))
    logical=FULL.logical_projection(gauss)
    total_norm=FULL.norm(total)
    first_zero=total_norm<=FULL.NONZERO_TOL
    hard={
        "six_surviving_prefix_packets_complete":len(by_index)==6,
        "twelve_outer_terms_continued":continued_terms==12,
        "twelve_other_terms_have_prior_pathwise_zero_certificate":len(ZERO_PREFIXES)*2==12,
        "final_trace_state_scalar":FULL.scalar_ok(scalar),
        "covariant_J0_projection_has_no_invalid_keys":len(map_diag.get("invalid_J0_covariant_keys",[]))==0,
        "single_source_node":len(source_nodes)==1,
        "single_logical_input":len(inputs)==1,
        "single_prefix_run":len(prefix_runs)==1,
        "single_prefix_head":len(prefix_heads)==1,
    }
    status=("FIRST_LOGICAL_COLUMN_ZERO_NO_GLOBAL_CONCLUSION" if first_zero
            else "FIRST_LOGICAL_COLUMN_NONZERO_NO_GLOBAL_CONCLUSION")
    return {
        "status":"reconstructed full signed raw Lorentzian first logical column",
        "passed":bool(all(hard.values())),
        "science_status":status,
        "source_node":next(iter(source_nodes)),
        "input_logical_basis_index":next(iter(inputs)),
        "Jmax":FULL.JMAX2/2,
        "upstream_nonzero_prefix_run_id":next(iter(prefix_runs)),
        "upstream_nonzero_prefix_head_sha":next(iter(prefix_heads)),
        "pathwise_zero_prefix_certificate":{
            "evidence_run_id":ZERO_PREFIX_EVIDENCE_RUN,
            "evidence_document":"PETER_WEYL_LORENTZIAN_MIDDLE_PREFIX_RESULT.md",
            "ordered_prefixes":[list(x) for x in ZERO_PREFIXES],
            "implied_zero_ordered_triples":12,
        },
        "continued_surviving_prefixes":pair_rows,
        "continued_ordered_triples":continued_terms,
        "total_ordered_triples_accounted_for":continued_terms+12,
        "raw_covariant_support":len(total),
        "raw_covariant_norm":total_norm,
        "raw_covariant_state_sha256":state_hash(total),
        "raw_scalar_diagnostics":scalar,
        "gauss_projection_support":len(gauss),
        "gauss_projection_norm":FULL.norm(gauss),
        "gauss_projection_state_sha256":state_hash(gauss),
        "gauss_projection_map_diagnostics":map_diag,
        "gauss_state":encode_state(gauss),
        "logical_return_support":len(logical),
        "logical_return_norm":FULL.norm(logical),
        "logical_return_state_sha256":state_hash(logical),
        "logical_return_state":encode_state(logical),
        "hard_integrity_checks":hard,
        "claim_boundary":"Exactly one preregistered logical input column at the frozen cutoff/order. Even a nonzero column is not a global H_L result; even a zero column does not imply the remaining 31 columns vanish. No physical projector or cosmological claim follows.",
        "next_use":"The serialized Gauss output is eligible as microscopic H_L column data for future E/L master-pencil overlap construction once the required column family and HDA/Dtarget production certificates are complete.",
        "errors":[],
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--directory",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    files,packets=load_packets(args.directory)
    out=aggregate(packets,files)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    slim={k:v for k,v in out.items() if k not in ("gauss_state","logical_return_state")}
    print(json.dumps(slim,indent=2))
    return 0 if out.get("passed") else 1


if __name__=="__main__":
    raise SystemExit(main())
