#!/usr/bin/env python3
"""Compare two Lorentzian sparse outgoing artifacts without physics inference.

Intended use: independently computed ungrouped and grouped-k continuations of
the same frozen Lorentzian prefix/column.  The gate compares the decoded sparse
state on the union support and reports max-absolute and l2-relative differences,
plus invariant metadata.  It never promotes an artifact to P_phys.
"""
from __future__ import annotations

import argparse, json, math
from pathlib import Path

ABS_TOL=5.0e-9
REL_TOL=5.0e-8


def decode(rows):
    out={}
    for r in rows:
        if "Kother" in r:
            key=(tuple(r["spins"]),tuple(r["Kother"]),int(r["J2"]),int(r["M2"]),int(r["K12"]),int(r["K34"]))
        elif "K_labels" in r:
            key=(tuple(r["spins"]),tuple(r["K_labels"]))
        else:
            raise ValueError("unsupported sparse row schema")
        out[key]=out.get(key,0j)+complex(float(r["amp"][0]),float(r["amp"][1]))
    return out


def norm(state):
    return math.sqrt(sum(abs(z)**2 for z in state.values()))


def load(path):
    d=json.loads(Path(path).read_text(encoding="utf-8"))
    rows=d.get("state")
    if rows is None and "complete_gauss_outgoing_column" in d:
        rows=d["complete_gauss_outgoing_column"]["state"]
    if rows is None: raise ValueError(f"{path}: no sparse state")
    return d,decode(rows)


def run(a_path,b_path):
    A,a=load(a_path); B,b=load(b_path)
    keys=set(a)|set(b)
    delta={k:a.get(k,0j)-b.get(k,0j) for k in keys}
    delta={k:z for k,z in delta.items() if z!=0j}
    na=norm(a); nb=norm(b); nd=norm(delta)
    maxabs=max((abs(z) for z in delta.values()),default=0.0)
    rel=nd/max(na,nb,1e-300)
    metadata_checks={
        "source_node_match":A.get("source_node")==B.get("source_node"),
        "input_index_match":A.get("input_logical_basis_index")==B.get("input_logical_basis_index"),
        "prefix_match":A.get("prefix")==B.get("prefix"),
        "Jmax_match":A.get("Jmax")==B.get("Jmax"),
        "habitat_hash_match":A.get("habitat_hash")==B.get("habitat_hash"),
        "boundary_domain_hash_match":A.get("boundary_domain_hash")==B.get("boundary_domain_hash"),
        "convention_hash_match":A.get("convention_hash")==B.get("convention_hash"),
    }
    checks={
        **metadata_checks,
        "max_abs_difference_below_tolerance":maxabs<ABS_TOL,
        "relative_l2_difference_below_tolerance":rel<REL_TOL,
    }
    return {
        "schema":"BQG_LORENTZIAN_SPARSE_EQUIVALENCE_V1",
        "passed":bool(all(checks.values())),
        "absolute_tolerance":ABS_TOL,"relative_tolerance":REL_TOL,
        "support_a":len(a),"support_b":len(b),"union_support":len(keys),"difference_support":len(delta),
        "norm_a":na,"norm_b":nb,"difference_norm":nd,"relative_l2_difference":rel,"max_abs_difference":maxabs,
        "checks":checks,
        "claim_boundary":"Numerical equivalence check between two implementations of the same frozen sparse microscopic action only; no projector or cosmological claim follows."
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--a",type=Path,required=True); ap.add_argument("--b",type=Path,required=True); ap.add_argument("--output",type=Path,required=True)
    x=ap.parse_args(); out=run(x.a,x.b); x.output.parent.mkdir(parents=True,exist_ok=True); x.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8"); print(json.dumps(out,indent=2)); return 0 if out["passed"] else 1

if __name__=="__main__": raise SystemExit(main())
