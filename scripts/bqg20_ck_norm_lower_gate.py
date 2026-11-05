#!/usr/bin/env python3
"""Exactly measured LOWER bound on complete CK operator norm from 20 true columns.
All 20 weighted source vectors are SHA checked. Real outgoing columns A|e_j> are
formed by dividing by original complex input coefficient a_j; eigenvalues of
their Gram give a strict lower bound on ||A||. It cannot upper-bound remainder.
"""
import argparse,base64,hashlib,json
from pathlib import Path
import numpy as np
import bqg_node1_ka_partial_gram as audit

INDICES=tuple(sorted(audit.STUDIED+(59,105,113,112,57,52,24)))
EXTRA=(59,105,113,112,57,52,24)
def run(root,repo):
    for idx in EXTRA:
        src=repo/"research/node1_fullstate_gate_2026_10_31/log_recovered_ka"/f"key_{idx}.json"
        d=json.loads(src.read_text(encoding="utf-8"));m=d["metadata"]
        raw=base64.b64decode(d["npz_base64"],validate=True)
        assert hashlib.sha256(raw).hexdigest()==m["sha256"]
        stem=root/f"ka_000_part{idx:03d}_of120"
        stem.with_suffix(".npz").write_bytes(raw)
        stem.with_suffix(".json").write_text(json.dumps(m))
    meta=json.loads((root/"kb_000.json").read_text())
    assert meta["sha256"]==audit.checksum(root/"kb_000.npz")
    with np.load(root/"kb_000.npz",allow_pickle=False) as kb:
        coeffs=np.asarray(kb["amp"])
    assert len(coeffs)==120
    v=[]
    for i in INDICES:
        col, m=audit.vector(root,i)
        assert m["kb_sha256"]==meta["sha256"] and m["freeze"]==meta["freeze"]
        assert abs(coeffs[i])>0
        v.append({k:complex(val)/coeffs[i] for k,val in col.items()})
    G=np.zeros((len(INDICES),len(INDICES)),complex)
    for i,a in enumerate(v):
        for j in range(i,len(INDICES)):
            b=v[j]
            G[i,j]=sum(np.conjugate(a[k])*b[k] for k in a.keys()&b.keys())
            G[j,i]=G[i,j].conjugate()
    ev=np.linalg.eigvalsh(G)
    assert ev.min()>0
    lower=float(np.sqrt(ev.max()))
    one=float(max(np.sqrt(G[i,i].real) for i in range(len(INDICES))))
    residual=float(np.linalg.norm(coeffs[[i for i in range(120) if i not in INDICES]]))
    threshold=0.1*0.06524151188803094/residual
    result=dict(test="true 20 measured UNWEIGHTED CK columns",
      measured_CK_norm_lower_bound=lower,
      measured_single_column_norm_lower_bound=one,
      is_global_upper_bound=False,
      lower_bound_exceeds_global_10pct_threshold=bool(lower>threshold),
      global_bound_needed_for_10pct_relative_partial_output=threshold,
      remaining_input_norm=residual,
      physical_HDA="OPEN",physical_TT="OPEN",
      interpretation="If lower_bound > threshold, no sound bound ||CK|| <= threshold can be proven on the full space; must derive a REMAINDER-SUBSPACE bound instead.")
    return result
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True)
    p.add_argument("--repo",type=Path,default=Path("."));p.add_argument("--output",type=Path)
    a=p.parse_args();z=run(a.root,a.repo);s=json.dumps(z,indent=2)
    print("CK_NORM_LOWER_BOUND_START");print(s);print("CK_NORM_LOWER_BOUND_END")
    if a.output:a.output.write_text(s+"\n")
