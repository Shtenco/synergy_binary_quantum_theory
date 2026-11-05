#!/usr/bin/env python3
"""Exact G17 from the already recorded K5 Lorentzian source snapshots."""
import json,base64,hashlib
from pathlib import Path
import bqg_node1_ka_partial_gram as audit

root=Path("data")
for idx in (59,105,113,24):
    source=Path("research/node1_fullstate_gate_2026_10_31/log_recovered_ka")/f"key_{idx}.json"
    packet=json.loads(source.read_text())
    meta=packet["metadata"]
    buf=base64.b64decode(packet["npz_base64"],validate=True)
    assert hashlib.sha256(buf).hexdigest()==meta["sha256"]
    assert meta["part_index"]==idx and meta["parts"]==120
    assert meta["rows"]==1 and meta["channel"]==[0,0,0]
    stem=root/f"ka_000_part{idx:03d}_of120"
    stem.with_suffix(".npz").write_bytes(buf)
    stem.with_suffix(".json").write_text(json.dumps(meta))
    vec,checked=audit.vector(root,idx)
    print("SHA_VERIFIED",idx,len(vec),checked["norm"],flush=True)
selected=tuple(sorted(audit.STUDIED+(59,105,113,24)))
result=audit.calc(root,study=selected)
assert result["completed_KA_input_keys"]==17
assert result["Gram_rank_numpy"]>=13
assert .39<result["input_KB_squared_weight_covered"]<.40
print("G17_FINAL_START",flush=True)
print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
print("G17_FINAL_END",flush=True)
