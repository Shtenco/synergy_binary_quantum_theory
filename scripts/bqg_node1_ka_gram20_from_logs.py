#!/usr/bin/env python3
"""Exact SHA-checked 20-source K5 Lorentzian Gram; requires archived real outputs."""
import base64,hashlib,json,argparse
from pathlib import Path
import bqg_node1_ka_partial_gram as audit

EXTRA=(59,105,113,112,57,52,24)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--data",type=Path,required=True)
    p.add_argument("--repo-root",type=Path,default=Path("."))
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    assert a.data.is_dir(),"Missing genuine archived baseline outputs"
    for idx in EXTRA:
        source=a.repo_root/"research"/"node1_fullstate_gate_2026_10_31"/"log_recovered_ka"/f"key_{idx}.json"
        if not source.exists():raise RuntimeError("BLOCKED: missing genuine output "+str(source))
        pkt=json.loads(source.read_text(encoding="utf-8"))
        meta=pkt["metadata"]
        buf=base64.b64decode(pkt["npz_base64"],validate=True)
        assert hashlib.sha256(buf).hexdigest()==meta["sha256"],idx
        assert meta["part_index"]==idx and meta["parts"]==120 and meta["rows"]==1
        assert meta["stage"]=="ka_part" and meta["channel"]==[0,0,0]
        stem=a.data/f"ka_000_part{idx:03d}_of120"
        stem.with_suffix(".npz").write_bytes(buf)
        stem.with_suffix(".json").write_text(json.dumps(meta),encoding="utf-8")
        audit.vector(a.data,idx)
    keys=tuple(sorted(audit.STUDIED+EXTRA))
    assert len(keys)==20 and len(set(keys))==20
    result=audit.calc(a.data,study=keys)
    assert result["completed_KA_input_keys"]==20
    assert result["Gram_rank_numpy"]>=16
    print("G20_VERIFIED_START")
    print(json.dumps(result,indent=2,ensure_ascii=False))
    print("G20_VERIFIED_END")
    if a.output:a.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
if __name__=="__main__":main()
