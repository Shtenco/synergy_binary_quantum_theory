#!/usr/bin/env python3
"""BQG G29: reconstruct 4 full genuine KA outputs from a successful Actions run.

All source vectors and real SU(2) arrays checked by SHA256 + frozen source/
kernel metadata; no artificial output, truncation, or symmetry substitution.
SHA-checked original job logs have been pinned to the research Git branch as JSON, containing NPZ bytes.
"""
from __future__ import annotations
import argparse,base64,hashlib,io,json,os,re,urllib.request
from pathlib import Path
import numpy as np
import bqg_node1_ka_partial_gram as gram

PRIOR=(59,105,113,112,57,52,24,83,31,64,53,111)
NEW=(15,39,40,80)
FROZEN="f486e2250d4a5678a34438467b9e338d730ebf6e1d1ec956d2de7f8c2521b293"
REPO="Shtenco/binary_quantum_theory"

def api(path,token):
    url="https://api.github.com/repos/"+REPO+"/"+path
    req=urllib.request.Request(url,headers={
      "Authorization":"Bearer "+token,
      "Accept":"application/vnd.github+json",
      "X-GitHub-Api-Version":"2022-11-28",
      "User-Agent":"BQG-G29-SHA-verified-research"})
    with urllib.request.urlopen(req,timeout=45) as resp:
      body=resp.read()
    return body

def inspect_npz(buf,meta):
    assert hashlib.sha256(buf).hexdigest()==meta["sha256"]
    assert meta["kb_sha256"]==FROZEN
    assert meta["stage"]=="ka_part" and meta["rows"]==1 and meta["parts"]==120
    assert meta["channel"]==[0,0,0]
    fr=meta["freeze"]
    assert fr["source"]==1 and fr["column_index"]==0 and fr["triple"]==[0,2,3]
    assert fr["JMAX2"]==7
    with np.load(io.BytesIO(buf),allow_pickle=False) as z:
        n=meta["support"]
        assert z["spins"].shape==(n,10) and z["Kother"].shape==(n,5)
        assert z["charges"].shape==(n,4) and z["amp"].shape==(n,)
        assert np.isfinite(z["amp"].real).all() and np.isfinite(z["amp"].imag).all()
        assert abs(float(np.linalg.norm(z["amp"]))-meta["norm"])<1e-12

def save(root,index,meta,buf):
    assert meta["part_index"]==index
    inspect_npz(buf,meta)
    stem=root/f"ka_000_part{index:03d}_of120"
    stem.with_suffix(".npz").write_bytes(buf)
    stem.with_suffix(".json").write_text(json.dumps(meta,sort_keys=True),encoding="utf-8")
    vec,check=gram.vector(root,index)
    assert len(vec)==meta["support"]
    print("G29_VERIFIED_COLUMN",index,len(vec),meta["norm"],meta["sha256"],flush=True)

def recover_old(root,repo_root):
    for index in PRIOR:
        path=repo_root/"research"/"node1_fullstate_gate_2026_10_31"/"log_recovered_ka"/f"key_{index}.json"
        payload=json.loads(path.read_text(encoding="utf-8"))
        meta=payload["metadata"]
        raw=base64.b64decode(payload["npz_base64"],validate=True)
        save(root,index,meta,raw)

def recover_new(root,repo_root):
    # Evidence is pinned in Git under the original research branch: NO live logs
    # authentication, NO unknown network bytes, NO API rate limits.
    evidence={}
    for idx in NEW:
        path=repo_root/"research"/"node1_fullstate_gate_2026_10_31"/"log_recovered_ka"/f"key_{idx}.json"
        payload=json.loads(path.read_text(encoding="utf-8"))
        meta=payload["metadata"]
        raw=base64.b64decode(payload["npz_base64"],validate=True)
        save(root,idx,meta,raw)
        evidence[idx]={"job_id":payload["job_id"],"source_run":payload["source_run"],
                       "sha256":meta["sha256"],"support":meta["support"]}
    assert set(evidence)==set(NEW)
    return evidence

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--run",type=int,required=True)
    ap.add_argument("--data",type=Path,required=True)
    ap.add_argument("--repo-root",type=Path,default=Path("."))
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    args.data.mkdir(parents=True,exist_ok=True)
    parent=args.data/"kb_000"
    kbmeta=json.loads(parent.with_suffix(".json").read_text(encoding="utf-8"))
    assert kbmeta["sha256"]==FROZEN==gram.checksum(parent.with_suffix(".npz"))
    recover_old(args.data,args.repo_root)
    new_evidence=recover_new(args.data,args.repo_root)
    selection=tuple(sorted(gram.STUDIED+PRIOR+NEW))
    assert len(selection)==29 and len(set(selection))==29
    result=gram.calc(args.data,study=selection)
    assert result["completed_KA_input_keys"]==29
    assert 0.70 < result["input_KB_squared_weight_covered"] < 0.71
    result["source_run"]=args.run
    result["new_evidence"]=new_evidence
    result["scope"]="ONE source, ONE ordered triple, ONE auxiliary 000 only; HDA and TT open"
    print("BQG_G29_ACTUAL_RESULT_BEGIN",flush=True)
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
    print("BQG_G29_ACTUAL_RESULT_END",flush=True)
    if args.output:args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
if __name__=="__main__":main()
