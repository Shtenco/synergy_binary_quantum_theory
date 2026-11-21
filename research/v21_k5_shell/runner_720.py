#!/usr/bin/env python3
"""V21: exact original K5 first-shell runner, 720 inputs x 5 node operators.

Requires a complete V20 canonical_columns directory. Every completed input
has 5 full canonical H_v images. Never zero-fill an uncomputed column.
"""
from __future__ import annotations
import argparse, collections, concurrent.futures as futures, hashlib
import importlib.util, json, os, time
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
BACKEND=Path(os.environ.get("BQG_K5_BACKEND", str(REPO/"scripts/k5_peter_weyl_safe_hda_column.py")))
INPUT=Path(os.environ.get("BQG_V20_COLUMNS", "/tmp/bqg_v20/canonical_columns"))
OUTPUT=Path(os.environ.get("BQG_V21_OUTPUT", "/tmp/bqg_v21"))
SHA256="596c4179cb2c8af20b41c18b3ab87b711c3850c5c5e05b1bc91a1f87c9f11cc4"
GIT_BLOB="e7e1e0705c4b833dcee167a2036ef5f8929a4388"

def key(x):return tuple(x[0]),tuple(x[1])
def portable(x):return [list(x[0]),list(x[1])]
def parse(rows):
    out={}
    for r in rows:
        k=key(r["key"])
        if k in out:raise ValueError("duplicate K5 output")
        z=r.get("coefficient",r.get("amplitude"))
        out[k]=complex(float.fromhex(z[0]),float.fromhex(z[1]))
    return out
def save(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    t=p.with_suffix(".partial")
    t.write_text(json.dumps(obj,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    t.replace(p)
def backend():
    raw=BACKEND.read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if hashlib.sha256(raw).hexdigest()!=SHA256 or blob!=GIT_BLOB:
        raise RuntimeError("Original GitHub K5 SHA mismatch")
    spec=importlib.util.spec_from_file_location("k5_v21",BACKEND)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

def shell():
    mod=backend();original=mod.basis_full_jhalf()
    if len(original)!=32:raise RuntimeError("Expected 32 original K5 states")
    outputs=set()
    for i,k in enumerate(original):
        for v in range(5):
            p=INPUT/f"source_{i:02d}_vertex_{v}_wall5.json"
            d=json.loads(p.read_text())
            if d["backend_sha256"]!=SHA256 or d["source_index"]!=i or d["node"]!=v or key(d["input"])!=k:
                raise RuntimeError(f"V20 provenance failure: {p}")
            outputs.update(parse(d["rows"]))
    if len(outputs)!=720 or set(original).intersection(outputs):
        raise RuntimeError("Expected 720 distinct, new first-shell states")
    c=collections.Counter(sum(s==0 for s in q[0]) for q in outputs)
    if dict(c)!={0:320,1:240,2:120,3:40}:raise RuntimeError("Wrong first-shell topology")
    ordered=sorted(outputs,key=lambda q:(-sum(s==0 for s in q[0]),sum(q[0]),repr(q)))
    if any(sum(q[0])%2!=1 for q in ordered):raise RuntimeError("Parity violation")
    catalog={"canonical_blob":GIT_BLOB,"sources":[portable(q) for q in ordered]}
    f=OUTPUT/"shell_catalog.json"
    if f.exists() and json.loads(f.read_text())!=catalog:raise RuntimeError("Catalog changed")
    if not f.exists():save(f,catalog)
    return ordered

def target(i,v):return OUTPUT/"columns"/f"shell_{i:03d}_node_{v}_wall5.json"
def check(p,i,v,state):
    d=json.loads(p.read_text())
    if d["canonical_sha256"]!=SHA256 or d["canonical_git_blob"]!=GIT_BLOB or d["index"]!=i or d["node"]!=v or key(d["input"])!=state:
        raise ValueError(f"Corrupt/incompatible checkpoint: {p}")
    c=parse(d["rows"])
    if len(c)!=d["support"] or abs(sum(abs(a)**2 for a in c.values())-d["norm2"])>1e-8:
        raise ValueError("Checkpoint completeness failure")
    return d

def compute(item):
    i,state=item;mod=backend()
    results=[]
    for v in range(5):
        p=target(i,v)
        if p.exists():
            d=check(p,i,v,state)
        else:
            st=time.monotonic()
            r=mod.prune_state(mod.apply_H_cached_state({state:1.+0j},v,5),1e-10)
            if any(sum(k[0])%2!=0 for k in r):raise RuntimeError("H_v selection failed")
            d={"index":i,"node":v,"input":portable(state),"canonical_sha256":SHA256,
               "canonical_git_blob":GIT_BLOB,"support":len(r),
               "norm2":sum(abs(a)**2 for a in r.values()),
               "seconds":time.monotonic()-st,
               "rows":[{"key":portable(k),"coefficient":[float(a.real).hex(),float(a.imag).hex()]}
                       for k,a in sorted(r.items(),key=lambda z:repr(z[0]))]}
            save(p,d)
        results.append(d["support"])
    return (i,results)

def progress(allkeys):
    ready=[];n=0
    for i,state in enumerate(allkeys):
        found=0
        for v in range(5):
            p=target(i,v)
            if p.exists():
                check(p,i,v,state);found+=1;n+=1
        if found==5:ready.append(i)
    stat={"required_inputs":720,"required_columns":3600,
          "complete_inputs":len(ready),"saved_columns":n,
          "complete_input_indices":ready,
          "missing_input_indices":[i for i in range(720) if i not in set(ready)],
          "fully_computed":len(ready)==720 and n==3600,
          "physical_projector_computed":False}
    save(OUTPUT/"progress.json",stat)
    return stat

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--start",type=int,default=0)
    p.add_argument("--stop",type=int,default=720)
    p.add_argument("--max-new",type=int,default=720)
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--report-only",action="store_true")
    a=p.parse_args()
    states=shell()
    if not a.report_only:
        todo=[]
        for i in range(a.start,min(a.stop,720)):
            if all(target(i,v).exists() for v in range(5)):continue
            todo.append((i,states[i]))
            if len(todo)==a.max_new:break
        with futures.ProcessPoolExecutor(max_workers=a.workers) as pool:
            jobs=[pool.submit(compute,x) for x in todo]
            for f in futures.as_completed(jobs):
                i,s=f.result();print("COMPLETE",i,s,flush=True)
    print(json.dumps(progress(states),ensure_ascii=False))
if __name__=="__main__":main()
