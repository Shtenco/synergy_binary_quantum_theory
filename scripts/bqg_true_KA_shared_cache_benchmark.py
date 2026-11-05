#!/usr/bin/env python3
"""True K5 KA CK: benchmark shared in-process SU2 recoupling cache.

Uses ORIGINAL frozen KB state and only two already SHA-verified small columns
(0,12). No operator shortcuts. Shared-vs-isolated executions both must match
the historical full outgoing complex vectors.
"""
import argparse, json, math, sys, time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import bqg_lorentzian_checkpointed_triple as Base
import bqg_node1_ka_partial_gram as Gram
import peter_weyl_lorentzian_logical_projection_gate as LP
import peter_weyl_lorentzian_ordered_triple_gate as RAW

def run(mode):
    root=Path("checkpoint")
    m=json.loads((root/"kb_000.json").read_text(encoding="utf-8"))
    kb,_=Base.load_result(root,"kb",(0,0,0),m["freeze"])
    items=sorted(kb.items(),key=lambda kv:repr(kv[0]))
    assert len(items)==120 and m["freeze"]["source"]==1
    assert m["freeze"]["triple"]==[0,2,3] and m["freeze"]["JMAX2"]==7
    LP.JMAX2=Base.JMAX2
    selected=(0,12)
    elapsed=[];errors=[];supports=[]
    t0=time.perf_counter()
    restore=None
    if mode=="shared":restore,_=LP.install_sine_cached_stack()
    try:
        for idx in selected:
            restore_one=None
            if mode=="isolated":restore_one,_=LP.install_sine_cached_stack()
            try:
                t1=time.perf_counter()
                out,diag=RAW.KCOMP.C_K_component({items[idx][0]:items[idx][1]},1,0,0,0,Base.JMAX2)
                seconds=time.perf_counter()-t1
            finally:
                if restore_one:restore_one()
            src,old=Gram.vector(root,idx)
            archive={(k[0],k[1],*k[2]):v for k,v in src.items()}
            keys=set(archive)|set(out)
            err=math.sqrt(sum(abs(out.get(k,0)-archive.get(k,0))**2 for k in keys))
            ref=math.sqrt(sum(abs(v)**2 for v in archive.values()))
            assert err/max(ref,1e-30)<1e-8,("mismatch",mode,idx,err)
            elapsed.append(seconds);errors.append(err/max(ref,1e-30))
            supports.append([len(out),len(archive)])
            print("TRUE_CK_CACHE_BENCH_KEY",idx,mode,seconds,len(out),errors[-1],flush=True)
    finally:
        if restore:restore()
    return dict(mode=mode,status="PASS",indices=list(selected),individual_seconds=elapsed,
      total_time_seconds=time.perf_counter()-t0,
      direct_vs_SHA_archived_relative_errors=errors,
      supports=supports,
      note="Measured only two real small sources. Cache performance for high-spin, different sources is not established.")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=("shared","isolated"),required=True)
    a=ap.parse_args()
    print("BQG_REAL_CK_CACHE_BENCH_BEGIN")
    print(json.dumps(run(a.mode),indent=2))
    print("BQG_REAL_CK_CACHE_BENCH_END")
