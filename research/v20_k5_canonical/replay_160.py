#!/usr/bin/env python3
"""V20 canonical K5 cold replayer. Recomputes all 160 columns; refuses partial results.
The full 160-column data and standalone executable archive are attached in the chat.
"""
import argparse,hashlib,importlib.util,json
from pathlib import Path
import numpy as np

BLOB="e7e1e0705c4b833dcee167a2036ef5f8929a4388"
SHA256="596c4179cb2c8af20b41c18b3ab87b711c3850c5c5e05b1bc91a1f87c9f11cc4"

def load(path):
    raw=path.read_bytes()
    sha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if sha!=BLOB or hashlib.sha256(raw).hexdigest()!=SHA256:
        raise ValueError("Canonical GitHub K5 backend mismatch")
    spec=importlib.util.spec_from_file_location("k5_v20",path)
    k=importlib.util.module_from_spec(spec);spec.loader.exec_module(k)
    return k

def dump_atomic(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    tmp.replace(path)

def decode(rows):
    out={}
    for d in rows:
        key=(tuple(d["key"][0]),tuple(d["key"][1]))
        z=d["coefficient"]
        if key in out:raise ValueError("Duplicate output key")
        out[key]=complex(float.fromhex(z[0]),float.fromhex(z[1]))
    return out

def compute(k,folder):
    src=k.basis_full_jhalf()
    if len(src)!=32:raise ValueError("Expected 32 K5 sources")
    cols={}
    for i,state in enumerate(src):
        for v in range(5):
            path=folder/f"source_{i:02d}_vertex_{v}_wall5.json"
            if path.exists():
                data=json.loads(path.read_text())
                if data["backend_sha256"]!=SHA256 or data["source_index"]!=i or data["node"]!=v or data["input"]!=[list(state[0]),list(state[1])]:
                    raise ValueError(f"Invalid checkpoint {path}")
                col=decode(data["rows"])
            else:
                col=k.prune_state(k.apply_H_cached_state({state:1.+0j},v,5),1e-10)
                rows=[{"key":[list(q[0]),list(q[1])],
                       "coefficient":[float(a.real).hex(),float(a.imag).hex()]}
                      for q,a in sorted(col.items(),key=lambda it:repr(it[0]))]
                dump_atomic(path,{"backend_sha256":SHA256,
                    "source_index":i,"node":v,"input":[list(state[0]),list(state[1])],
                    "rows":rows})
                print(f"{i}:{v} complete {len(rows)}",flush=True)
            cols[i,v]=col
    if len(cols)!=160:raise RuntimeError("Incomplete 160-column data")
    outputs=set(src)
    for col in cols.values():outputs.update(col)
    order=sorted(outputs,key=repr)
    ix={q:j for j,q in enumerate(order)}
    C=np.zeros((5*len(order),32),complex)
    for i in range(32):
        for v in range(5):
            for q,a in cols[i,v].items():
                C[v*len(order)+ix[q],i]=a
    M=C.conj().T@C
    ev=np.linalg.eigvalsh(M)
    return {"input_count":32,"column_count":160,
            "output_union_with_originals":len(outputs),
            "rank":int(np.sum(ev>1e-10)),
            "min_eigenvalue":float(ev[0]),"max_eigenvalue":float(ev[-1]),
            "trace":float(np.trace(M).real),"physical_projector":False}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--backend",type=Path,default=Path("scripts/k5_peter_weyl_safe_hda_column.py"))
    p.add_argument("--cache",type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(compute(load(a.backend),a.cache),indent=2))
