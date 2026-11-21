#!/usr/bin/env python3
"""Export one complete Euclidean Peter-Weyl outgoing constraint column.

The output schema matches the complete Gauss-state serialization used by the
Lorentzian epsilon aggregate and can be consumed directly by
bqg_constraint_master_assembler_gate.py.

The microscopic action is evaluated once before sparse pruning.  The payload
records the exact discarded-vector norm from that same action so any later
near-zero/master calculation can bound the numerical effect of the declared
prune threshold rather than treating the threshold itself as an error bar.
"""
from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW


def encode(state):
    rows=[]
    for (spins,Ks),amp in sorted(state.items(),key=lambda kv:repr(kv[0])):
        rows.append({"spins":[int(x) for x in spins],"K_labels":[int(x) for x in Ks],"amp":[float(complex(amp).real),float(complex(amp).imag)]})
    return rows


def norm2(state): return float(sum(abs(z)**2 for z in state.values()))


def difference(a,b):
    keys=set(a)|set(b);out={}
    for k in keys:
        z=complex(a.get(k,0j))-complex(b.get(k,0j))
        if z!=0j: out[k]=z
    return out


def run(node,input_index,jmax2=5,prune=1e-8):
    basis=PW.basis_full_jhalf()
    if not (0<=node<5): raise ValueError("node must be 0..4")
    if not (0<=input_index<len(basis)): raise ValueError("bad input index")
    initial=basis[input_index]
    raw=PW.apply_H_cached_state({initial:1+0j},node,jmax2)
    out=PW.prune_state(dict(raw),prune)
    discarded=difference(raw,out)
    raw_n2=norm2(raw);ret_n2=norm2(out);disc_n2=norm2(discarded)
    raw_norm=raw_n2**0.5;ret_norm=ret_n2**0.5;disc_norm=disc_n2**0.5
    pyth=abs(raw_n2-ret_n2-disc_n2)/max(raw_n2,1e-300)
    max_disc=max((abs(z) for z in discarded.values()),default=0.0)
    rows=encode(out); payload={
        "status":"complete Euclidean Peter-Weyl outgoing constraint column",
        "passed":bool(pyth<2e-12 and max_disc<=prune*(1+1e-12)),
        "family":"E","node":node,"input_index":input_index,"input_K_labels":list(initial[1]),
        "Jmax":jmax2/2,"prune_threshold":prune,
        "complete_gauss_outgoing_column":{"support":len(out),"norm":ret_norm,"basis":"Peter-Weyl Gauss basis (spins,K_labels)","state":rows},
        "pruning_audit":{
            "raw_support":len(raw),"retained_support":len(out),"discarded_support":len(discarded),
            "raw_norm":raw_norm,"retained_norm":ret_norm,"discarded_norm":disc_norm,
            "relative_discarded_norm":disc_norm/max(raw_norm,1e-300),
            "max_discarded_amplitude":max_disc,
            "pythagorean_relative_error":pyth,
            "interpretation":"discarded_norm is measured from the same unpruned microscopic action; no second H_E evaluation is used"
        },
        "column_sha256":hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(",",":")).encode()).hexdigest(),
        "claim_boundary":"Serialized microscopic H_E column plus measured pruning error norm only; no projector or physical observable is inferred."
    }
    return payload


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--node',type=int,required=True); ap.add_argument('--input-index',type=int,required=True); ap.add_argument('--jmax2',type=int,default=5); ap.add_argument('--prune',type=float,default=1e-8); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args()
    out=run(a.node,a.input_index,a.jmax2,a.prune); a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8'); print(json.dumps({k:v for k,v in out.items() if k!='complete_gauss_outgoing_column'},indent=2)); return 0 if out['passed'] else 1

if __name__=='__main__': raise SystemExit(main())
