#!/usr/bin/env python3
"""Real two-input linearity gate for frozen K5 Lorentzian KA at node1.
Independent DIRECT two-key C(K) action vs sum of TWO already completed full
complex outputs. No synthetic operator; exact same Jmax2, source and kernels.
"""
import json, math, sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import bqg_lorentzian_checkpointed_triple as Base
import bqg_node1_ka_partial_gram as Gram
import peter_weyl_lorentzian_logical_projection_gate as LP
import peter_weyl_lorentzian_ordered_triple_gate as RAW

ROOT=Path("checkpoint")
IDS=(0,12)
def norm(x):return math.sqrt(sum(abs(v)**2 for v in x.values()))
meta=json.loads((ROOT/"kb_000.json").read_text(encoding="utf-8"))
kb,_=Base.load_result(ROOT,"kb",(0,0,0),meta["freeze"])
rows=sorted(kb.items(),key=lambda kv:repr(kv[0]))
assert len(rows)==120 and meta["freeze"]["source"]==1
assert meta["freeze"]["triple"]==[0,2,3] and meta["freeze"]["JMAX2"]==7
mixed={rows[i][0]:rows[i][1] for i in IDS}
assert len(mixed)==2
saved={}
for idx in IDS:
    vec,one=Gram.vector(ROOT,idx)
    assert one["freeze"]==meta["freeze"]
    assert one["kb_sha256"]==meta["sha256"] and one["rows"]==1
    # Gram.vector returns a serialized 3-part key (spins,Kother,charges4),
    # while the genuine CK engine uses a FLAT 6-part covariant key.
    # The old comparison matched disjoint dictionaries and produced sqrt(2).
    native={ (k[0],k[1],*k[2]):v for k,v in vec.items() }
    assert all(len(k)==6 and len(k[0])==10 and len(k[1])==5 for k in native)
    assert len(native)==len(vec)
    saved[idx]=native
reference={}
for idx in IDS:
    for k,v in saved[idx].items():reference[k]=reference.get(k,0j)+v
LP.JMAX2=Base.JMAX2
restore,_=LP.install_sine_cached_stack()
try:
    direct,diag=RAW.KCOMP.C_K_component(mixed,1,0,0,0,Base.JMAX2)
finally:restore()
delta={k:direct.get(k,0j)-reference.get(k,0j) for k in set(direct)|set(reference)}
error=norm(delta)
scale=max(norm(direct),norm(reference),1e-30)
shared=set(direct)&set(reference)
overlap=sum(np.conjugate(direct[k])*reference[k] for k in shared)
report=dict(status="PASS" if error/scale<1e-8 else "FAIL",
    test="independent REAL C_K two-source coherent linearity",
    source=1,channel="000",input_indices=list(IDS),
    direct_support=len(direct),summed_separate_support=len(reference),
    direct_norm=norm(direct),reference_norm=norm(reference),
    absolute_difference=error,relative_difference=error/scale,
    exact_shared_keys=len(shared),
    complex_overlap=[float(overlap.real),float(overlap.imag)],
    overlap_relative_to_norms=float(abs(overlap)/max(norm(direct)*norm(reference),1e-30)),
    direct_diagnostics=diag,
    frozen_KB_sha256=meta["sha256"],
    important="Live 6-field CK keys are normalized against archived 3-field serialized keys. Two real inputs only; no full HDA.")
print("BQG_TWO_KEY_REAL_LINEARITY_START",flush=True)
print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)
print("BQG_TWO_KEY_REAL_LINEARITY_END",flush=True)
if report["status"]!="PASS":raise SystemExit(2)
