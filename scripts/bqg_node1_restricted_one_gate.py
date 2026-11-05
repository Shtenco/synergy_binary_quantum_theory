#!/usr/bin/env python3
"""True SU(2) held-out *remaining* column gate, no synthetic model amplitudes.

Phase SELECT: verify frozen KB; choose largest-weight uncomputed basis key and
write index to $GITHUB_ENV. Record original KB as SHA-verifiable log evidence
so future work survives GitHub Actions artifact expiration.
Phase VERIFY: SHA-check complete KA outgoing state; divide by original input
coefficient; lower-bound ||C(K) Q_remaining||; optionally compute G21 from
all 20 archived genuine columns + one new genuine column.

A basis Hilbert norm certificate requires the established orthonormality of
the exact source Peter-Weyl key system. This test does not certify full HDA.
"""
from __future__ import annotations

import argparse, base64, hashlib, json, os, sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import bqg_node1_ka_partial_gram as gram
import bqg_node1_residual_operator_bound_gate as residual

EXTRA=(59,105,113,112,57,52,24)
DONE=tuple(residual.DONE)
FROZEN=residual.KB_SHA
MEASURED_PARTIAL_NORM=residual.MEASURED_OUT_NORM
EPS=0.10
OLD_REMAINING_NORM=0.08986860634523096
OLD_THRESHOLD=EPS*MEASURED_PARTIAL_NORM/OLD_REMAINING_NORM

def kb_checked(root):
    stem=root/"kb_000"
    meta=json.loads(stem.with_suffix(".json").read_text(encoding="utf-8"))
    raw=stem.with_suffix(".npz").read_bytes()
    assert hashlib.sha256(raw).hexdigest()==meta["sha256"]==FROZEN
    assert meta["stage"]=="kb" and meta["channel"]==[0,0,0]
    fr=meta["freeze"]
    assert fr["source"]==1 and fr["column_index"]==0 and fr["triple"]==[0,2,3]
    with np.load(stem.with_suffix(".npz"),allow_pickle=False) as src:
        amp=np.asarray(src["amp"],np.complex128)
        spins=np.asarray(src["spins"])
        assert len(amp)==len(spins)==120
        assert np.isfinite(amp.real).all() and np.isfinite(amp.imag).all()
    return amp,spins,meta,raw

def pick(root):
    amp,spins,meta,raw=kb_checked(root)
    remaining=[idx for idx in range(120) if idx not in DONE]
    index=max(remaining,key=lambda i:(float(abs(amp[i])**2),-i))
    weight=float(abs(amp[index])**2/sum(abs(amp)**2))
    remaining_norm=float(np.linalg.norm(amp[remaining]))
    report={
        "status":"FROZEN_KB_INPUT_REMAINDER_KEY_SELECTED",
        "index":index,
        "source_node":1,"input_column":0,"triple":[0,2,3],"auxiliary_channel":"000",
        "relative_input_weight_squared":weight,
        "input_coefficient_real":float(amp[index].real),
        "input_coefficient_imag":float(amp[index].imag),
        "max_spin_of_input_edges":float(np.max(spins[index])/2),
        "measured_KB_SHA256":meta["sha256"],
        "num_remain":len(remaining),
        "remainder_input_norm":remaining_norm,
        "required_upper_bound_for_10pct_via_restricted_norm":OLD_THRESHOLD,
        "note":"Selection by greatest not-yet-computed |a_i|^2, not a post hoc amplitude fit."}
    os.environ.get("GITHUB_ENV") and Path(os.environ["GITHUB_ENV"]).open("a").write(f"BQG_HELDOUT_INDEX={index}\n")
    print("BQG_INPUT_PRECHECK",json.dumps(report,separators=(",",":")),flush=True)
    print("BQG_KB_PARENT_BEGIN",flush=True)
    print("BQG_KB_PARENT_METADATA",json.dumps(meta,separators=(",",":")),flush=True)
    print("BQG_KB_PARENT_BASE64",base64.b64encode(raw).decode("ascii"),flush=True)
    print("BQG_KB_PARENT_END",flush=True)
    return report

def reconstruct_previous(root,repo):
    for i in EXTRA:
        source=repo/"research"/"node1_fullstate_gate_2026_10_31"/"log_recovered_ka"/f"key_{i}.json"
        d=json.loads(source.read_text(encoding="utf-8"))
        meta=d["metadata"]
        raw=base64.b64decode(d["npz_base64"],validate=True)
        assert hashlib.sha256(raw).hexdigest()==meta["sha256"]
        assert meta["part_index"]==i and meta["parts"]==120 and meta["rows"]==1
        stem=root/f"ka_000_part{i:03d}_of120"
        stem.with_suffix(".npz").write_bytes(raw)
        stem.with_suffix(".json").write_text(json.dumps(meta),encoding="utf-8")

def verify(root,repo,index):
    amp,spins,kbmeta,raw=kb_checked(root)
    assert index not in DONE and 0<=index<120
    stem=root/f"ka_000_part{index:03d}_of120"
    meta=json.loads(stem.with_suffix(".json").read_text(encoding="utf-8"))
    buffer=stem.with_suffix(".npz").read_bytes()
    assert hashlib.sha256(buffer).hexdigest()==meta["sha256"]
    assert meta["kb_sha256"]==FROZEN and meta["freeze"]==kbmeta["freeze"]
    assert meta["channel"]==[0,0,0] and meta["stage"]=="ka_part"
    assert meta["part_index"]==index and meta["parts"]==120 and meta["rows"]==1
    outgoing,checked=gram.vector(root,index)
    assert outgoing and abs(amp[index])>0
    unweighted_norm=float(np.linalg.norm([x/amp[index] for x in outgoing.values()]))
    remaining=[k for k in range(120) if k not in DONE]
    assert index in remaining
    rem_prev=float(np.linalg.norm(amp[remaining]))
    assert abs(rem_prev-OLD_REMAINING_NORM)<1e-12
    rem_next=float(np.linalg.norm(amp[[k for k in remaining if k!=index]]))
    output={
        "status":"REAL_HELDOUT_Q_COLUMN_SHA_VERIFIED",
        "source_run":"current workflow","new_source_index":index,
        "new_weighted_column_support":len(outgoing),
        "new_weighted_column_norm":meta["norm"],
        "unweighted_CK_column_norm":unweighted_norm,
        "certified_Q_operator_norm_lower_bound":unweighted_norm,
        "bound_required_for_10pct_certificate_before":OLD_THRESHOLD,
        "restricted_norm_10pct_method_excluded":bool(unweighted_norm>OLD_THRESHOLD),
        "remaining_source_count_after":len(remaining)-1,
        "remaining_input_norm_after":rem_next,
        "status_scope":"One held-out exact finite-cutoff Peter-Weyl source; NOT physical full Lorentzian HDA",
        "important":"If threshold excluded, actual coherent remainder may still be small; this is NO GO only for operator-norm sufficient condition."
    }
    print("BQG_HELDOUT_COLUMN_RESULT",json.dumps(output,separators=(",",":")),flush=True)
    print("BQG_HELDOUT_NUMERICAL_EVIDENCE_BEGIN",index,flush=True)
    print("BQG_HELDOUT_METADATA",json.dumps(meta,separators=(",",":")),flush=True)
    print("BQG_HELDOUT_NPZ_BASE64",base64.b64encode(buffer).decode("ascii"),flush=True)
    print("BQG_HELDOUT_NUMERICAL_EVIDENCE_END",index,flush=True)
    reconstruct_previous(root,repo)
    selection=tuple(sorted(DONE+(index,)))
    G=gram.calc(root,study=selection)
    assert len(selection)==21 and G["completed_KA_input_keys"]==21
    output["G21_Genuine_Peter_Weyl"]=G
    print("BQG21_VERIFIED_RESULT_BEGIN",flush=True)
    print(json.dumps(output,indent=2,ensure_ascii=False),flush=True)
    print("BQG21_VERIFIED_RESULT_END",flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mode",choices=("select","verify"),required=True)
    p.add_argument("--checkpoint-dir",type=Path,required=True)
    p.add_argument("--repo-root",type=Path,default=Path("."))
    p.add_argument("--index",type=int)
    a=p.parse_args()
    if a.mode=="select":pick(a.checkpoint_dir)
    else:
        if a.index is None:raise SystemExit("VERIFY requires --index")
        verify(a.checkpoint_dir,a.repo_root,a.index)
if __name__=="__main__":main()
