#!/usr/bin/env python3
"""BQG residual bound gate: exact KB input residual + CONDITIONAL factorized operator bound.
No invented full-operator norms. Only report a certified relative error after
explicitly supplied upper bounds for the REAL same-cutoff V and HE_sine.
"""
import argparse, json, math, hashlib
from pathlib import Path
import numpy as np

DONE=(0,1,2,3,4,5,12,13,16,17,24,32,33,52,57,58,59,105,112,113)
KB_SHA="f486e2250d4a5678a34438467b9e338d730ebf6e1d1ec956d2de7f8c2521b293"
MEASURED_OUT_NORM=0.06524151188803094

def run(kb_path,meta_path,v_bound=None,he_bound=None,contractive=None):
    b=kb_path.read_bytes()
    metadata=json.loads(meta_path.read_text(encoding="utf-8"))
    assert hashlib.sha256(b).hexdigest()==KB_SHA==metadata["sha256"]
    assert metadata["stage"]=="kb" and metadata["channel"]==[0,0,0]
    assert metadata["freeze"]["source"]==1 and metadata["freeze"]["triple"]==[0,2,3]
    with np.load(kb_path,allow_pickle=False) as archive:
        amps=np.asarray(archive["amp"],np.complex128)
    assert len(amps)==120 and np.isfinite(amps.real).all() and np.isfinite(amps.imag).all()
    done=np.zeros(120,bool);done[list(DONE)]=True
    norm_full=float(np.linalg.norm(amps))
    remaining=float(np.linalg.norm(amps[~done]))
    covered=float(np.linalg.norm(amps[done])**2/norm_full**2)
    result={
      "status":"INPUT_RESIDUAL_EXACT__OUTPUT_CERTIFICATE_PENDING",
      "num_processed":int(done.sum()),"num_remaining":int((~done).sum()),
      "kb_norm":norm_full,"kb_residual_norm":remaining,
      "kb_input_weight_covered":covered,
      "actual_partial_output_norm":MEASURED_OUT_NORM,
      "inequality":"||delta_out|| <= ||C_a(K) Q_remaining|| * ||Q_remaining psi_KB||",
      "certified_CK_norm_upper":None,"certified_absolute_error":None,
      "certified_relative_to_partial_output":None,
      "conditions":[
        "V and H_E_sine upper bounds for identical K5 graph, spin cutoff and domains",
        "Holonomy h unitary or truncated representation contractive",
        "All intertwiners orthonormal and intermediate projection maps contractive",
        "The implemented CK operator matches K=[V,H_E_sine] and CK=K-hKh^-1 on common habitat",
        "Primitive charge leak audit unresolved and must NOT be hidden by this bound"
      ],
      "physical_HDA":"OPEN","physical_TT":"OPEN"}
    for error in (0.01,0.05,0.1):
        result[f"CK_operator_norm_needed_for_{int(error*100)}pct_output_accuracy"]=(error*MEASURED_OUT_NORM/remaining)
    if v_bound is not None or he_bound is not None:
        if v_bound is None or he_bound is None or contractive is not True:
            raise SystemExit("BLOCKED: require BOTH independently certified V/HE_sine bounds and --contractive-proven")
        if not all(math.isfinite(a) and a>=0 for a in (v_bound,he_bound)):
            raise SystemExit("BLOCKED: invalid operator norm bound")
        ck=4*v_bound*he_bound
        result["certified_CK_norm_upper"]=ck
        result["certified_absolute_error"]=ck*remaining
        result["certified_relative_to_partial_output"]=ck*remaining/MEASURED_OUT_NORM
        result["status"]="CONDITIONAL_NORM_BOUND_AFTER_EXPLICIT_USER_SUPPLIED_CERTIFICATES"
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--checkpoint-dir",type=Path,required=True)
    p.add_argument("--v-bound",type=float)
    p.add_argument("--he-bound",type=float)
    p.add_argument("--contractive-proven",action="store_true")
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    data=run(a.checkpoint_dir/"kb_000.npz",a.checkpoint_dir/"kb_000.json",a.v_bound,a.he_bound,a.contractive_proven)
    raw=json.dumps(data,ensure_ascii=False,indent=2)
    print(raw)
    if a.output:a.output.write_text(raw+"\n",encoding="utf-8")
if __name__=="__main__":main()
