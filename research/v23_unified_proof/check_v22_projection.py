#!/usr/bin/env python3
"""BQG V23: evaluate genuine finite Gram evidence from V22, never infer full physics."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
STEPS=[0,1,4,8,16,24,32,41,49,65,69]
SOURCE_SHA="596c4179cb2c8af20b41c18b3ab87b711c3850c5c5e05b1bc91a1f87c9f11cc4"
def evaluate(root):
    root=Path(root)
    f=root/"V21_RESTRICTED_MASTER.npz"
    with np.load(f,allow_pickle=False) as z:M=z["M"]
    status=json.loads((root/"V22_STATUS.json").read_text())
    if M.shape!=(101,101) or status.get("canonical_original_sha256")!=SOURCE_SHA:
        raise ValueError("wrong V22 provenance or matrix")
    if status["complete_first_shell_sources"]!=69:raise ValueError("wrong V22 coverage")
    if np.linalg.norm(M-M.conj().T)>1e-10:raise ValueError("not Hermitian")
    eig=np.linalg.eigvalsh(M)
    if eig[0]<=0:raise ValueError("unexpected numerical nullity")
    if abs(eig[0]-status["restricted_Gram_min_eigenvalue"])>1e-7:
        raise ValueError("numerical data mismatch")
    history=[]
    for a,b in zip(STEPS,STEPS[1:]):
        A=M[:32+a,:32+a];B=M[:32+b,:32+b]
        I=np.eye(32+b,32+a,dtype=complex)
        R=B@I-I@A
        gap=min(np.linalg.eigvalsh(A)[0],np.linalg.eigvalsh(B)[0])
        if gap<=0:raise ValueError("nonpositive gap")
        rn=float(np.linalg.norm(R,2))
        history.append({"old_shell":a,"new_shell":b,"residual_norm":rn,
                        "minimum_positive_gap":float(gap),"epsilon":rn/gap})
    if np.linalg.norm(M[:32,32:])>1e-10:raise ValueError("parity block violated")
    result={"K5_backend_sha256":SOURCE_SHA,"v22_matrix_sha256":hashlib.sha256(f.read_bytes()).hexdigest(),
        "full_matrix_dim":101,"complete_first_shell_inputs":69,"required":720,
        "finite_master_rank":101,"minimum_eigenvalue":float(eig[0]),
        "maximum_eigenvalue":float(eig[-1]),"trivial_zero_projectors":True,
        "even_to_odd_master_norm":float(np.linalg.norm(M[:32,32:])),
        "nested_coordinate_inclusion_epsilon":history,
        "observed_finite_sum_epsilon":float(sum(x["epsilon"] for x in history)),
        "coordinates_are_true_BQG_refinement":False,
        "full_Lorentzian_constraints_verified":False,
        "sum_over_all_physical_refinements_proved":False,
        "nonzero_physical_anchor_proved":False,
        "physical_connected_history_proved":False,
        "interacting_physical_TT_pole_proved":False,
        "microscopic_Einstein_IR_proved":False,
        "note":"Finite 101D K5 Gram only. No physical BQG projector is inferred."}
    return result
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--v22-data",type=Path,required=True,help="Extracted V22 JSON/NPZ data directory")
    p.add_argument("--out",type=Path)
    a=p.parse_args();out=evaluate(a.v22_data);s=json.dumps(out,ensure_ascii=False,indent=2)
    if a.out:a.out.write_text(s+"\n")
    print(s)
