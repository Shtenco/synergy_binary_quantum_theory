#!/usr/bin/env python3
"""V18 parity-resolved K5 Schur/Feshbach finite-shell research audit.
Expected inputs: V16 32x32 Gram NPZ and V17 36x36 NPZ with 'states',
plus the V17 metadata JSON. See README_RU.md.
"""
import argparse, ast, json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

MIRROR_SHA="1353c9fd937d5dbc83615fd041dc8df5da3895a2"
OFFICIAL_SHA="e7e1e0705c4b833dcee167a2036ef5f8929a4388"

def audit(v16,v17,status):
    meta=json.loads(Path(status).read_text(encoding="utf-8"))
    if meta.get("backend_git_blob_sha") != MIRROR_SHA:
        raise ValueError("Unrecognized K5 research backend")
    if meta.get("canonical_sha_confirmed"):
        raise ValueError("Contradictory provenance")
    with np.load(v16,allow_pickle=False) as d: base=d["M"].copy()
    with np.load(v17,allow_pickle=False) as d:
        M=d["M"].copy()
        keys=[ast.literal_eval(str(s)) for s in d["states"]]
    if M.shape != (36,36) or base.shape != (32,32) or len(keys)!=36:
        raise ValueError("Invalid finite input dimensions")
    if np.linalg.norm(M-M.conj().T)>1e-10:
        raise ValueError("Master is not Hermitian")
    even=[i for i,k in enumerate(keys) if sum(k[0])%2==0]
    odd=[i for i,k in enumerate(keys) if sum(k[0])%2==1]
    if even[:32]!=list(range(32)) or len(even)!=33 or len(odd)!=3:
        raise ValueError("Unexpected parity basis")
    E=M[np.ix_(even,even)];O=M[np.ix_(odd,odd)]
    if np.linalg.norm(M[np.ix_(even,odd)])>1e-10:
        raise ValueError("Parity selection broken")
    A=E[:32,:32];b=E[:32,32];c=float(E[32,32].real)
    if np.linalg.norm(A-base)>5e-12:
        raise ValueError("V16 and V17 Gram blocks disagree")
    a0=float(np.linalg.eigvalsh(A)[0])
    e0=float(np.linalg.eigvalsh(E)[0])
    o0=float(np.linalg.eigvalsh(O)[0])
    if not 0<e0<a0:raise ValueError("Non-positive or unordered finite spectrum")
    f=lambda z: float((c-z-np.vdot(b,np.linalg.solve(A-z*np.eye(32),b))).real)
    froot=float(brentq(f,0.0,a0-1e-8,xtol=1e-13))
    checks=[]
    for z in (0.,3.,5.):
        whole=np.linalg.inv(E-z*np.eye(33))[:32,:32]
        F=A-z*np.eye(32)-np.outer(b,b.conj())/(c-z)
        err=float(np.linalg.norm(whole-np.linalg.inv(F)))
        checks.append({"z":z,"error":err})
        if err>1e-11:raise AssertionError("Feshbach resolvent failed")
    if abs(froot-e0)>3e-10:raise AssertionError("Schur eigenvalue root failed")
    all_e=np.linalg.eigvalsh(M)
    if np.linalg.norm(np.sort(all_e)-np.sort(np.r_[np.linalg.eigvalsh(E),np.linalg.eigvalsh(O)]))>1e-10:
        raise AssertionError("Parity spectral union mismatch")
    r={
       "experiment":"finite noncanonical K5 parity/Schur audit",
       "canonical_backend_verified":False,
       "research_mirror_blob":MIRROR_SHA,"official_github_blob":OFFICIAL_SHA,
       "input_dimensions":{"original":32,"even":len(even),"odd":len(odd),"total":36},
       "original_min":a0,"even_min":e0,"odd_min":o0,"full_min":float(all_e[0]),
       "even_min_shift":a0-e0,
       "second_step_diagonal":c,
       "coupling_norm":float(np.linalg.norm(b)),
       "coupling_nonzero_entries":int((abs(b)>1e-12).sum()),
       "static_schur":float((c-np.vdot(b,np.linalg.solve(A,b))).real),
       "feshbach_root":froot,"feshbach_root_error":abs(froot-e0),
       "resolvent_checks":checks,"parity_cross_norm":float(np.linalg.norm(M[np.ix_(even,odd)])),
       "numeric_nullity":int(np.sum(all_e<=max(all_e[-1],1.)*1e-10)),
       "full_physical_projector":False,
       "lorentzian_HDA_verified":False,
       "claim_boundary":"Selected 36D Euclidean mirror only. No full shell or continuum."
    }
    return r

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--v16-npz",type=Path,required=True)
    p.add_argument("--v17-npz",type=Path,required=True)
    p.add_argument("--v17-json",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    d=audit(a.v16_npz,a.v17_npz,a.v17_json)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(d,indent=2,ensure_ascii=False))
