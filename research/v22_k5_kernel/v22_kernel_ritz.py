#!/usr/bin/env python3
"""V22 canonical K5 nested-kernel evidence — finite Ritz diagnostic, NOT P_phys.

Supply the BQG_V21 folder extracted from the published V22 archive.
Build Gram again from canonical original and complete shell input columns,
then report variational Ritz minima and eigenvector overlaps.
"""
import argparse, hashlib, importlib.util, json, sys
from pathlib import Path
import numpy as np

EXPECTED_SHA256 = "596c4179cb2c8af20b41c18b3ab87b711c3850c5c5e05b1bc91a1f87c9f11cc4"
EXPECTED_GIT_BLOB = "e7e1e0705c4b833dcee167a2036ef5f8929a4388"

def assemble(root):
    root = Path(root).resolve()
    canon = root / "backend" / "k5_peter_weyl_safe_hda_column.py"
    b=canon.read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
    if hashlib.sha256(b).hexdigest()!=EXPECTED_SHA256 or blob!=EXPECTED_GIT_BLOB:
        raise ValueError("Canonical K5 backend SHA mismatch")
    sys.path.insert(0,str(root))
    import v21_shell720 as v
    import v21_analyze as a
    if Path(v.CANON).resolve()!=canon:
        raise ValueError("Unexpected active K5 backend path")
    meta=a.run()  # independently re-assembles all complete input images
    with np.load(root/"V21_RESTRICTED_MASTER.npz",allow_pickle=False) as f:
        M=f["M"].copy()
    if M.shape!=(32+len(meta["fully_computed_first_shell_indices"]),)*2:
        raise ValueError("Incomplete/mismatched column coverage")
    return meta,M

def analyze(meta,M,steps=(1,4,8,16,24,32,41,49,65,69)):
    n=M.shape[0]-32
    if M.shape!=(n+32,n+32) or n != len(meta["fully_computed_first_shell_indices"]):
        raise ValueError("Coverage mismatch")
    if not np.all(np.isfinite(M)) or np.linalg.norm(M-M.conj().T)>1e-9:
        raise ValueError("Non-Hermitian/infinite Gram")
    cross=float(np.linalg.norm(M[:32,32:]))
    if cross>1e-10:
        raise ValueError("Parity forbidden block nonzero")
    G=M[32:,32:];hist=[];prior=None
    for k in sorted(set(x for x in steps if 1<=x<=n)|{n}):
        w,U=np.linalg.eigh(G[:k,:k])
        u=U[:,0]
        if w[0]<-1e-9:raise ValueError("Gram is not PSD")
        hist.append({
            "shell_inputs":k,"lowest_Gram_eigenvalue":float(w[0]),
            "constraint_norm":float(np.sqrt(max(float(w[0]),0))),
            "Ritz_residual":float(np.linalg.norm(G[:k,:k]@u-w[0]*u)),
            "overlap_with_previous":float(abs(np.vdot(prior,u[:len(prior)]))) if prior is not None else None,
            "participation_ratio":float(1/np.sum(abs(u)**4))
        })
        prior=u
    if any(hist[i]["lowest_Gram_eigenvalue"]>hist[i-1]["lowest_Gram_eigenvalue"]+1e-10
           for i in range(1,len(hist))):
        raise ValueError("Ritz variational monotonicity violated")
    w=np.linalg.eigvalsh(M)
    rank=int(np.count_nonzero(w>1e-10*max(float(w[-1]),1.)))
    return {
      "initial_sources":32,"complete_first_shell_sources":n,"required_first_shell_sources":720,
      "verified_complete_first_shell_columns":5*n,
      "restricted_rank":rank,"restricted_nullity":len(w)-rank,
      "restricted_min_eigenvalue":float(w[0]),"even_min_eigenvalue":float(np.linalg.eigvalsh(M[:32,:32])[0]),
      "parity_cross_block_norm":cross,
      "nested_ritz_history":hist,
      "one_H_first_shell_safe_spin_wall_2j":4,
      "next_layer_generic_safe_spin_wall_2j":6,
      "physically_normalizable_kernel_found":False,
      "full_shell_kernel_convergence_proven":False,
      "warning":"Stable finite Ritz vector != convergent full physical projector."
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--v21-root",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    meta,M=assemble(args.v21_root)
    d=analyze(meta,M)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(d,indent=2,ensure_ascii=False))

if __name__=="__main__":main()
