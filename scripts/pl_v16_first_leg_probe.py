#!/usr/bin/env python3
"""True PL 16-cell Lorentzian first auxiliary leg; production backend only.

This tool computes one first-leg C(V) or C(K) amplitude vector. It does not
certify a completed Lorentzian pair, physical HDA or an exact continuum limit.
"""
from __future__ import annotations
import argparse, hashlib, json, math, os, time, traceback
from pathlib import Path
import numpy as np

def identify(n):
    if n not in (0,1,4,5,16,18,20,22,24,26,28,30):
        raise ValueError("NOT_A_CANONICAL_FIRST_LEG_INDEX")
    mode = "adjoint" if n >= 16 else "forward"
    order = "q" if (n//8)%2 else "p"
    i,j,k = (n//4)%2,(n//2)%2,n%2
    return mode,order,i,j,k

def save_charged(path, state, ne, nv):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    rows=sorted(state.items())
    sp=np.asarray([key[0] for key,z in rows],dtype=np.int16).reshape(-1,ne)
    ko=np.asarray([key[1] for key,z in rows],dtype=np.int16).reshape(-1,nv)
    charges=np.asarray([[key[2],key[3],key[4],key[5]] for key,z in rows],dtype=np.int16).reshape(-1,4)
    amp=np.asarray([z for key,z in rows],dtype=np.complex128)
    if not np.isfinite(amp).all():raise ValueError("NONFINITE_CHARGED_AMPLITUDE")
    np.savez_compressed(path,spins=sp,K_other=ko,charges=charges,amp=amp)
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--index",required=True,type=int)
    ap.add_argument("--output-prefix",required=True,type=Path)
    args=ap.parse_args()
    start=time.monotonic()
    outprefix=args.output_prefix
    outprefix.parent.mkdir(parents=True,exist_ok=True)
    try:
        mode,order,i,j,k=identify(args.index)
        from pl_dual_complex import DualComplex,seed_16cell_boundary
        from pl_peter_weyl_euclidean import PLPeterWeylEuclidean
        from pl_lorentzian_physical_pair_worker_v4 import pair_indices
        from pl_covariant_backend import install_pl_graph
        from tetrahedral_volume_backend import install_tetrahedral_volume_backend
        import peter_weyl_lorentzian_logical_projection_gate as LP
        D=DualComplex(seed_16cell_boundary());G=PLPeterWeylEuclidean(D)
        p,q,ip,iq=pair_indices(D,0,0,0)
        perm=p if order=="p" else q
        a,b,c=(D.neighbor[(0,t)] for t in perm)
        seed=((1,)*len(G.EDGES),(0,)*D.n_tets)
        with install_tetrahedral_volume_backend():
            with install_pl_graph(G):
                oldJ=LP.JMAX2;LP.JMAX2=7
                restore,caches=LP.install_sine_cached_stack()
                try:
                    psi=LP.CV.gauss_to_covariant({seed:1.+0j},0)
                    if mode=="forward":
                        state,leak=LP.RAW.COMP.C_volume_component(psi,0,c,k,i,7)
                        diag={"CV_complete_basis_leakage":float(leak)}
                        signature=["CV",int(c),int(k),int(i)]
                    else:
                        state,detail=LP.RAW.KCOMP.C_K_component(psi,0,a,j,i,7)
                        diag={str(name):float(val) for name,val in detail.items()
                              if isinstance(val,(int,float,np.integer,np.floating))}
                        signature=["CK",int(a),int(j),int(i)]
                finally:
                    restore();LP.JMAX2=oldJ
        if not all(math.isfinite(v) for v in diag.values()):
            raise ValueError("NONFINITE_LEAKAGE")
        digest=save_charged(outprefix.with_suffix(".npz"),state,len(G.EDGES),D.n_tets)
        result={"status":"REAL_FIRST_LORENTZIAN_AUX_LEG_ONLY",
                "passed":True,"full_pair_completed":False,"HDA_proved":False,
                "source_commit":os.getenv("GITHUB_SHA","UNKNOWN"),
                "graph":"PL-16cell","source_node":0,"omit":0,"cycle":0,
                "index":args.index,"mode":mode,"ordering":order,
                "aux_indices":[i,j,k],"operator_signature":signature,
                "Jmax2":7,"state_support":len(state),
                "state_norm":float(math.sqrt(sum(abs(z)**2 for z in state.values()))),
                "diagnostics":diag,"npz_sha256":digest,
                "runtime_seconds":time.monotonic()-start}
        exit_code=0
    except Exception as e:
        result={"status":"FIRST_LEG_COMPUTE_FAILED_NOT_SCIENCE_DATA",
                "passed":False,"index":args.index,
                "error_type":type(e).__name__,"error":str(e),
                "traceback":traceback.format_exc(),"runtime_seconds":time.monotonic()-start}
        exit_code=2
    outprefix.with_suffix(".json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2),flush=True)
    return exit_code
if __name__=="__main__":
    raise SystemExit(main())
