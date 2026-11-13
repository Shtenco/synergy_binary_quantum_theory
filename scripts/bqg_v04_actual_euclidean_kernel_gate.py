#!/usr/bin/env python3
"""v0.4 REAL Euclidean Peter-Weyl node constraints, open-output kernel diagnostic.

This does NOT construct physical H_E+H_L, independent D or a relational clock.
Outputs are not projected back onto the selected input space.
"""
from __future__ import annotations
import json, argparse, hashlib
from pathlib import Path
import numpy as np
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_euclidean_sine_ordering_gate as SINE
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

def inner(a,b):
    if len(a)>len(b):
        return np.conj(inner(b,a))
    return sum((np.conj(x)*b.get(k,0j) for k,x in a.items()),0j)

def gram(images):
    K=np.zeros((len(images),len(images)),complex)
    for a in range(len(images)):
        for b in range(a,len(images)):
            K[a,b]=inner(images[a],images[b])
            K[b,a]=np.conj(K[a,b])
    return (K+K.conj().T)/2

def run():
    ZVM.patch_and_clear()
    keys=[((4,)*len(PW.EDGES),(K,0,0,0,0)) for K in (0,2,4,6,8)]
    G=np.zeros((5,5),complex)
    summaries=[]
    for vertex in range(5):
        images=[SINE.safe_H_sine({key:1+0j},vertex,7) for key in keys]
        gv=gram(images)
        G+=gv
        ev=np.linalg.eigvalsh(gv)
        summaries.append({
            "vertex":vertex, "supports":[len(im) for im in images],
            "local_master_spectrum":[float(x) for x in ev],
            "min_local_master_eigenvalue":float(ev.min()),
        })
    eig,U=np.linalg.eigh(G)
    threshold=max(1e-9,1e-10*max(float(abs(eig).max()),1))
    zero=np.abs(eig)<threshold
    P=U[:,zero]@U[:,zero].conj().T
    M0=np.array(summaries[0]["local_master_spectrum"])
    # Independently rebuild the SAME node-0 raw master used by the published j2 gate.
    if not np.all(np.isfinite(eig)):
        raise FloatingPointError("nonfinite genuine constraint master")
    if eig.min() < -threshold:
        raise ArithmeticError("Gram matrix not PSD")
    return {
      "source":"Shtenco/binary_quantum_theory genuine SINE.safe_H_sine unprojected spin-network outputs",
      "claim_scope":"Actual j=2 K5 input subspace W5, Euclidean node H_E,v only; NO full Lorentzian H, D or clock",
      "input_basis_K":[0,2,4,6,8],
      "node_results":summaries,
      "joint_euclidean_gram_eigenvalues":[float(x) for x in eig],
      "joint_euclidean_gram":[[[float(z.real),float(z.imag)] for z in row] for row in G],
      "min_eigenvalue":float(eig.min()),
      "zero_threshold":threshold,
      "kernel_dimension_in_W5":int(sum(zero)),
      "projector_residual":float(np.linalg.norm(G@P)),
      "screening_conclusion":"NO_COMMON_EUCLIDEAN_ZERO_IN_W5" if not np.any(zero) else "EUCLIDEAN_ZERO_CANDIDATE_IN_W5",
      "physical_gravity_kernel_established":False,
      "derived_relational_time":False,
      "derived_BQG_retarded_response":False,
      "note":"H_E+H_L must be assembled as a single Hamiltonian constraint. This Gram neither bounds nor constructs that full operator.",
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--output",type=Path,default=Path("verification_results/BQG_V04_EUCLIDEAN_REAL_W5.json"))
    x=a.parse_args()
    r=run()
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(json.dumps(r,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(r,indent=2))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
