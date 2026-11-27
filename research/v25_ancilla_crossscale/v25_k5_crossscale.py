#!/usr/bin/env python3
"""BQG V25: canonical K5 j=1/2 -> j=1 finite local cross-scale audit.

This builds true q=2 old+ancilla SU(2) singlet bilinear fusion, an entangled
logical isometry, and exact one-hit Euclidean K5 Gram operators including all
external output keys. This is NOT an isometry of the full graph-changing
Lorentzian BQG habitat; output-space embeddings remain to be constructed.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

HERE=Path(__file__).resolve().parent
SCRIPTS=HERE.parents[1]/"scripts"
sys.path.insert(0,str(SCRIPTS))
import k5_peter_weyl_safe_hda_column as K5

ORIGINAL_BLOB="e7e1e0705c4b833dcee167a2036ef5f8929a4388"
ORIGINAL_SHA256="596c4179cb2c8af20b41c18b3ab87b711c3850c5c5e05b1bc91a1f87c9f11cc4"
CHAR={(1,1,1,1):2,(2,1,1):0,(2,2):2,(3,1):-1,(4,):0}

def check_canonical():
    raw=Path(K5.__file__).read_bytes()
    digest=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\x00"+raw).hexdigest()
    if digest!=ORIGINAL_BLOB or hashlib.sha256(raw).hexdigest()!=ORIGINAL_SHA256:
        raise RuntimeError(f"wrong upstream K5 blob={digest}")

def cycles(p):
    seen=set();sizes=[]
    for a in range(4):
        if a in seen:continue
        x=a;n=0
        while x not in seen:
            seen.add(x);x=p[x];n+=1
        sizes.append(n)
    return tuple(sorted(sizes,reverse=True))

def basis(s):
    return [K5.intertwiner_tensor_cached((s,)*4,k) for k in range(0,2*s+1,2)]

def representation(bas,p):
    return np.array([[np.vdot(a,np.transpose(b,p)) for b in bas] for a in bas])

def p22(bas):
    P=sum((CHAR[cycles(p)]*representation(bas,p) for p in itertools.permutations(range(4))),
          np.zeros((len(bas),len(bas)),complex))*2/24
    return (P+P.conj().T)/2

def microscopic():
    half=basis(1);one=basis(2)
    W=np.column_stack([[0,1,0],[2/3,0,-math.sqrt(5)/3]]).astype(complex)
    B=np.zeros((3,2,2),complex)
    for i,m in enumerate((2,0,-2)):
        for j,a in enumerate((1,-1)):
            for k,b in enumerate((1,-1)):
                B[i,j,k]=K5.cg2(1,1,2,a,b,m)
    F=np.zeros((3,4),complex)
    for i,old in enumerate(half):
        for k,anc in enumerate(half):
            state=np.einsum("Aae,Bbf,Ccg,Ddh,abcd,efgh->ABCD",
                            B,B,B,B,old,anc,optimize=True)
            F[:,2*i+k]=[np.vdot(x,state) for x in one]
    F=p22(one)@F
    J=math.sqrt(2)*F.conj().T@W
    covariance=max(np.linalg.norm(np.kron(representation(half,p),representation(half,p))@J
                                  -J@representation(half,p))
                   for p in itertools.permutations(range(4)))
    checks={
        "J_isometry_residual":float(np.linalg.norm(J.conj().T@J-np.eye(2))),
        "J_S4_covariance_residual":float(covariance),
        "FJ_minus_W_over_sqrt2":float(np.linalg.norm(F@J-W/math.sqrt(2))),
        "F_singular_values":np.linalg.svd(F,compute_uv=False).tolist(),
        "full_BQG_isometry_constructed":False}
    assert max(checks["J_isometry_residual"],checks["J_S4_covariance_residual"],
               checks["FJ_minus_W_over_sqrt2"])<1e-10
    return W,F,J,checks

def master(s2,output):
    cols=[]
    ks=list(range(0,2*s2+1,2))
    for k in ks:
        src=((s2,)*10,(k,0,0,0,0))
        for v in range(5):
            image=K5.prune_state(K5.apply_H_cached_state({src:1+0j},v,5),1e-10)
            cols.append(image)
            print(f"j2={s2} K={k} node={v} nnz={len(image)}",flush=True)
            file=output/f"source_spin2_{s2}_K_{k}_node_{v}.json"
            file.write_text(json.dumps({
                "canonical_backend_sha256":ORIGINAL_SHA256,"node":v,
                "input":[list(src[0]),list(src[1])],
                "rows":[{"key":[list(q[0]),list(q[1])],
                         "real":float(a.real).hex(),"imag":float(a.imag).hex()}
                        for q,a in sorted(image.items(),key=lambda x:repr(x[0]))]},
                ensure_ascii=False)+"\n")
    M=np.zeros((len(ks),len(ks)),complex)
    for i in range(len(ks)):
        for j in range(i,len(ks)):
            val=0j
            for v in range(5):
                a=cols[5*i+v];b=cols[5*j+v]
                val+=sum(x.conjugate()*b.get(q,0j) for q,x in a.items())
            M[i,j]=val
            M[j,i]=val.conjugate()
    return M

def twirl(M,bas):
    out=sum((representation(bas,p).conj().T@M@representation(bas,p)
            for p in itertools.permutations(range(4))),np.zeros_like(M))/24
    return (out+out.conj().T)/2

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    check_canonical()
    a.out_dir.mkdir(parents=True,exist_ok=True)
    W,F,J,checks=microscopic()
    A=master(1,a.out_dir);B=master(2,a.out_dir)
    P=W@W.conj().T
    L=(np.eye(3)-P)@B@W
    D=W.conj().T@B@W-A
    R=B@W-W@A
    assert np.linalg.norm(R-L-W@D)<1e-11
    opt=minimize_scalar(lambda s:np.linalg.norm(B@W-s*W@A,2),
                        method="bounded",bounds=(.01,5.))
    At=twirl(A,basis(1));Bt=twirl(B,basis(2))
    st=float((W.conj().T@Bt@W)[0,0].real/At[0,0].real)
    result={
        "canonical_backend_git_blob":ORIGINAL_BLOB,
        "status":"finite j=1/2 to j=1 local Euclidean K5 cross-spin audit",
        "canonical_columns":25,
        "microscopic":checks,
        "Mhalf_eigenvalues":np.linalg.eigvalsh(A).tolist(),
        "Mone_eigenvalues":np.linalg.eigvalsh(B).tolist(),
        "R_norm":float(np.linalg.norm(R,2)),
        "L_leakage_norm":float(np.linalg.norm(L,2)),
        "D_internal_norm":float(np.linalg.norm(D,2)),
        "R_split_error":float(np.linalg.norm(R-L-W@D)),
        "optimal_positive_energy_scale":float(opt.x),
        "residual_after_optimal_energy_rescaling":float(opt.fun),
        "artificial_S4_twirl_scale":st,
        "artificial_S4_twirl_residual":float(np.linalg.norm(Bt@W-st*W@At,2)),
        "physical_projector":False,
        "full_Lorentzian_constraint_family":False,
        "true_full_habitat_intertwiner":False,
        "interacting_TT_and_Einstein_IR":False}
    (a.out_dir/"V25_RESULT.json").write_text(json.dumps(result,indent=2)+"\n")
    np.savez_compressed(a.out_dir/"V25_MATRICES.npz",W=W,F=F,J=J,Mhalf=A,Mone=B,R=R,L=L,D=D)
    print(json.dumps(result,indent=2))
if __name__=="__main__":main()
