#!/usr/bin/env python3
"""V26: one-hit graph-sector obstruction to naive BQG refinement.

Read all 25 CANONICAL V25 K5 complex output columns from V26 ZIP,
verify the original backend, and recompute source/target graph-sector
Gram matrices. This is a local Euclidean test, NOT a physical projector.
"""
import argparse, collections, hashlib, json, math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

SHA256="596c4179cb2c8af20b41c18b3ab87b711c3850c5c5e05b1bc91a1f87c9f11cc4"
BLOB="e7e1e0705c4b833dcee167a2036ef5f8929a4388"

def sqrt_hermitian(a):
    w,u=np.linalg.eigh((a+a.conj().T)/2)
    if w[0]<-1e-9:raise ValueError("non-positive sector Gram")
    return (u*np.sqrt(np.maximum(w,0)))@u.conj().T

def fidelity(a,b):
    root=sqrt_hermitian(a)
    c=root@b@root
    w=np.linalg.eigvalsh((c+c.conj().T)/2)
    return float(np.sqrt(np.maximum(w,0)).sum())

def grams(root,s):
    # Orthogonality is by both graph quantum state and constraint node v.
    d=s+1
    section=collections.defaultdict(lambda:[{} for _ in range(d)])
    for i,k in enumerate(range(0,2*s+1,2)):
        for v in range(5):
            file=root/f"results/canonical_columns/j2_{s}_K_{k}_node_{v}.json"
            j=json.loads(file.read_text())
            if j["backend_sha256"]!=SHA256 or j["node"]!=v or j["jmax2"]!=5:
                raise ValueError(f"bad column identity {file}")
            if j["input"]!=[[s]*10,[k,0,0,0,0]]:
                raise ValueError(f"bad source state {file}")
            for r in j["rows"]:
                state=(tuple(r["spin"]),tuple(r["K"]))
                graph=tuple(t for t,spin in enumerate(state[0]) if spin==0)
                q=(v,state)
                if q in section[graph][i]:raise ValueError("duplicate amplitude")
                a=complex(float.fromhex(r["real"]),float.fromhex(r["imag"]))
                if not (np.isfinite(a.real) and np.isfinite(a.imag)):
                    raise ValueError("nonfinite amplitude")
                section[graph][i][q]=a
    G={}
    for graph,cols in section.items():
        m=np.zeros((d,d),complex)
        for i in range(d):
            for j in range(i,d):
                val=sum(a.conjugate()*cols[j].get(q,0) for q,a in cols[i].items())
                m[i,j]=val;m[j,i]=val.conjugate()
        G[graph]=m
    return G

def audit(root):
    root=Path(root)
    raw=(root/"upstream/k5_peter_weyl_safe_hda_column.py").read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if hashlib.sha256(raw).hexdigest()!=SHA256 or blob!=BLOB:
        raise ValueError("wrong original K5 backend")
    mats=np.load(root/"results/V25_MATRICES.npz")
    W=mats["W"];J=np.load(root/"results/V25_BILINEAR.npz")["J"]
    assert np.linalg.norm(J.conj().T@J-np.eye(2))<1e-11
    f=grams(root,1);c=grams(root,2)
    A=sum(f.values(),np.zeros((2,2),complex))
    Mc=sum(c.values(),np.zeros((3,3),complex))
    B=W.conj().T@Mc@W
    assert np.linalg.norm(A-mats["Mhalf"])<1e-10
    assert np.linalg.norm(Mc-mats["Mone"])<1e-10
    missing=sum((v for k,v in f.items() if k),np.zeros((2,2),complex))
    eig=np.linalg.eigvalsh(missing)
    def coarse(graph):
        return W.conj().T@c[graph]@W if graph in c else np.zeros((2,2),complex)
    sector_fid=sum(fidelity(f.get(k,np.zeros((2,2))),coarse(k)) for k in f.keys()|c.keys())
    total_fid=fidelity(A,B)
    trA=float(np.trace(A).real);trB=float(np.trace(B).real)
    def best(fid):
        t=fid/trA
        return dict(amplitude_scale=t,
                    minimum_frobenius_residual=math.sqrt(max(0,trB-fid*fid/trA)))
    o=minimize_scalar(lambda s:np.linalg.norm(B-s*A,2),
                      method="bounded",bounds=(0,4),options={"xatol":1e-14})
    out={
        "scope":"canonical K5 finite 25-column j=1/2 to j=1 local Euclidean outputs",
        "source_graph_sectors":len(f),"coarse_graph_sectors":len(c),
        "source_absent_edge_norm2":float(np.trace(missing).real),
        "source_absent_edge_fraction":float(np.trace(missing).real/trA),
        "target_absent_edge_norm2":sum(float(np.trace(v).real) for k,v in c.items() if k),
        "source_missing_edge_Gram_eigenvalues":eig.tolist(),
        "graph_preserving_t_equals_one_operator_lower_bound":float(math.sqrt(max(eig))),
        "graph_preserving":best(sector_fid),
        "graph_mixing_allowed":best(total_fid),
        "Gram_best_energy_scale_s":float(o.x),
        "unavoidable_Gram_operator_defect":float(o.fun),
        "real_full_habitat_isometry":False,
        "full_physical_BQG_one_hit_residual":False,
        "physical_projector_interacting_TT_EH":False,
    }
    if len(f)!=47 or len(c)!=1 or out["target_absent_edge_norm2"]!=0:
        raise AssertionError("unexpected output graph structure")
    if out["unavoidable_Gram_operator_defect"]<=0.5:
        raise AssertionError("naive exact output isometry unexpectedly possible")
    return out

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--archive-root",type=Path,required=True)
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    d=audit(a.archive_root)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(d,indent=2,ensure_ascii=False))
