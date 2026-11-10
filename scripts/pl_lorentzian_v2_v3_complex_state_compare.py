#!/usr/bin/env python3
"""Independent complex-state comparison of direct V2 and orbit-reconstructed V3.

Inputs are two separately collected worker roots (each with term JSON/NPZ pairs)
and the two collector state bundles. No phase alignment is applied to the primary
comparison: a global or term-dependent phase error is a real amplitude mismatch.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np

VERSION = "tetrahedral-charged-volume-v2"
MODES = ("forward", "adjoint")
FIELDS = ("forward", "adjoint", "S")
TOL = 1e-8


def load_sparse(path: Path):
    with np.load(path, allow_pickle=False) as d:
        spins=np.asarray(d["spins"]); Ks=np.asarray(d["Ks"]); amp=np.asarray(d["amp"])
    if spins.ndim!=2 or Ks.ndim!=2 or amp.ndim!=1 or spins.shape[1]!=32 or Ks.shape[1]!=16:
        raise ValueError(f"malformed sparse state arrays: {path}")
    if not (len(spins)==len(Ks)==len(amp)):
        raise ValueError(f"misaligned sparse state arrays: {path}")
    if not np.issubdtype(spins.dtype,np.integer) or not np.issubdtype(Ks.dtype,np.integer):
        raise ValueError(f"non-integer state labels: {path}")
    if not np.isfinite(amp.real).all() or not np.isfinite(amp.imag).all():
        raise ValueError(f"non-finite amplitude: {path}")
    out={}
    for s,k,a in zip(spins,Ks,amp):
        key=(tuple(map(int,s)),tuple(map(int,k)))
        if key in out: raise ValueError(f"duplicate state key {key} in {path}")
        out[key]=complex(a)
    return out


def load_bundle(path: Path):
    with np.load(path, allow_pickle=False) as d:
        out={}
        for name in FIELDS:
            sk=np.asarray(d[name+"_spins"]); kk=np.asarray(d[name+"_Ks"]); aa=np.asarray(d[name+"_amp"])
            if sk.ndim!=2 or kk.ndim!=2 or aa.ndim!=1 or sk.shape[1]!=32 or kk.shape[1]!=16:
                raise ValueError(f"malformed {name} bundle arrays")
            if not (len(sk)==len(kk)==len(aa)): raise ValueError(f"misaligned {name} bundle arrays")
            if not np.isfinite(aa.real).all() or not np.isfinite(aa.imag).all():
                raise ValueError(f"non-finite {name} bundle amplitudes")
            state={}
            for s,k,a in zip(sk,kk,aa):
                key=(tuple(map(int,s)),tuple(map(int,k)))
                if key in state: raise ValueError(f"duplicate {name} key {key}")
                state[key]=complex(a)
            out[name]=state
    return out


def compare(a,b):
    keys=set(a)|set(b)
    na=math.sqrt(sum(abs(a.get(k,0j))**2 for k in keys))
    nb=math.sqrt(sum(abs(b.get(k,0j))**2 for k in keys))
    diff=math.sqrt(sum(abs(a.get(k,0j)-b.get(k,0j))**2 for k in keys))
    # Diagnostic only: the primary error above deliberately preserves absolute phase.
    overlap=sum(a.get(k,0j).conjugate()*b.get(k,0j) for k in keys)
    phase=overlap/abs(overlap) if abs(overlap)>1e-300 else 1+0j
    aligned=math.sqrt(sum(abs(a.get(k,0j)-phase*b.get(k,0j))**2 for k in keys))
    return {
        "support_a":len(a),"support_b":len(b),"support_equal":set(a)==set(b),
        "norm_a":na,"norm_b":nb,
        "relative_complex_state_error":diff/max(na,1e-300),
        "max_absolute_amplitude_error":max((abs(a.get(k,0j)-b.get(k,0j)) for k in keys),default=0.0),
        "optimal_global_phase_diagnostic":[float(phase.real),float(phase.imag)],
        "phase_aligned_relative_error_diagnostic":aligned/max(na,1e-300),
    }


def one_root(root: Path, label: str, expected_counts: dict):
    rows=[]
    files=list(root.rglob("term_*.json"))
    if len(files)!=48: raise RuntimeError(f"{label}: expected exactly 48 term JSON files, got {len(files)}")
    seen=set()
    for jp in files:
        meta=json.loads(jp.read_text(encoding="utf-8"))
        mode=meta.get("mode"); idx=meta.get("index")
        if mode not in MODES or isinstance(idx,bool) or not isinstance(idx,int) or not 0<=idx<24:
            raise RuntimeError(f"{label}: invalid term identity in {jp}")
        ident=(mode,idx)
        if ident in seen: raise RuntimeError(f"{label}: duplicate term identity {ident}")
        seen.add(ident)
        if meta.get("operator_version")!=VERSION: raise RuntimeError(f"{label}: wrong operator version in {jp}")
        npz=jp.with_suffix(".npz")
        if not npz.is_file(): raise RuntimeError(f"{label}: missing paired NPZ {npz}")
        digest=hashlib.sha256(npz.read_bytes()).hexdigest()
        if meta.get("state_sha256") and digest!=meta["state_sha256"]:
            raise RuntimeError(f"{label}: state SHA mismatch in {jp}")
        state=load_sparse(npz)
        execution=meta.get("execution_version","direct-v2-default")
        if execution not in expected_counts: raise RuntimeError(f"{label}: unexpected execution provenance {execution}")
        rows.append({"mode":mode,"index":idx,"execution_version":execution,"state":state})
    expected={(m,i) for m in MODES for i in range(24)}
    if seen!=expected: raise RuntimeError(f"{label}: incomplete term orbit; missing={sorted(expected-seen)}")
    counts={}
    for row in rows: counts[row["execution_version"]]=counts.get(row["execution_version"],0)+1
    if counts!=expected_counts: raise RuntimeError(f"{label}: execution provenance counts {counts}, expected {expected_counts}")
    return {(r["mode"],r["index"]):r for r in rows},counts


def run(v2_root: Path,v3_root: Path,v2_bundle: Path,v3_bundle: Path):
    v2,c2=one_root(v2_root,"V2",{"direct-v2-default":48})
    v3,c3=one_root(v3_root,"V3",{"direct-corrected-v2":12,"pairing-stabilizer-reconstructed-v3":36})
    term_rows=[]
    for ident in sorted(v2):
        d=compare(v2[ident]["state"],v3[ident]["state"])
        d.update({"mode":ident[0],"index":ident[1],"v3_execution_version":v3[ident]["execution_version"]})
        term_rows.append(d)
    b2=load_bundle(v2_bundle); b3=load_bundle(v3_bundle)
    bundle_rows={name:compare(b2[name],b3[name]) for name in FIELDS}
    max_term=max((r["relative_complex_state_error"] for r in term_rows),default=math.inf)
    max_bundle=max(r["relative_complex_state_error"] for r in bundle_rows.values())
    checks={
        "V2_has_48_direct_terms":c2=={"direct-v2-default":48},
        "V3_has_12_direct_and_36_reconstructed":c3=={"direct-corrected-v2":12,"pairing-stabilizer-reconstructed-v3":36},
        "all_48_term_supports_identical":all(r["support_equal"] for r in term_rows),
        "all_48_complex_term_states_agree":max_term<TOL,
        "forward_adjoint_and_S_bundle_supports_identical":all(r["support_equal"] for r in bundle_rows.values()),
        "forward_adjoint_and_S_complex_states_agree":max_bundle<TOL,
    }
    return {
        "status":"independent V2-direct versus V3-orbit complex-state comparison",
        "passed":all(checks.values()),
        "science_status":"AMPLITUDE_IMPLEMENTATION_CROSSCHECK_NOT_PHYSICAL_GR_PROOF",
        "tolerance":TOL,"checks":checks,
        "V2_execution_counts":c2,"V3_execution_counts":c3,
        "max_term_relative_complex_state_error":max_term,
        "max_bundle_relative_complex_state_error":max_bundle,
        "term_comparisons":term_rows,"bundle_comparisons":bundle_rows,
        "phase_policy":"No phase alignment is used for PASS. Phase-aligned errors are diagnostic only."
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--v2-root",type=Path,required=True)
    p.add_argument("--v3-root",type=Path,required=True)
    p.add_argument("--v2-bundle",type=Path,required=True)
    p.add_argument("--v3-bundle",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    out=run(a.v2_root,a.v3_root,a.v2_bundle,a.v3_bundle)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\\n",encoding="utf-8")
    print(json.dumps(out,indent=2,sort_keys=True))
    return 0 if out["passed"] else 1

if __name__=="__main__": raise SystemExit(main())
