#!/usr/bin/env python3
"""Reusable BQG Euclidean/Lorentzian/full-Dirac constraint-master assembler.

Expensive microscopic actions are serialized as sparse outgoing columns. Given
one column for every declared domain basis vector and every separately labelled
normal constraint, this module assembles

    M_H(lambda) = M_EE + lambda M_EL + lambda^2 M_LL.

If explicit target tangential/diffeomorphism columns are supplied, it also adds

    M_D = sum_I D_I^dagger D_I,
    M_full = M_H + M_D.

A spectral physical projector is emitted only when ALL conditions hold:

1. the serialized domain is declared complete for the finite regulated habitat;
2. E/L/D columns act in the declared common habitat and master convention;
3. a matching BQG_QUANTUM_HDA_RESIDUAL_CERTIFICATE_V1 proves the actual
   graph-changing quantum residual [H_k,H_l]-D_target on that same packet.

Explicit D_target columns are necessary when those tangential constraints belong
in M_full, but they are NOT by themselves a proof of HDA closure.

Restricted boundary/Krylov domains, legacy target-only certificates, and
HDA-uncertified complete masters are diagnostics and never emit P_phys.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

LEGACY_DTARGET_SCHEMA = "BQG_DTARGET_HDA_CERTIFICATE_V1"
QUANTUM_HDA_SCHEMA = "BQG_QUANTUM_HDA_RESIDUAL_CERTIFICATE_V1"


def sparse_inner(a: dict, b: dict) -> complex:
    if len(a) > len(b):
        return np.conj(sparse_inner(b, a))
    return sum(np.conj(x) * b.get(k, 0j) for k, x in a.items())


def gram(images: list[dict]) -> np.ndarray:
    n=len(images); G=np.zeros((n,n),complex)
    for i in range(n):
        for j in range(i,n):
            z=sparse_inner(images[i],images[j]); G[i,j]=z; G[j,i]=np.conj(z)
    return G


def cross_gram(left: list[dict], right: list[dict]) -> np.ndarray:
    if len(left)!=len(right): raise ValueError("left/right domain dimensions differ")
    n=len(left); G=np.zeros((n,n),complex)
    for i in range(n):
        for j in range(n): G[i,j]=sparse_inner(left[i],right[j])
    return G


def spectral_audit(M: np.ndarray, rtol=1e-10):
    M=0.5*(M+M.conj().T)
    ev,U=np.linalg.eigh(M)
    scale=max(float(np.max(np.abs(ev))),1.0); thr=rtol*scale
    zero=np.abs(ev)<=thr; pos=ev[ev>thr]; Q=U[:,zero]
    return {"matrix":M,"eigenvalues":ev,"rank_tolerance":thr,
            "rank":int(np.sum(ev>thr)),"nullity":int(np.sum(zero)),
            "smallest_positive":float(np.min(pos)) if pos.size else None,
            "condition_number_on_support":float(np.max(pos)/np.min(pos)) if pos.size else None,
            "Q0":Q,"P0":Q@Q.conj().T}


def assemble(node_E, node_L, lam):
    nodes=sorted(set(node_E)|set(node_L))
    if set(node_E)!=set(node_L): raise ValueError("E/L node labels do not match")
    if not nodes: raise ValueError("empty normal constraint family")
    n=len(node_E[nodes[0]])
    if any(len(node_E[v])!=n or len(node_L[v])!=n for v in nodes):
        raise ValueError("each E/L node must contain every domain column")
    MEE=np.zeros((n,n),complex); MLL=np.zeros_like(MEE); MEL=np.zeros_like(MEE); rows=[]
    for v in nodes:
        EE=gram(node_E[v]); LL=gram(node_L[v]); X=cross_gram(node_E[v],node_L[v]); mixed=X+X.conj().T
        MEE+=EE; MLL+=LL; MEL+=mixed
        rows.append({"node":v,"EE_trace":float(np.trace(EE).real),"LL_trace":float(np.trace(LL).real),
                     "EL_cross_frobenius_norm":float(np.linalg.norm(X)),"mixed_frobenius_norm":float(np.linalg.norm(mixed))})
    M=MEE+lam*MEL+lam*lam*MLL
    return MEE,MEL,MLL,0.5*(M+M.conj().T),rows


def assemble_dtarget(d_images: dict[str,list[dict]], dim: int):
    MD=np.zeros((dim,dim),complex); rows=[]
    for label in sorted(d_images):
        imgs=d_images[label]
        if len(imgs)!=dim or any(x is None for x in imgs):
            raise ValueError(f"D_target {label} does not contain every domain column")
        G=gram(imgs); MD+=G
        rows.append({"constraint":label,"trace":float(np.trace(G).real),"frobenius_norm":float(np.linalg.norm(G))})
    return 0.5*(MD+MD.conj().T),rows


def hash_arrays(*arrs):
    h=hashlib.sha256()
    for A in arrs:
        x=np.ascontiguousarray(A); h.update(str(x.shape).encode()); h.update(x.view(np.uint8))
    return h.hexdigest()


def decode_state_rows(rows):
    state={}
    for r in rows:
        if "K_labels" in r:
            key=(tuple(int(x) for x in r["spins"]),tuple(int(x) for x in r["K_labels"]))
        elif "Kother" in r:
            key=(tuple(int(x) for x in r["spins"]),tuple(int(x) for x in r["Kother"]),
                 int(r["J2"]),int(r["M2"]),int(r["K12"]),int(r["K34"]))
        elif "basis_index" in r:
            key=(int(r["basis_index"]),)
        else:
            raise ValueError("unrecognized sparse-state row schema")
        z=complex(float(r["amp"][0]),float(r["amp"][1])); state[key]=state.get(key,0j)+z
    return state


def load_column(path: Path):
    data=json.loads(path.read_text(encoding="utf-8"))
    if "complete_gauss_outgoing_column" in data:
        rows=data["complete_gauss_outgoing_column"]["state"]
    elif "state" in data:
        rows=data["state"]
    else:
        raise ValueError(f"{path}: no complete_gauss_outgoing_column/state")
    return decode_state_rows(rows),data


def audit_hda_certificate(cert: dict, manifest: dict, master_hash: str):
    schema=cert.get("schema")
    if schema == LEGACY_DTARGET_SCHEMA:
        return {
            "schema": schema,
            "schema_ok": False,
            "legacy_target_only_schema": True,
            "quantum_habitat_residual_certified": False,
            "valid_for_this_master": False,
            "reason": "legacy Dtarget target-side certificate cannot authorize P_phys; actual HH-Dtarget residual certificate required",
        }

    schema_ok=(schema==QUANTUM_HDA_SCHEMA)
    quantum_ok=bool(cert.get("quantum_habitat_residual_certified",False))
    pre=cert.get("preconditions",{})
    no_fit=bool(pre.get("lambda_was_not_fitted",False)) and not bool(
        cert.get("norm_squared_polynomial",{}).get("lambda_fit_used_for_certification",True)
    )
    cert_hash=cert.get("hashes",{})
    expected={
        "habitat_hash": str(manifest.get("habitat_hash","")),
        "domain_hash": str(manifest.get("domain_hash","")),
        "constraint_packet_hash": str(master_hash),
        "convention_hash": str(manifest.get("convention_hash","")),
    }
    expected_present=all(bool(v) for v in expected.values())
    hash_match={name: (str(cert_hash.get(name,""))==value) for name,value in expected.items()}
    coverage_ok=bool(pre.get("domain_complete",False))
    numerical_ok=bool(pre.get("numerical_controlled",False))
    refinement_ok=bool(pre.get("refinement_pass",False))
    valid=bool(
        schema_ok and quantum_ok and no_fit and expected_present
        and all(hash_match.values()) and coverage_ok and numerical_ok and refinement_ok
    )
    return {
        "schema": schema,
        "schema_ok": schema_ok,
        "legacy_target_only_schema": False,
        "quantum_habitat_residual_certified": quantum_ok,
        "lambda_was_not_fitted": no_fit,
        "expected_hashes_present": expected_present,
        "hash_matches": hash_match,
        "coverage_complete": coverage_ok,
        "numerical_controlled": numerical_ok,
        "refinement_pass": refinement_ok,
        "valid_for_this_master": valid,
        "reason": "matched actual quantum HDA residual certificate" if valid else "quantum HDA certificate does not satisfy same-master fail-closed contract",
    }


def load_hda_certificate(base:Path, manifest:dict, master_hash):
    rel=manifest.get("quantum_hda_residual_certificate") or manifest.get("dtarget_hda_certificate")
    if not rel:
        return None,{"present":False,"valid_for_this_master":False,"reason":"no quantum HDA residual certificate supplied"}
    p=(base/str(rel)).resolve(); cert=json.loads(p.read_text(encoding="utf-8"))
    audit=audit_hda_certificate(cert,manifest,master_hash)
    audit.update({
        "present":True,
        "path":str(rel),
        "certificate_sha256":hashlib.sha256(p.read_bytes()).hexdigest(),
    })
    return cert,audit


def assemble_manifest(manifest_path: Path):
    m=json.loads(manifest_path.read_text(encoding="utf-8")); base=manifest_path.parent
    dim=int(m["domain_dimension"]); nodes=[int(x) for x in m["nodes"]]; lam=float(m["lambda_L"])
    domain_complete=bool(m.get("domain_complete",False)); domain_label=m.get("domain_label")
    E={v:[None]*dim for v in nodes}; L={v:[None]*dim for v in nodes}; hashes=[]
    for row in m["columns"]:
        fam=str(row["family"]); v=int(row["node"]); i=int(row["input_index"])
        if fam not in ("E","L") or v not in E or not (0<=i<dim): raise ValueError(f"bad E/L column descriptor {row}")
        p=(base/row["path"]).resolve(); st,_=load_column(p); hashes.append(hashlib.sha256(p.read_bytes()).hexdigest())
        target=E if fam=="E" else L
        if target[v][i] is not None: raise ValueError(f"duplicate {fam},node={v},input={i}")
        target[v][i]=st
    missing=[]
    for fam,target in (("E",E),("L",L)):
        for v in nodes:
            for i,x in enumerate(target[v]):
                if x is None: missing.append([fam,v,i])
    if missing: raise ValueError(f"missing E/L columns ({len(missing)}): {missing[:12]}")

    MEE,MEL,MLL,MH,rows=assemble(E,L,lam)

    D={}
    for row in m.get("dtarget_columns",[]):
        label=str(row["constraint"]); i=int(row["input_index"])
        if not (0<=i<dim): raise ValueError(f"bad D_target input index {row}")
        if label not in D: D[label]=[None]*dim
        if D[label][i] is not None: raise ValueError(f"duplicate D_target,{label},input={i}")
        p=(base/row["path"]).resolve(); st,_=load_column(p); hashes.append(hashlib.sha256(p.read_bytes()).hexdigest())
        D[label][i]=st
    MD,drows=assemble_dtarget(D,dim) if D else (np.zeros((dim,dim),complex),[])
    explicit_dtarget=bool(D)
    Mfull=0.5*(MH+MD+(MH+MD).conj().T)
    pencil_hash=hash_arrays(MEE,MEL,MLL,MD)
    cert,cert_audit=load_hda_certificate(base,m,pencil_hash)
    hda_closed=bool(cert_audit["valid_for_this_master"])
    projector_allowed=bool(domain_complete and hda_closed)
    a=spectral_audit(Mfull)
    mixed_ratio=float(np.linalg.norm(MEL)/max(np.linalg.norm(MEE)+lam*lam*np.linalg.norm(MLL),1e-300))

    if not domain_complete:
        status="RESTRICTED_DOMAIN_MASTER_DIAGNOSTIC"
    elif not hda_closed and explicit_dtarget:
        status="COMPLETE_FINITE_FULL_DIRAC_MASTER_HDA_UNCERTIFIED"
    elif not hda_closed:
        status="COMPLETE_DOMAIN_MASTER_HDA_UNCERTIFIED"
    elif explicit_dtarget:
        status="COMPLETE_FINITE_FULL_DIRAC_MASTER_QUANTUM_HDA_CERTIFIED"
    else:
        status="COMPLETE_FINITE_NORMAL_MASTER_QUANTUM_HDA_CERTIFIED"

    result={
        "status":status,"passed":True,"domain_label":domain_label,"domain_dimension":dim,"domain_complete":domain_complete,
        "habitat_hash":m.get("habitat_hash"),"domain_hash":m.get("domain_hash"),"convention_hash":m.get("convention_hash"),
        "nodes":nodes,"lambda_L":lam,"master_pencil_hash":pencil_hash,
        "column_file_hashes":hashes,"per_node":rows,"dtarget_per_constraint":drows,
        "explicit_dtarget_columns_in_master":explicit_dtarget,"quantum_hda_closed":hda_closed,
        "quantum_hda_certificate_audit":cert_audit,"mixed_block_relative_norm":mixed_ratio,
        "normal_master_frobenius_norm":float(np.linalg.norm(MH)),"dtarget_master_frobenius_norm":float(np.linalg.norm(MD)),
        "spectrum":{"rank":a["rank"],"nullity":a["nullity"],"rank_tolerance":a["rank_tolerance"],
                    "eigenvalue_min":float(np.min(a["eigenvalues"])),"eigenvalue_max":float(np.max(a["eigenvalues"])),
                    "smallest_positive":a["smallest_positive"],"condition_number_on_support":a["condition_number_on_support"],
                    "eigenvalues":[float(x) for x in a["eigenvalues"]]},
        "physical_projector_emitted":projector_allowed,
        "claim_boundary":"P_phys is emitted only for a complete finite habitat with a matching actual quantum HH-Dtarget residual certificate. Explicit Dtarget columns alone do not certify HDA closure."
    }
    if projector_allowed:
        P=a["P0"]
        result["physical_projector"]=[[[float(z.real),float(z.imag)] for z in row] for row in P]
    return result


def matrix_columns(A):
    out=[]
    for j in range(A.shape[1]): out.append({(i,):complex(A[i,j]) for i in range(A.shape[0]) if abs(A[i,j])>1e-14})
    return out


def selftest():
    rng=np.random.default_rng(5092026); dim=6; E={};L={};dE={};dL={}
    for v in range(3):
        A=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5)); A=.5*(A+A.conj().T)+(2+.25*v)*np.eye(5)
        B=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5)); B=.5*(B+B.conj().T)+(1.3+.17*v)*np.eye(5)
        Ae=np.zeros((dim,dim),complex); Bl=np.zeros_like(Ae); Ae[1:,1:]=A; Bl[1:,1:]=B
        dE[v]=Ae;dL[v]=Bl;E[v]=matrix_columns(Ae);L[v]=matrix_columns(Bl)
    lam=.73; MEE,MEL,MLL,MH,_=assemble(E,L,lam)
    direct=sum((dE[v]+lam*dL[v]).conj().T@(dE[v]+lam*dL[v]) for v in range(3)); a=spectral_audit(MH)
    e0=np.zeros(dim);e0[0]=1
    algebra_ok=np.linalg.norm(MH-direct)<2e-10 and a["nullity"]==1 and np.linalg.norm(a["P0"]-np.outer(e0,e0))<2e-10

    # Restricted-domain negative control: compressed nullity can miss a full zero sector.
    Mfull=np.diag([1.,0.]); B=np.array([[1.],[0.]]); ar=spectral_audit(B.T@Mfull@B); af=spectral_audit(Mfull)
    restricted_ok=ar["nullity"]==0 and af["nullity"]==1

    # Explicit D_target is a constraint contribution, not an HDA proof.
    D0=np.diag([0.,1.,1.,1.,1.,1.]); MD,_=assemble_dtarget({"D_target_control":matrix_columns(D0)},dim)
    afd=spectral_audit(MH+MD)
    explicit_d_spectrum_ok=afd["nullity"]==1 and np.linalg.norm(afd["P0"]-np.outer(e0,e0))<2e-10

    master_hash=hash_arrays(MEE,MEL,MLL,MD)
    manifest={
        "habitat_hash":"hab:test","domain_hash":"dom:test","convention_hash":"conv:test"
    }
    legacy={
        "schema":LEGACY_DTARGET_SCHEMA,
        "quantum_habitat_residual_certified":True,
        "certified_for_physical_projector":True,
    }
    legacy_audit=audit_hda_certificate(legacy,manifest,master_hash)
    legacy_rejected=not legacy_audit["valid_for_this_master"]

    new_cert={
        "schema":QUANTUM_HDA_SCHEMA,
        "quantum_habitat_residual_certified":True,
        "hashes":{
            "habitat_hash":"hab:test",
            "domain_hash":"dom:test",
            "constraint_packet_hash":master_hash,
            "convention_hash":"conv:test",
        },
        "preconditions":{
            "domain_complete":True,
            "numerical_controlled":True,
            "refinement_pass":True,
            "lambda_was_not_fitted":True,
        },
        "norm_squared_polynomial":{"lambda_fit_used_for_certification":False},
    }
    new_audit=audit_hda_certificate(new_cert,manifest,master_hash)
    new_authorized=new_audit["valid_for_this_master"]
    wrong_hash=json.loads(json.dumps(new_cert)); wrong_hash["hashes"]["domain_hash"]="dom:wrong"
    wrong_hash_rejected=not audit_hda_certificate(wrong_hash,manifest,master_hash)["valid_for_this_master"]

    fail_closed_ok=explicit_d_spectrum_ok and legacy_rejected and new_authorized and wrong_hash_rejected

    return {"status":"BQG E/L/D_target outgoing-column master assembler regression",
            "passed":bool(algebra_ok and restricted_ok and fail_closed_ok),
            "formula":"M_full=M_EE+lambda M_EL+lambda^2 M_LL+sum_I D_I^dagger D_I",
            "complete_normal_habitat_algebra_control":{"passed":bool(algebra_ok),"normal_nullity":a["nullity"],"normal_master_hash":hash_arrays(MEE,MEL,MLL)},
            "restricted_domain_negative_control":{"passed":bool(restricted_ok),"restricted_nullity":ar["nullity"],"full_nullity":af["nullity"]},
            "hda_fail_closed_control":{
                "passed":bool(fail_closed_ok),
                "explicit_Dtarget_spectrum_is_well_defined":bool(explicit_d_spectrum_ok),
                "explicit_Dtarget_alone_authorized":False,
                "legacy_target_only_certificate_rejected":bool(legacy_rejected),
                "new_quantum_residual_certificate_authorized":bool(new_authorized),
                "mismatched_hash_certificate_rejected":bool(wrong_hash_rejected),
            },
            "production_rule":"P_phys requires domain_complete=true AND a same-habitat/same-domain/same-master/same-convention BQG_QUANTUM_HDA_RESIDUAL_CERTIFICATE_V1. Explicit Dtarget columns alone never certify [H,H]=Dtarget."}


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("--manifest",type=Path); ap.add_argument("--output",type=Path,default=Path("verification_results/BQG_CONSTRAINT_MASTER_ASSEMBLER.json")); a=ap.parse_args()
    out=assemble_manifest(a.manifest) if a.manifest else selftest(); a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8"); print(json.dumps({k:v for k,v in out.items() if k!="physical_projector"},indent=2)); return 0 if out["passed"] else 1

if __name__=="__main__": raise SystemExit(main())
