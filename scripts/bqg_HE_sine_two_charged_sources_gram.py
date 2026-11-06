#!/usr/bin/env python3
"""Exact SHA verified 2x2 Gram from original full HE sine complex SU2 outputs.
Not physical Lorentzian constraint or HDA: positive Euclidean Gram alone proves
no common physical kernel on the whole untruncated Hilbert space.
"""
import base64,gzip,hashlib,json,math
from pathlib import Path
import numpy as np
R=Path("research/node1_fullstate_gate_2026_10_31/he_sine_spin_refinement")
FILES=[R/"sine_source1_target0_jmax2_3.json",R/"sine_source1_target0_jmax2_3_basis1.json"]
def load(p):
    packet=json.loads(p.read_text(encoding="utf-8"))
    meta=packet["metadata"]
    raw=gzip.decompress(base64.b64decode(packet["gzip_base64"],validate=True))
    assert hashlib.sha256(raw).hexdigest()==meta["complex_output_SHA256"]
    data=json.loads(raw)
    assert data["source"]==1 and data["target"]==0 and data["Jmax2"]==3
    v={}
    for spins,labels,re,im in data["states"]:
        key=json.dumps([spins,labels],separators=(",",":"))
        assert key not in v
        v[key]=complex(re,im)
    c=complex(*data["basis_amp"])
    assert c!=0
    assert abs(np.sqrt(sum(abs(x)**2 for x in v.values()))-meta["output_norm"])<1e-12
    return v,data,c
def main():
    inputs=[load(p) for p in FILES]
    keyA=inputs[0][1]["basis_key"];keyB=inputs[1][1]["basis_key"]
    assert keyA!=keyB and keyA[0]!=keyB[0],"Two input basis states not independent"
    # Unweight by genuine complex input coefficients to obtain an ON source basis.
    columns=[{k:z/c for k,z in out.items()} for out,d,c in inputs]
    overlap=np.zeros((2,2),complex)
    for i,x in enumerate(columns):
        for j,y in enumerate(columns):
            overlap[i,j]=sum(np.conj(x[k])*y[k] for k in x.keys()&y.keys())
    assert np.linalg.norm(overlap-overlap.conj().T)<1e-12
    eig=np.linalg.eigvalsh(overlap)
    assert eig[0]>=-1e-12
    cov={k:inputs[0][0].get(k,0j)+inputs[1][0].get(k,0j) for k in set(inputs[0][0])|set(inputs[1][0])}
    coherent_norm=math.sqrt(sum(abs(x)**2 for x in cov.values()))
    a=inputs[0][0];b=inputs[1][0]
    expected=sum(abs(x)**2 for x in a.values())+sum(abs(x)**2 for x in b.values())
    res=dict(status="ACTUAL_EUCLIDEAN_HE_SINE_TWO_CHARGED_SOURCE_GRAM",
       input_sources=2,input_complex_coefficients=[[c.real,c.imag] for _,_,c in inputs],
       weighted_output_norms=[math.sqrt(sum(abs(x)**2 for x in v.values())) for v,_,_ in inputs],
       coherent_weighted_sum_norm=coherent_norm,
       signed_interference_squared_norm=coherent_norm**2-expected,
       union_outgoing_SU2_keys=len(cov),
       unweighted_two_by_two_Gram=[[[float(x.real),float(x.imag)] for x in row] for row in overlap],
       Gram_eigenvalues=eig.tolist(),
       smallest_singular_value=float(np.sqrt(max(0,eig[0]))),
       largest_singular_value=float(np.sqrt(max(0,eig[-1]))),
       finite_subspace_HE_kernel_empty=bool(eig[0]>1e-12),
       proof_scope="two actual charged inverse-holonomy basis vectors, full 12 H_E sine, Gauss/projected; NOT full Lorentzian physical master, not independent D/HDA/TT.")
    print("BQG_EUCLIDEAN_HE_TWO_SOURCE_GRAM_BEGIN")
    print(json.dumps(res,indent=2))
    print("BQG_EUCLIDEAN_HE_TWO_SOURCE_GRAM_END")
if __name__=="__main__":main()
