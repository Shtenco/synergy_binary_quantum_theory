#!/usr/bin/env python3
"""Independent exact complex-wavefunction convergence at Jmax2=3,5,7.
Reads immutable source-file evidence from original GitHub Actions, verifies
SHA256 of DECOMPRESSED complete output and all complex SU2 basis labels.
"""
import base64,gzip,hashlib,json,math
from pathlib import Path
BASE=Path("research/node1_fullstate_gate_2026_10_31/he_sine_spin_refinement")
def read(j):
    p=BASE/f"sine_source1_target0_jmax2_{j}.json"
    outer=json.loads(p.read_text(encoding="utf-8"))
    raw=gzip.decompress(base64.b64decode(outer["gzip_base64"],validate=True))
    m=outer["metadata"]
    assert hashlib.sha256(raw).hexdigest()==m["complex_output_SHA256"]
    z=json.loads(raw)
    assert z["source"]==1 and z["target"]==0 and z["Jmax2"]==j
    assert len(z["states"])==m["output_support"]
    vals={}
    for spins,labels,re,im in z["states"]:
        key=json.dumps([spins,labels],separators=(",",":"))
        assert key not in vals
        vals[key]=complex(re,im)
    assert abs(math.sqrt(sum(abs(a)**2 for a in vals.values()))-m["output_norm"])<1e-12
    return vals,z,m

def norm(a):
    return math.sqrt(sum(abs(x)**2 for x in a.values()))
def difference(a,b):
    return {key:a.get(key,0j)-b.get(key,0j) for key in set(a)|set(b)}
def main():
    states={j:read(j) for j in (3,5,7)}
    assert states[3][1]["basis_key"]==states[5][1]["basis_key"]==states[7][1]["basis_key"]
    result={"status":"REAL_FULL_12_ORIENTATION_SINE_HE_REFINEMENT_ACTUAL_COMPLEX_OUTPUT"}
    for lo,hi in ((3,5),(5,7),(3,7)):
        a=states[lo][0];b=states[hi][0];delta=difference(a,b)
        result[f"delta_Jmax2_{lo}_{hi}"]={
          "norm":norm(delta),"relative_to_higher_norm":norm(delta)/max(norm(b),1e-30),
          "different_full_SU2_output_keys":sum(abs(v)>1e-12 for v in delta.values()),
          "support_lower":len(a),"support_higher":len(b),
          "common_keys":len(set(a)&set(b)),
          "lower_norm":norm(a),"higher_norm":norm(b)}
    result["scope"]="ONE charged K5 physical-source basis, same full 12 oriented HE sine under different spin truncations; NO complete full HL/HDA."
    print("BQG_HE_REFINEMENT_COMPARISON_BEGIN")
    print(json.dumps(result,indent=2))
    print("BQG_HE_REFINEMENT_COMPARISON_END")
    assert result["delta_Jmax2_3_5"]["relative_to_higher_norm"]<1e-9
    assert result["delta_Jmax2_5_7"]["relative_to_higher_norm"]<1e-9
if __name__=="__main__":main()
