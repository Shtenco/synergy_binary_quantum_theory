#!/usr/bin/env python3
"""Genuine full 12-orientation H_E^sine K5 spin-cutoff refinement source.

Preserve EXACT projected complex spin-network amplitude keys and SHA256 as
compressed JSON in GitHub Actions log. This is not yet the unprojected
post-sum leakage; Jmax refinement checks SAME physical source.
"""
import argparse,base64,gzip,hashlib,json,math,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_covariant_K_sine_composition_gate as SINE
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

def run(jmax,source,target):
    ZVM.patch_and_clear()
    initial=PW.basis_full_jhalf()[0]
    state,leak,wrong,_=CK.inverse_outer_complete(initial,source,target,0,0,jmax)
    assert state and leak<1e-10
    key=max(state,key=lambda k:(abs(state[k]),str(k)))
    amp=state[key]
    t=time.perf_counter()
    projected,vol,primitive=SINE.complete_HE_sine(key,source,jmax,(source,target))
    full={k:amp*v for k,v in projected.items()}
    wrong_s,ws=CK.wrong_charge_fraction(full,source)
    wrong_t,wt=CK.wrong_charge_fraction(full,target)
    assert wrong_s<1e-10 and wrong_t<1e-10
    records=[]
    for (spins,labels),val in sorted(full.items(),key=lambda kv:str(kv[0])):
        records.append([list(map(int,spins)),[list(label) for label in labels],float(val.real),float(val.imag)])
    original=json.dumps({"basis_key":[list(key[0]),[list(q) for q in key[1]]],
                         "basis_amp":[float(amp.real),float(amp.imag)],
                         "source":source,"target":target,"Jmax2":jmax,
                         "states":records},separators=(",",":"),ensure_ascii=True).encode()
    compact=gzip.compress(original,compresslevel=9,mtime=0)
    sha=hashlib.sha256(original).hexdigest()
    info={
      "status":"GENUINE_FULL_12_ORIENTATION_HE_SINE_REFINEMENT",
      "source":source,"target":target,"Jmax2":jmax,
      "input_basis_key":str(key),"input_complex_amp":[float(amp.real),float(amp.imag)],
      "output_support":len(full),"output_norm":math.sqrt(sum(abs(z)**2 for z in full.values())),
      "primitive_charge_leak_max":primitive,
      "primitive_volume_leak_max":vol,
      "projected_wrong_endpoint_charge_source":wrong_s,
      "projected_wrong_endpoint_charge_target":wrong_t,
      "complex_output_SHA256":sha,
      "original_encoded_bytes":len(original),"compressed_bytes":len(compact),
      "elapsed_seconds":time.perf_counter()-t,
      "scope":"full 12 oriented H_E sine on ONE input, charged projector applied; no raw unprojected charge closure for Jmax>3, no physical HDA/TT"
    }
    print("BQG_HE_REFINEMENT_BEGIN",flush=True)
    print("BQG_HE_REFINEMENT_METADATA",json.dumps(info,separators=(",",":")),flush=True)
    print("BQG_HE_REFINEMENT_GZIP_BASE64",base64.b64encode(compact).decode(),flush=True)
    print("BQG_HE_REFINEMENT_END",flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--Jmax2",type=int,required=True)
    p.add_argument("--source",type=int,default=1)
    p.add_argument("--target",type=int,default=0)
    a=p.parse_args()
    run(a.Jmax2,a.source,a.target)
