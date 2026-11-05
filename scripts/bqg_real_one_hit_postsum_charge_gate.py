#!/usr/bin/env python3
"""Real SU2 one-hit coherent charge projector audit on K5 spin networks.

Prerequisite, NOT the full post-summed HE_sine/H_L charge-closure certificate:
evaluate exactly ||sum_br raw_br||^2 by local SU(2) tensor-product overlaps,
and compare with the same sum projected into complete charged/uncharged
intertwiner basis. This is a small independent sanity test of the method
used later on the full coherent 24-term Lorentzian action.
"""
from __future__ import annotations
import argparse,itertools,json,math,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_covariant_K_leg_gate as CK

def rawnorm2(branches):
    # Hilbert product of factorized local Peter-Weyl tensors on common spin support.
    res=0j
    for a in branches:
        sa,ta,ca=a
        for b in branches:
            sb,tb,cb=b
            if sa!=sb:continue
            overlap=complex(np.conjugate(ca)*cb)
            for A,B in zip(ta,tb):overlap*=np.vdot(A,B)
            res+=overlap
    if abs(res.imag)>1e-8*max(1,abs(res.real)):
        raise ValueError("Unprojected norm not real: incompatible local basis")
    if res.real< -1e-9:raise ValueError("Negative norm")
    return max(0,float(res.real))

def projected(branches,charged):
    total={};primitive_leaks=[]
    for br in branches:
        p,leak=CK.project_branch_complete_charges(br,tuple(charged))
        primitive_leaks.append(leak)
        CK.add_state(total,p,tol=0.)
    return total,primitive_leaks

def test(source=1,target=0):
    assert target in PW.NEIG[source]
    init=PW.basis_full_jhalf()[0]
    original=PW.initial_factorized_oriented(init)
    charged=(source,target)
    pure,gauss_leak=projected([original],charged)
    if gauss_leak[0]>1e-8:
        raise AssertionError("Gauss input fails exact charged projector on source-target nodes")
    hit=list(PW.apply_hit_branch(original,target,source,0,0,7))
    if not hit:raise RuntimeError("No genuine K5 hit branches")
    raw2=rawnorm2(hit)
    projected_states,leaks=projected(hit,charged)
    proj2=sum(abs(z)**2 for z in projected_states.values())
    if proj2>raw2+1e-8*max(1,raw2):
        raise RuntimeError("Projection increased norm; inner product/conventions incompatible")
    ratio=math.sqrt(max(0.,raw2-proj2)/raw2) if raw2>1e-30 else None
    return dict(test="REAL_K5_ONE_HIT_POSTSUM_PROJECTOR_CONTROL",
      source_node=source,target_node=target,Jmax=3.5,
      gauss_pure_primitive_charge_leak=gauss_leak[0],
      number_real_primitive_hit_branches=len(hit),
      primitive_charge_leak_max=max(leaks),
      primitive_charge_leak_count_gt_1e_8=sum(x>1e-8 for x in leaks),
      raw_coherent_hit_norm2=raw2,projected_coherent_hit_norm2=proj2,
      coherent_one_hit_relative_charge_leak=ratio,
      distinct_projected_state_keys=len(projected_states),
      status="ONE_HIT_CONTROL_ONLY",
      caveat="This is one actual SU2 holonomy hit, NOT the complete full HE_sine branch sum, and NOT physical full HDA or TT.")
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=int,default=1)
    p.add_argument("--target",type=int,default=0)
    a=p.parse_args()
    r=test(a.source,a.target)
    print("BQG_REAL_ONE_HIT_CHARGE_BEGIN",flush=True)
    print(json.dumps(r,indent=2),flush=True)
    print("BQG_REAL_ONE_HIT_CHARGE_END",flush=True)
