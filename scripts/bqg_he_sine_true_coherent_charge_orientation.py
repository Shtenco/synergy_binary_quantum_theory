#!/usr/bin/env python3
"""Raw-versus-projected SU2 charged closure for an ACTUAL complete oriented H_E^sine term.

Computes both adjoint halves and ALL 32 T-sequences for ONE of 12 genuine
oriented source-node 1 HH terms. Unlike primitive-max leakage, it coherently
sums BEFORE measuring norm. Raw factorized tensors are contracted by
per-spin Gram blocks; the same ORIGINAL projection code determines P R.
Fails closed on truncation, unsupported source, negative norm, or basis mismatch.

IMPORTANT: One full oriented H_E term, not all 12, not 24-term H_L or physical HDA.
"""
from __future__ import annotations
import argparse,collections,json,math,time,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

def raw_gram_norm2(branches,block=64):
    # A product tensor has weight c and orthogonal spin-sector key.
    grouped=collections.defaultdict(list)
    for spins,tensors,c in branches:
        if abs(c):grouped[spins].append((tensors,complex(c)))
    total=0j
    max_group=max(map(len,grouped.values()),default=0)
    for spins,group in grouped.items():
        n=len(group)
        c=np.asarray([z[1] for z in group],complex)
        # cache exact untruncated raw local SU(2) amplitudes, common spin sector
        arrays=[np.asarray([np.ravel(x[0][v]) for x in group],complex) for v in PW.VERT]
        for i in range(0,n,block):
            I=slice(i,min(i+block,n))
            for j in range(0,n,block):
                J=slice(j,min(j+block,n))
                overlap=np.ones((I.stop-I.start,J.stop-J.start),complex)
                for X in arrays:
                    overlap*=np.conjugate(X[I]) @ X[J].T
                total+=np.einsum("i,ij,j->",np.conjugate(c[I]),overlap,c[J])
    if abs(total.imag)>2e-7*max(1.,abs(total.real)):
        raise AssertionError(f"raw norm imaginary: {total}")
    if total.real< -1e-9*max(1.,abs(total.real)):
        raise AssertionError(f"raw norm negative: {total}")
    return max(0.,float(total.real)),len(grouped),max_group

def choose_actual_source(source,target,Jmax2):
    initial=PW.basis_full_jhalf()[0]
    state,leak,wrong,weights=CK.inverse_outer_complete(initial,source,target,0,0,Jmax2)
    if not state:raise RuntimeError("Genuine inverse holonomy produced no charged input states")
    key=max(state,key=lambda k:(abs(state[k]),str(k)))
    assert len(key)==2 and key[1][source][0]=="C" and key[1][target][0]=="C"
    return key,state[key],len(state),float(leak),float(wrong)

def run(source,target,orient,Jmax2):
    ZVM.patch_and_clear()
    assert source!=target and target in PW.NEIG[source]
    specs=PW.oriented_specs(source)
    assert len(specs)==12 and 0<=orient<len(specs)
    key,amp,inv_size,invleak,wrong=choose_actual_source(source,target,Jmax2)
    sign,(v,a,b,c)=specs[orient]
    original=CK.branch_from_key(key,amp)
    charged=(source,target)
    projected={}
    raw=[]
    primitive_max=0.
    primitive_nonzero=0
    seq_total=0;raw_total=0;primitive_sum_squared=0.
    timer=time.perf_counter()
    # This is exactly the original sine prefactor for each of the 32 T terms,
    # with both adjoints, for one real oriented HH contribution.
    for adj in (False,True):
        pref=(.5j if adj else -.5j)*sign
        for coef,seq0 in PW.T_sequences(v,a,b,c):
            seq=PW.adjoint_sequence(seq0) if adj else seq0
            incoming,vleak=CK.apply_sequence_to_branch(original,seq,source,Jmax2)
            seq_total+=1
            for br in incoming:
                raw.append((br[0],br[1],br[2]*pref*coef))
                primitive_sum_squared += abs(br[2]*pref*coef)**2 * math.prod(float(np.vdot(T,T).real) for T in br[1])
                p,bleak=CK.project_branch_complete_charges(br,charged)
                CK.add_state(projected,p,scale=pref*coef,tol=0.0)
                primitive_max=max(primitive_max,bleak)
                if bleak>1e-8:primitive_nonzero+=1
                raw_total+=1
    if seq_total!=64:raise RuntimeError(f"not all 64 actual signed contributions visited: {seq_total}")
    raw2,groups,maxgroup=raw_gram_norm2(raw)
    projected2=sum(abs(x)**2 for x in projected.values())
    if projected2>raw2+max(1e-10,raw2*1e-8):
        raise AssertionError(f"True projection increased norm: raw={raw2}, projected={projected2}")
    leak2=max(0.,raw2-projected2)
    e=math.sqrt(leak2/max(raw2,1e-300)) if raw2>1e-25 and raw2>1e-12*primitive_sum_squared else None
    return {
      "status":"REAL_ONE_ORIENTATION_FULL_SINE_HE_PRESUM_CHARGE_TEST",
      "tested_operator":"ONE of 12 oriented actual H_E sine terms, full 32 T-sequences and both adjoints",
      "source":source,"target":target,"orientation_index":orient,
      "orientation_nodes":[v,a,b,c],"orientation_sign":sign,
      "Jmax2":Jmax2,"source_index":"max amplitude in genuine inverse holonomy state",
      "source_complex_amplitude":[float(amp.real),float(amp.imag)],
      "inverse_input_support":inv_size,"inverse_charge_projection_leak":invleak,
      "inverse_wrong_charge_fraction":wrong,
      "complete_signed_sequences":seq_total,"raw_factorized_branches":raw_total,
      "spin_sector_count":groups,"max_raw_branches_in_one_spin_sector":maxgroup,
      "maximum_primitive_charge_leak":primitive_max,
      "num_primitive_branches_charge_leak_gt_1e-8":primitive_nonzero,
      "true_unprojected_coherent_squared_norm":raw2,
      "true_projected_coherent_squared_norm":projected2,
      "postsum_relative_charge_leak":e,
      "sum_individual_raw_branch_squared_norm":primitive_sum_squared,
      "coherent_to_incoherent_squared_ratio":raw2/max(primitive_sum_squared,1e-300),
      "numerical_stability":"CANCELLATION_DOMINATED_NOT_CERTIFIED" if e is None else "MEASURABLE_COHERENT_OUTPUT",
      "elapsed_s":time.perf_counter()-timer,
      "caveat":"ONE orientation and ONE source basis state, NOT full H_E, not a physical HDA or TT test",
    }
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=int,default=1)
    p.add_argument("--target",type=int,default=0)
    p.add_argument("--orientation",type=int,default=0)
    p.add_argument("--Jmax2",type=int,default=3)
    a=p.parse_args()
    z=run(a.source,a.target,a.orientation,a.Jmax2)
    print("BQG_RAW_CHARGE_EXACT_BEGIN",flush=True)
    print(json.dumps(z,indent=2,ensure_ascii=False),flush=True)
    print("BQG_RAW_CHARGE_EXACT_END",flush=True)
