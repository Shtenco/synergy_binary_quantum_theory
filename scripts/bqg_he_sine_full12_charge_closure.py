#!/usr/bin/env python3
"""True full signed H_E^sine on one real charged K5 spin-network state.

Runs ALL 12 oriented sine contributions (64 signed primitive terms each) using
the original CK/SINE kernel. No projected-free equivalence is assumed:
this calculates the full projected wavefunction. The separate raw-vs-projected
norm audits of all 12 orientations are loaded as independent data.

The orthogonal-projector triangle inequality supplies a *conditional*
full-sum leakage envelope, but a floating-point norm-difference near zero
is not a validated interval bound. Thus no false physical proof.
"""
import argparse,json,math,sys,collections,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import numpy as np
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_covariant_K_sine_composition_gate as SINE
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

AUDIT=Path("research/node1_fullstate_gate_2026_10_31/HE_SINE_12_ORIENTATIONS_TRUE_POSTSUM_CHARGE_20261106.json")
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=int,default=1)
    ap.add_argument("--target",type=int,default=0)
    ap.add_argument("--Jmax2",type=int,default=3)
    opt=ap.parse_args()
    ZVM.patch_and_clear()
    init=PW.basis_full_jhalf()[0]
    input_state,inv_leak,wrong,_=CK.inverse_outer_complete(init,opt.source,opt.target,0,0,opt.Jmax2)
    assert input_state
    key=max(input_state,key=lambda k:(abs(input_state[k]),str(k)))
    coeff=input_state[key]
    t=time.perf_counter()
    projected,volume_leak,primitive_charge_leak=SINE.complete_HE_sine(key,opt.source,opt.Jmax2,(opt.source,opt.target))
    projected={k:coeff*a for k,a in projected.items()}
    norm2=float(sum(abs(a)**2 for a in projected.values()))
    by_spins=collections.Counter()
    for (sp,labels),amp in projected.items():
        by_spins[tuple(sp)] += abs(amp)**2
        assert labels[opt.source][0]=="C" and labels[opt.target][0]=="C"
        assert all(labels[v][0]=="G" for v in PW.VERT if v not in (opt.source,opt.target))
    # HE gauge-covariance implies J expected unchanged at both charged endpoints,
    # but we MEASURE it on the projected answer; it is separate from raw projection leakage.
    wrong_source,weights_s=CK.wrong_charge_fraction(projected,opt.source,expected_J2=1)
    wrong_target,weights_t=CK.wrong_charge_fraction(projected,opt.target,expected_J2=1)
    a=json.loads(AUDIT.read_text(encoding="utf-8"))
    assert len(a["rows"])==12
    near_zero=[x["index"] for x in a["rows"] if x["numerical_stability"]!="MEASURABLE_COHERENT_OUTPUT"]
    robust=[x for x in a["rows"] if x["numerical_stability"]=="MEASURABLE_COHERENT_OUTPUT"]
    if len(robust)+len(near_zero)!=12:raise RuntimeError("not all genuine orientations audited")
    # This is a diagnostic uncertainty scale, not an interval-arithmetic proof.
    # Each coherent norm squared is computed through factorized exact SU2 inner
    # products, with floating point cancellation floor from incoherent sum.
    eps_sq=1e-12
    err_est=[]
    for x in a["rows"]:
        raw=x["raw2"]
        proj=x["projected2"]
        diff=max(0.,raw-proj)
        floor=eps_sq*(raw+x["incoherent2"])
        err_est.append(math.sqrt(diff+floor))
    leakage_absolute_diagnostic=sum(err_est)
    leak_ratio_diagnostic=leakage_absolute_diagnostic/max(math.sqrt(norm2),1e-30)
    res=dict(
      status="ACTUAL_FULL_12_ORIENTATION_PROJECTED_SINE_HE_FROM_ONE_GENUINE_CHARGED_K5_SOURCE",
      Jmax2=opt.Jmax2,source=opt.source,target=opt.target,
      original_input_basis=repr(key),
      original_input_complex_coefficient=[float(coeff.real),float(coeff.imag)],
      input_support_from_inverse=len(input_state),inverse_basis_leak=inv_leak,
      complete_signed_primitive_sequence_count=12*2*len(PW.T_sequences(*PW.oriented_specs(opt.source)[0][1])),
      H_E_sine_projected_support=len(projected),
      H_E_sine_projected_coherent_norm=math.sqrt(norm2),
      H_E_sine_projected_coherent_norm_squared=norm2,
      H_E_sine_output_spin_sectors=len(by_spins),
      H_E_sine_wrong_endpoint_charge_source_fraction=wrong_source,
      H_E_sine_wrong_endpoint_charge_target_fraction=wrong_target,
      H_E_sine_source_charge_weights=weights_s,
      H_E_sine_target_charge_weights=weights_t,
      measured_max_primitive_basis_leak=primitive_charge_leak,
      measured_max_primitive_volume_leak=volume_leak,
      independently_verified_individual_orientations=len(a["rows"]),
      numerically_cancellation_dominated_orientations=near_zero,
      one_orientation_postsum_residual_numeric_envelope=leakage_absolute_diagnostic,
      full_relative_charge_residual_estimate_NOT_CERTIFIED=leak_ratio_diagnostic,
      numeric_envelope_floor_factor=eps_sq,
      proof_scope="The raw-vs-projected identity of individual orientations implies the full linear sum is in the common range of P only in exact arithmetic. Here checked numerically, not interval-certified. No physical HDA, TT or group averaging.",
      elapsed_s=time.perf_counter()-t,
    )
    print("BQG_FULL_HE_SINE_CHARGE_BEGIN",flush=True)
    print(json.dumps(res,indent=2,ensure_ascii=False),flush=True)
    print("BQG_FULL_HE_SINE_CHARGE_END",flush=True)
if __name__=="__main__":main()
