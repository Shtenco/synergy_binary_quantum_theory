#!/usr/bin/env python3
"""BQG V26: reproduce K5 H_E 12-orientation conditioning from real archived logs.

Only individual norm summaries are available, not the 12 full coherent state
vectors. Therefore this script CANNOT compute full H_E charge leakage.
"""
from __future__ import annotations
import json,math,sys
from pathlib import Path
from fractions import Fraction

ROOT=Path(__file__).resolve().parent
SOURCE_BLOB='078c311c74c0131d09ecacefe57b9eb7735b30ea'
SOURCE_RUN=42136642363

def load():
    s=json.loads((ROOT/'K5_CHARGE_12_SCALARS.json').read_text())
    assert s['source_git_blob_sha']==SOURCE_BLOB and s['upstream_run']==SOURCE_RUN
    rows=s['rows']
    if len(rows)!=12 or [row['index'] for row in rows]!=list(range(12)):
        raise ValueError('INCOMPLETE_ORIENTATION_CENSUS')
    if any(row['branch_count']!=6144 or row['raw2']<0 or row['projected2']<0 or row['incoherent2']<0 for row in rows):
        raise ValueError('INCOMPLETE_BRANCHES_OR_INVALID_NORMS')
    return rows

def diagnostics(rows):
    norm2=sum((Fraction(str(r['raw2'])) for r in rows),Fraction(0))
    incoh2=sum((Fraction(str(r['incoherent2'])) for r in rows),Fraction(0))
    summary=[]
    for r in rows:
        raw=r['raw2'];incoh=r['incoherent2'];projected=r['projected2']
        ratio=None if raw==0 else incoh/raw
        # Ratio of incoherent to coherent SQUARED norms: a cancellation/conditioning
        # diagnostic, NOT a forward-error bound for any particular BLAS call.
        bits=None if ratio is None else math.log2(ratio)
        summary.append({
            'index':r['index'],'raw2':raw,'projected2':projected,
            'incoherent2':incoh,'stability':r['stability'],
            'incoherent_to_coherent_squared_ratio':ratio,
            'scale_separation_bits':bits,
            'impossible_exact_projector_norm_if_values_exact':projected>raw,
            'possible_not_certified_zero':raw==0,
            'projector_subtraction_defect':raw-projected,
        })
    return {'source':'Actual 12 original-oriented SU2 H_E^sine scalar log; no complete cross-orientation vectors',
            'source_git_blob_sha':SOURCE_BLOB,'upstream_run':SOURCE_RUN,
            'sum_raw2_from_printed_decimal_exact':str(norm2),
            'sum_incoherent2_from_printed_decimal_exact':str(incoh2),
            'rows':summary,'noncertified': [r['index'] for r in summary if r['stability']=='CANCELLATION_DOMINATED_NOT_CERTIFIED'],
            'full_charge_projection_proved':False,
            'all_12_coherent_output_vectors_available':False}

def minimum_consistent_uniform_eta(rows):
    """Required uniform absolute error on each of the two squared norms.

    If a COMMON true orthogonal projector exists, true_raw2>=true_proj2.
    Repairing any printed projected2>raw2 discrepancy therefore needs at least
    half that gap in the two independently rounded norm-squares.
    """
    return max(0.,*(0.5*(r['projected2']-r['raw2']) for r in rows))

def bound_absolute_charge_leak(rows, absolute_error_each):
    """Conditional bound from a COMMON orthogonal charge projector.

    Assume for each i BOTH true squared norms are within eta of their printed
    values; assume one common P for all 12 orientations and same Hilbert norm.
    Then ||(I-P)sum_i psi_i|| <= sum_i sqrt(max(0,raw_i^2-proj_i^2)+2eta).
    eta is a USER-PROVIDED ASSUMPTION, NOT measured SU(2) accuracy.
    """
    eta=float(absolute_error_each)
    if not math.isfinite(eta) or eta <0:
        raise ValueError('INVALID_SQUARED_NORM_ABSOLUTE_ERROR_ASSUMPTION')
    if eta < minimum_consistent_uniform_eta(rows):
        raise ValueError('NO_COMMON_ORTHOGONAL_PROJECTOR_WITH_ASSUMED_ACCURACY')
    return math.fsum(math.sqrt(max(0.,r['raw2']-r['projected2'])+2*eta) for r in rows)

def required_uniform_eta_for_absolute_charge_leak(target,rows):
    """Bisection of sufficient uniform squared-norm absolute error budget.

    If observed max(0,raw2-proj2) terms already exceed target, no eta>=0
    can make this *particular* conservative triangle inequality certify target.
    """
    if not(math.isfinite(target) and target>0):raise ValueError('BAD_TARGET')
    lo=minimum_consistent_uniform_eta(rows)
    if bound_absolute_charge_leak(rows,lo)>target:return None
    hi=1.
    for _ in range(250):
        mid=(lo+hi)/2
        if bound_absolute_charge_leak(rows,mid)>target:hi=mid
        else:lo=mid
    return lo

def main():
    rows=load(); d=diagnostics(rows)
    eta_floor=minimum_consistent_uniform_eta(rows)
    d['minimum_uniform_squared_norm_error_to_reconcile_projector_monotonicity']=eta_floor
    d['best_attainable_absolute_leakage_upper_from_this_triangle_method']=bound_absolute_charge_leak(rows,eta_floor)
    d['conditional_absolute_charge_leak_bounds']={f'eta_{e:g}':bound_absolute_charge_leak(rows,e)
        for e in [eta_floor,1e-12,1e-14,1e-15]}
    d['required_uniform_abs_error_in_each_squared_norm_for_absolute_target']={
        f'target_{e:g}':required_uniform_eta_for_absolute_charge_leak(e,rows)
        for e in [1e-4,1e-6,1e-7,1e-8,1e-9]}
    d['true_global_charged_leak_unknown']=True
    (ROOT/'V26_CHARGE_DIAGNOSTICS.json').write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in d.items() if k!='rows'},indent=2,ensure_ascii=False))
if __name__=='__main__':main()
