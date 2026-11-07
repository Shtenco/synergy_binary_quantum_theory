#!/usr/bin/env python3
"""V21 fail-closed C1 acceptance checks. No independent SU(2) amplitudes."""
from __future__ import annotations
import math

BOUND=1e-8
PHYSICAL_LEAK_KEYS=(
    'CV_complete_basis_leakage',
    'outer_complete_basis_leakage',
    'internal_volume_sector_leakage',
    'CK_outer_complete_basis_leakage',
    'CK_internal_volume_sector_leakage'
)
CHARGE_KEYS=('complete_charge_basis_leakage','CK_complete_charge_basis_leakage')

def validate_physical_stage(diag, state, stage, threshold=BOUND):
    if not isinstance(diag,dict) or not isinstance(state,dict):
        raise ValueError('BAD_STAGE_OBJECT')
    if stage not in (1,2,3):raise ValueError('BAD_STAGE')
    for key,value in diag.items():
        if not isinstance(value,(int,float)) or not math.isfinite(value):
            raise ValueError('NONFINITE_OR_BAD_LEAKAGE')
        if value<0:raise ValueError('NEGATIVE_LEAKAGE')
    relevant=[float(diag[k]) for k in PHYSICAL_LEAK_KEYS if k in diag]
    if not relevant:raise ValueError('MISSING_PHYSICAL_LEAKAGE_DIAGNOSTIC')
    maximum=max(relevant)
    if maximum>=threshold:
        raise ValueError('PHYSICAL_BASIS_LEAKAGE_EXCEEDS_PRODUCTION_THRESHOLD')
    if not all(isinstance(z,(complex,int,float)) for z in state.values()):
        raise ValueError('BAD_AMPLITUDE_TYPE')
    if not all(math.isfinite(z.real) and math.isfinite(z.imag) for z in state.values()):
        raise ValueError('NONFINITE_AMPLITUDE')
    charge=max((float(diag[k]) for k in CHARGE_KEYS if k in diag),default=None)
    return {
        'stage_physical_basis_pass':True,
        'physical_leakage_max':maximum,
        'charge_leakage_max_unprojected':charge,
        'charge_covariance_verified':False,
        'stage_zero_support':len(state)==0,
        'stage_nonempty_not_full_pair':True,
        'threshold':threshold,
        'scope':'intermediate precursor only; no full physical Lorentzian state'
    }
