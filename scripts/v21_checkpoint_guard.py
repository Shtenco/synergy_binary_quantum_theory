#!/usr/bin/env python3
"""Fail-closed V21 stage acceptance. This validates diagnostics, not physics closure."""
from __future__ import annotations
import math

BOUND = 1e-8
CV_KEYS = ("CV_complete_basis_leakage",)
CK_KEYS = ("outer_complete_basis_leakage",
           "internal_volume_sector_leakage",
           "complete_charge_basis_leakage")


def _finite_nonnegative(value):
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(float(value)) and float(value) >= 0.0)


def validate_physical_stage(diag, state, stage, operator_signature,
                            threshold=BOUND):
    if not isinstance(diag, dict) or not isinstance(state, dict):
        raise ValueError("BAD_STAGE_OBJECT")
    if (not isinstance(stage, int) or isinstance(stage, bool)
            or stage not in (1, 2, 3)):
        raise ValueError("BAD_STAGE")
    if (not isinstance(threshold, (int, float)) or isinstance(threshold, bool)
            or not math.isfinite(float(threshold)) or threshold <= 0):
        raise ValueError("BAD_THRESHOLD")
    if (not isinstance(operator_signature, (list, tuple))
            or len(operator_signature) != 4
            or operator_signature[0] not in ("CV", "CK")
            or any(not isinstance(value, int) or isinstance(value, bool)
                   for value in operator_signature[1:])):
        raise ValueError("BAD_OPERATOR_SIGNATURE")

    for key, value in diag.items():
        if not isinstance(key, str) or not _finite_nonnegative(value):
            raise ValueError("NONFINITE_OR_BAD_LEAKAGE")

    required = CV_KEYS if operator_signature[0] == "CV" else CK_KEYS
    missing = [key for key in required if key not in diag]
    if missing:
        raise ValueError("MISSING_OPERATION_SPECIFIC_DIAGNOSTIC:" + ",".join(missing))
    physical_keys = CV_KEYS if operator_signature[0] == "CV" else CK_KEYS[:2]
    maximum = max(float(diag[key]) for key in physical_keys)
    if maximum >= float(threshold):
        raise ValueError("PHYSICAL_BASIS_LEAKAGE_EXCEEDS_PRODUCTION_THRESHOLD")

    if not state:
        raise ValueError("EXACT_ZERO_STAGE_OUTPUT")
    for amplitude in state.values():
        if not isinstance(amplitude, (complex, int, float)):
            raise ValueError("BAD_AMPLITUDE_TYPE")
        if not (math.isfinite(complex(amplitude).real)
                and math.isfinite(complex(amplitude).imag)):
            raise ValueError("NONFINITE_AMPLITUDE")
    state_norm = 0.0
    for amplitude in state.values():
        state_norm = math.hypot(state_norm, abs(complex(amplitude)))
    if not math.isfinite(state_norm):
        raise ValueError("NONFINITE_STATE_NORM")
    if state_norm == 0.0:
        raise ValueError("EXACT_ZERO_STAGE_OUTPUT")

    charge = (float(diag["complete_charge_basis_leakage"])
              if operator_signature[0] == "CK" else None)
    return {
        "stage_physical_basis_pass": True,
        "operation": operator_signature[0],
        "operator_signature": list(operator_signature),
        "required_diagnostics": list(required),
        "physical_leakage_max": maximum,
        "charge_leakage_max_unprojected": charge,
        "charge_covariance_verified": False,
        "stage_zero_support": False,
        "stage_nonempty_not_full_pair": True,
        "threshold": float(threshold),
        "scope": "nonzero intermediate precursor only; no full physical Lorentzian state"
    }
