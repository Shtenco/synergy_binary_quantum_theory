#!/usr/bin/env python3
"""Strict, dependency-light validators for V2 term acceptance artifacts."""
from __future__ import annotations
import math

REQUIRED_TERM_DIAGNOSTICS = (
    "CV_complete_basis_leakage",
    "CK_outer_complete_basis_leakage",
    "CK_internal_volume_sector_leakage",
    "CK_complete_charge_basis_leakage",
)


def validate_term_diagnostics(diag, accepted2, rejected2):
    """Require every operation-specific leakage diagnostic and valid norms."""
    if not isinstance(diag, dict):
        raise ValueError("MISSING_TERM_DIAGNOSTICS")
    values = {}
    for key in REQUIRED_TERM_DIAGNOSTICS:
        if key not in diag:
            raise ValueError("MISSING_TERM_DIAGNOSTIC_" + key.upper())
        try:
            value = float(diag[key])
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("INVALID_TERM_DIAGNOSTIC_" + key.upper()) from exc
        if not math.isfinite(value) or value < 0.0:
            raise ValueError("INVALID_TERM_DIAGNOSTIC_" + key.upper())
        values[key] = value
    for name, value in (("ACCEPTED2", accepted2), ("REJECTED2", rejected2)):
        try:
            numeric = float(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("INVALID_PROJECTION_NORM_" + name) from exc
        if not math.isfinite(numeric) or numeric < 0.0:
            raise ValueError("INVALID_PROJECTION_NORM_" + name)
        values[name.lower()] = numeric
    return values
