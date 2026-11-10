import math
import unittest


REQUIRED_TERM_DIAGNOSTICS = (
    "CV_complete_basis_leakage",
    "CK_outer_complete_basis_leakage",
    "CK_internal_volume_sector_leakage",
    "CK_complete_charge_basis_leakage",
)


def validate_term_diagnostics(diag, accepted2, rejected2):
    """Isolated copy of the V2 worker's strict diagnostic validator."""
    if not isinstance(diag, dict):
        raise ValueError("MISSING_TERM_DIAGNOSTICS")
    values = {}
    for key in REQUIRED_TERM_DIAGNOSTICS:
        if key not in diag:
            raise ValueError("MISSING_TERM_DIAGNOSTIC_" + key.upper())
        raw = diag[key]
        if isinstance(raw, bool) or not isinstance(raw, (int, float)):
            raise ValueError("INVALID_TERM_DIAGNOSTIC_" + key.upper())
        value = float(raw)
        if not math.isfinite(value) or value < 0.0:
            raise ValueError("INVALID_TERM_DIAGNOSTIC_" + key.upper())
        values[key] = value
    for name, value in (("ACCEPTED2", accepted2), ("REJECTED2", rejected2)):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError("INVALID_PROJECTION_NORM_" + name)
        numeric = float(value)
        if not math.isfinite(numeric) or numeric < 0.0:
            raise ValueError("INVALID_PROJECTION_NORM_" + name)
        values[name.lower()] = numeric
    return values


class V2DiagnosticGuardAudit(unittest.TestCase):
    def setUp(self):
        self.good = {
            "CV_complete_basis_leakage": 0.0,
            "CK_outer_complete_basis_leakage": 0.0,
            "CK_internal_volume_sector_leakage": 0.0,
            "CK_complete_charge_basis_leakage": 0.0,
        }

    def test_missing_required_diagnostic_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "MISSING_TERM_DIAGNOSTIC_CV_COMPLETE_BASIS_LEAKAGE"):
            validate_term_diagnostics({}, 1.0, 0.0)

    def test_negative_physical_leakage_is_rejected(self):
        diag = dict(self.good, CV_complete_basis_leakage=-1e-9)
        with self.assertRaisesRegex(ValueError, "INVALID_TERM_DIAGNOSTIC_CV_COMPLETE_BASIS_LEAKAGE"):
            validate_term_diagnostics(diag, 1.0, 0.0)

    def test_nonfinite_charge_diagnostic_is_rejected(self):
        diag = dict(self.good, CK_complete_charge_basis_leakage=float("nan"))
        with self.assertRaisesRegex(ValueError, "INVALID_TERM_DIAGNOSTIC_CK_COMPLETE_CHARGE_BASIS_LEAKAGE"):
            validate_term_diagnostics(diag, 1.0, 0.0)

    def test_negative_projection_norm_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "INVALID_PROJECTION_NORM_REJECTED2"):
            validate_term_diagnostics(self.good, 1.0, -1.0)

    def test_boolean_leakage_diagnostic_is_rejected(self):
        diag = dict(self.good, CV_complete_basis_leakage=True)
        with self.assertRaisesRegex(ValueError, "INVALID_TERM_DIAGNOSTIC_CV_COMPLETE_BASIS_LEAKAGE"):
            validate_term_diagnostics(diag, 1.0, 0.0)

    def test_numeric_string_leakage_diagnostic_is_rejected(self):
        diag = dict(self.good, CV_complete_basis_leakage="0.0")
        with self.assertRaisesRegex(ValueError, "INVALID_TERM_DIAGNOSTIC_CV_COMPLETE_BASIS_LEAKAGE"):
            validate_term_diagnostics(diag, 1.0, 0.0)

    def test_boolean_projection_norm_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "INVALID_PROJECTION_NORM_ACCEPTED2"):
            validate_term_diagnostics(self.good, True, 0.0)

    def test_nonfinite_projection_norm_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "INVALID_PROJECTION_NORM_ACCEPTED2"):
            validate_term_diagnostics(self.good, float("inf"), 0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
