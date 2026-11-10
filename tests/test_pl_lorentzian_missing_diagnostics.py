import sys
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import peter_weyl_lorentzian_logical_projection_gate as LP
import pl_lorentzian_triple_worker as WORKER


class MissingLeakageDiagnosticTests(unittest.TestCase):
    def test_forward_stage_not_run_is_not_reported_as_zero(self):
        with patch.object(LP.CV, "gauss_to_covariant", return_value={"seed": 1 + 0j}), \
             patch.object(LP.RAW.COMP, "C_volume_component", return_value=({}, 0.0)):
            state, diag = LP.ordered_triple_state("seed", 0, 1, 2, 3)

        self.assertEqual(state, {})
        self.assertIn("CV_complete_basis_leakage", diag)
        self.assertNotIn("CK_outer_complete_basis_leakage", diag)
        self.assertNotIn("CK_internal_volume_sector_leakage", diag)
        self.assertNotIn("CK_complete_charge_basis_leakage", diag)

    def test_dagger_volume_stage_not_run_is_not_reported_as_zero(self):
        ck_diag = {
            "outer_complete_basis_leakage": 0.0,
            "internal_volume_sector_leakage": 0.0,
            "complete_charge_basis_leakage": 0.0,
        }
        with patch.object(LP.CV, "gauss_to_covariant", return_value={"seed": 1 + 0j}), \
             patch.object(LP.RAW.KCOMP, "C_K_component", return_value=({}, ck_diag)):
            state, diag = WORKER.ordered_dagger_state("seed", 0, 1, 2, 3)

        self.assertEqual(state, {})
        self.assertIn("CK_outer_complete_basis_leakage", diag)
        self.assertIn("CK_internal_volume_sector_leakage", diag)
        self.assertIn("CK_complete_charge_basis_leakage", diag)
        self.assertNotIn("CV_complete_basis_leakage", diag)


if __name__ == "__main__":
    unittest.main()
