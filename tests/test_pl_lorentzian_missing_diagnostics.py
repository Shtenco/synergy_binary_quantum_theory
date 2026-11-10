import sys
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import peter_weyl_lorentzian_logical_projection_gate as LP
import pl_lorentzian_triple_worker as WORKER


class MissingLeakageDiagnosticTests(unittest.TestCase):
    def test_aggregate_physical_gate_rejects_missing_charge_diagnostic(self):
        partial = {
            "CV_complete_basis_leakage": 0.0,
            "CK_outer_complete_basis_leakage": 0.0,
            "CK_internal_volume_sector_leakage": 0.0,
        }
        result = LP.leakage_summary(partial)
        self.assertFalse(result["complete"])
        self.assertIsNone(result["maximum"])
        self.assertEqual(result["missing"], ["CK_complete_charge_basis_leakage"])

    def test_aggregate_physical_gate_accepts_four_explicit_diagnostics(self):
        complete = {
            "CV_complete_basis_leakage": 1e-12,
            "CK_outer_complete_basis_leakage": 2e-12,
            "CK_internal_volume_sector_leakage": 3e-12,
            "CK_complete_charge_basis_leakage": 4e-12,
        }
        result = LP.leakage_summary(complete)
        self.assertTrue(result["complete"])
        self.assertEqual(result["maximum"], 4e-12)

    def test_one_term_missing_charge_cannot_be_masked_by_another_term(self):
        complete = {
            "CV_complete_basis_leakage": 1e-12,
            "CK_outer_complete_basis_leakage": 2e-12,
            "CK_internal_volume_sector_leakage": 3e-12,
            "CK_complete_charge_basis_leakage": 4e-12,
        }
        partial = dict(complete)
        partial.pop("CK_complete_charge_basis_leakage")
        aggregate = {}
        missing_ever = set()
        invalid_ever = set()
        LP.merge_leakage_max(aggregate, missing_ever, invalid_ever, complete)
        LP.merge_leakage_max(aggregate, missing_ever, invalid_ever, partial)
        for key in missing_ever | invalid_ever:
            aggregate.pop(key, None)
        result = LP.leakage_summary(aggregate)
        self.assertFalse(result["complete"])
        self.assertIn("CK_complete_charge_basis_leakage", result["missing"])

    def test_nan_inf_negative_and_boolean_leakage_are_rejected_before_max(self):
        for bad in (float("nan"), float("inf"), -1e-12, True):
            with self.subTest(value=bad):
                with self.assertRaises(ValueError):
                    LP.update_diag({}, "CK_complete_charge_basis_leakage", bad)
                with self.assertRaises(ValueError):
                    WORKER.update_diag({}, {"CK_complete_charge_basis_leakage": bad})

    def test_collector_requires_all_four_numeric_leakage_channels(self):
        meta = {
            "leakage_diagnostics": {
                "CV_complete_basis_leakage": 0.0,
                "CK_outer_complete_basis_leakage": 0.0,
                "CK_internal_volume_sector_leakage": 0.0,
                "CK_complete_charge_basis_leakage": 0.0,
            },
            "physical_acceptance_max_leakage": 0.0,
            "complete_charge_diagnostic": 0.0,
            "exact_empty_input_zero_fields": [],
            "exact_zero_ordered_term": False,
        }
        self.assertEqual(BASE.validate_worker_leakage_provenance(meta, "test"), [])
        del meta["leakage_diagnostics"]["CK_complete_charge_basis_leakage"]
        with self.assertRaisesRegex(RuntimeError, "complete four-channel"):
            BASE.validate_worker_leakage_provenance(meta, "missing-charge")

    def test_collector_rejects_unsubstantiated_analytic_zero(self):
        meta = {
            "leakage_diagnostics": {
                "CV_complete_basis_leakage": 0.0,
                "CK_outer_complete_basis_leakage": 0.0,
                "CK_internal_volume_sector_leakage": 0.0,
                "CK_complete_charge_basis_leakage": 0.0,
            },
            "physical_acceptance_max_leakage": 0.0,
            "complete_charge_diagnostic": 0.0,
            "exact_empty_input_zero_fields": ["CV_complete_basis_leakage"],
            "exact_zero_ordered_term": False,
        }
        with self.assertRaisesRegex(RuntimeError, "nonzero term"):
            BASE.validate_worker_leakage_provenance(meta, "false-empty-claim")

    def test_forward_stage_not_run_is_not_reported_as_zero(self):
        with patch.object(LP.CV, "gauss_to_covariant", return_value={"seed": 1 + 0j}), \
             patch.object(LP.RAW.COMP, "C_volume_component", return_value=({}, 0.0)):
            state, diag = LP.ordered_triple_state("seed", 0, 1, 2, 3)

        self.assertEqual(state, {})
        self.assertIn("CV_complete_basis_leakage", diag)
        self.assertNotIn("CK_outer_complete_basis_leakage", diag)
        self.assertNotIn("CK_internal_volume_sector_leakage", diag)
        self.assertIn("CK_complete_charge_basis_leakage", diag)
        self.assertEqual(set(diag["exact_empty_input_zero_fields"]), {
            "CK_outer_complete_basis_leakage",
            "CK_internal_volume_sector_leakage",
            "CK_complete_charge_basis_leakage",
        })

    def test_permutation_ledger_retains_empty_input_zero_provenance(self):
        fake_diag = {
            "CV_complete_basis_leakage": 0.0,
            "CK_outer_complete_basis_leakage": 0.0,
            "CK_internal_volume_sector_leakage": 0.0,
            "CK_complete_charge_basis_leakage": 0.0,
            "exact_empty_input_zero_fields": [
                "CK_outer_complete_basis_leakage",
                "CK_internal_volume_sector_leakage",
                "CK_complete_charge_basis_leakage",
            ],
        }
        with patch.object(LP.PW, "NEIG", {0: (1, 2, 3, 4)}), \
             patch.object(LP, "ordered_triple_state", return_value=({}, fake_diag)):
            _, rows, summary = LP.epsilon_sum_state(("seed",), 0)

        self.assertEqual(len(rows), 24)
        self.assertTrue(all(
            set(row["exact_empty_input_zero_fields"]) == set(fake_diag["exact_empty_input_zero_fields"])
            for row in rows
        ))
        self.assertTrue(LP.leakage_summary(summary)["complete"])

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
        self.assertIn("CV_complete_basis_leakage", diag)
        self.assertEqual(diag["exact_empty_input_zero_fields"], ["CV_complete_basis_leakage"])


if __name__ == "__main__":
    unittest.main()
