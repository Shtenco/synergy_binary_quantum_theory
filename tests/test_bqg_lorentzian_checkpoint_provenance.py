from __future__ import annotations

import ast
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import bqg_lorentzian_checkpointed_triple as BASE
import bqg_lorentzian_ka_partition as PART


class LorentzianCheckpointProvenanceTests(unittest.TestCase):
    def test_source_hashes_include_transitive_operator_dependencies(self):
        hashes = BASE.source_hashes()
        required = {
            "bqg_lorentzian_checkpointed_triple.py",
            "bqg_lorentzian_ka_partition.py",
            "peter_weyl_lorentzian_K_block_gate.py",
            "peter_weyl_covariant_composition_gate.py",
            "peter_weyl_covariant_K_composition_gate.py",
            "peter_weyl_covariant_K_leg_gate.py",
            "peter_weyl_covariant_K_projection_audit_gate.py",
            "peter_weyl_euclidean_sine_ordering_gate.py",
            "charged_intertwiner_recoupling_gate.py",
        }
        self.assertTrue(required.issubset(hashes), sorted(required - set(hashes)))
        self.assertTrue(all(len(value) == 64 for value in hashes.values()))

    def test_expected_freeze_binds_schema_runtime_and_source(self):
        args = SimpleNamespace(
            source=1, column_index=0, a=0, b=2, c=3, coefficient=-1
        )
        freeze = PART.expected_freeze(args)
        self.assertEqual(freeze["checkpoint_schema"], BASE.CHECKPOINT_SCHEMA)
        self.assertEqual(freeze["runtime"], BASE.runtime_fingerprint())
        self.assertEqual(freeze["modules"], BASE.source_hashes())
        self.assertIn("bqg_lorentzian_ka_partition.py", freeze["modules"])

    def test_partition_key_digest_binds_membership_and_order(self):
        keys = [((1, 2), (0, 0)), ((1, 3), (0, 2))]
        self.assertNotEqual(PART.key_digest(keys), PART.key_digest(keys[:1]))
        self.assertNotEqual(PART.key_digest(keys), PART.key_digest(list(reversed(keys))))

    def test_partial_stage_diagnostics_are_preserved_without_zero_filling(self):
        partial = {"CV_complete_basis_leakage": 0.0}
        self.assertEqual(
            BASE.validate_leakage_map(partial, require_all=False),
            partial,
        )
        with self.assertRaisesRegex(ValueError, "missing required leakage"):
            BASE.validate_leakage_map(partial, require_all=True)

    def test_finish_requires_all_four_diagnostics(self):
        complete = {key: 0.0 for key in BASE.LEAK_FIELDS}
        self.assertEqual(
            BASE.validate_leakage_map(complete, require_all=True),
            complete,
        )
        complete.pop("CK_complete_charge_basis_leakage")
        with self.assertRaisesRegex(ValueError, "CK_complete_charge_basis_leakage"):
            BASE.validate_leakage_map(complete, require_all=True)

    def test_merging_empty_stages_does_not_synthesize_zero_diagnostics(self):
        merged = BASE.maxdiag({}, {})
        self.assertEqual(merged, {})
        with self.assertRaisesRegex(ValueError, "missing required leakage"):
            BASE.validate_leakage_map(merged, require_all=True)

    def test_exact_empty_input_zero_requires_explicit_provenance(self):
        meta = {
            "operator_input_support": 0,
            "operator_stage_applied": False,
            "exact_empty_input_zero_fields": [
                "CK_outer_complete_basis_leakage",
                "CK_internal_volume_sector_leakage",
                "CK_complete_charge_basis_leakage",
            ],
            "diagnostic_provenance_scope": "current_operator_stage",
            "stage_leakage_diagnostics": {
                "CK_outer_complete_basis_leakage": 0.0,
                "CK_internal_volume_sector_leakage": 0.0,
                "CK_complete_charge_basis_leakage": 0.0,
            },
        }
        result = BASE.validate_operator_provenance(meta)
        self.assertFalse(result["operator_stage_applied"])
        self.assertEqual(len(result["exact_empty_input_zero_fields"]), 3)

    def test_empty_input_provenance_cannot_be_claimed_for_nonempty_input(self):
        meta = {
            "operator_input_support": 4,
            "operator_stage_applied": True,
            "exact_empty_input_zero_fields": [
                "CK_outer_complete_basis_leakage",
                "CK_internal_volume_sector_leakage",
                "CK_complete_charge_basis_leakage",
            ],
            "diagnostic_provenance_scope": "current_operator_stage",
            "stage_leakage_diagnostics": {
                "CK_outer_complete_basis_leakage": 0.0,
                "CK_internal_volume_sector_leakage": 0.0,
                "CK_complete_charge_basis_leakage": 0.0,
            },
        }
        with self.assertRaisesRegex(ValueError, "provenance mismatch"):
            BASE.validate_operator_provenance(meta)

    def test_operator_skipped_flag_must_match_input_support(self):
        meta = {
            "operator_input_support": 0,
            "operator_stage_applied": True,
            "exact_empty_input_zero_fields": [
                "CK_outer_complete_basis_leakage",
                "CK_internal_volume_sector_leakage",
                "CK_complete_charge_basis_leakage",
            ],
            "diagnostic_provenance_scope": "current_operator_stage",
            "stage_leakage_diagnostics": {
                "CK_outer_complete_basis_leakage": 0.0,
                "CK_internal_volume_sector_leakage": 0.0,
                "CK_complete_charge_basis_leakage": 0.0,
            },
        }
        with self.assertRaisesRegex(ValueError, "execution flag mismatch"):
            BASE.validate_operator_provenance(meta)

    def test_missing_nan_and_negative_diagnostics_fail_closed(self):
        for value in (float("nan"), float("inf"), -1.0, True, "0"):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    BASE.validate_leakage_map(
                        {"CV_complete_basis_leakage": value},
                        require_all=False,
                    )

    def test_partition_worker_has_no_optimizable_assert_guards(self):
        tree = ast.parse(
            Path(PART.__file__).read_text(encoding="utf-8"),
            filename=PART.__file__,
        )
        self.assertFalse(any(isinstance(node, ast.Assert) for node in ast.walk(tree)))


if __name__ == "__main__":
    unittest.main()
