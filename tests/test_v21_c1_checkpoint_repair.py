#!/usr/bin/env python3
"""V21 structural regressions. Synthetic states test software only, never physics."""
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from v21_stage_logic import all_32_specs, first_leg_census, operations, predecessor_ok
from v21_checkpoint_guard import validate_physical_stage
from pl_v21_pair_collector import compare_states, read_stage, read_v2
from pl_v20_stage_worker import save_charged, load_charged


class V21C1Tests(unittest.TestCase):
    def setUp(self):
        self.st = {((1,) * 32, (0,) * 16, 0, 0, 1, 1): 1 + 2j}
        self.one = {
            "stage": 1, "passed": True, "source_commit": "a" * 40,
            "source_node": 0, "omit": 0, "cycle": 0, "mode": "forward",
            "order": "p", "aux": [0, 0, 1], "Jmax2": 7,
            "volume_backend": "prod", "volume_backend_sha256": "b" * 64,
            "seed_state_sha256": "c" * 64, "input_state_sha256": "c" * 64,
            "output_state_sha256": "d" * 64, "npz_sha256": "d" * 64,
            "guard": {"stage_physical_basis_pass": True},
            "operator_signature": ["CV", 4, 1, 0],
        }
        self.two = dict(self.one, stage=2, input_state_sha256="d" * 64,
                        output_state_sha256="e" * 64, npz_sha256="e" * 64,
                        order="q", aux=[0, 1, 1],
                        operator_signature=["CK", 2, 1, 1],
                        post_second_scalar_filter=True)
        self.three = dict(self.two, stage=3, input_state_sha256="e" * 64,
                          output_state_sha256="f" * 64, npz_sha256="f" * 64,
                          operator_signature=["CK", 3, 0, 1])

    def test_01_census32(self):
        self.assertEqual(len(list(all_32_specs((2, 3, 4)))), 32)
        self.assertEqual(first_leg_census((2, 3, 4)), (12, 32))

    def test_02_original_forward(self):
        self.assertEqual(operations((2, 3, 4), "forward", (0, 1, 1)),
                         [("CV", 4, 1, 0), ("CK", 3, 1, 1), ("CK", 2, 0, 1)])

    def test_03_original_adjoint(self):
        self.assertEqual(operations((2, 3, 4), "adjoint", (0, 1, 1)),
                         [("CK", 2, 1, 0), ("CK", 3, 1, 1), ("CV", 4, 0, 1)])

    def test_04_safe_reuse_first_with_direct_hash(self):
        self.assertTrue(predecessor_ok(self.one, self.two, 2, ["CV", 4, 1, 0]))

    def test_05_reject_wrong_first(self):
        with self.assertRaisesRegex(ValueError, "MISMATCH_FIRST_OPERATOR_SIGNATURE"):
            predecessor_ok(self.one, self.two, 2, ["CV", 4, 0, 1])

    def test_06_stage3_requires_aux(self):
        with self.assertRaisesRegex(ValueError, "MISMATCH_aux"):
            predecessor_ok(self.two, dict(self.three, aux=[1, 1, 1]),
                           3, None, ["CK", 2, 1, 1])

    def test_07_requires_scalar_filter(self):
        with self.assertRaisesRegex(ValueError, "STAGE2_FILTER_NOT_APPLIED"):
            predecessor_ok(dict(self.two, post_second_scalar_filter=False),
                           self.three, 3, None, ["CK", 2, 1, 1])

    def test_08_reject_backend_mismatch(self):
        with self.assertRaisesRegex(ValueError, "MISMATCH_volume_backend_sha256"):
            predecessor_ok(self.one, dict(self.two, volume_backend_sha256="f" * 64),
                           2, ["CV", 4, 1, 0])

    def test_09_reject_broken_direct_checkpoint_chain(self):
        with self.assertRaisesRegex(ValueError, "BROKEN_DIRECT_CHECKPOINT_CHAIN"):
            predecessor_ok(self.one, dict(self.two, input_state_sha256="9" * 64),
                           2, ["CV", 4, 1, 0])

    def test_10_reject_output_npz_digest_disagreement(self):
        with self.assertRaisesRegex(ValueError, "PREDECESSOR_OUTPUT_DIGEST_MISMATCH"):
            predecessor_ok(dict(self.one, npz_sha256="9" * 64), self.two,
                           2, ["CV", 4, 1, 0])

    def test_11_cv_requires_cv_diagnostic(self):
        guard = validate_physical_stage(
            {"CV_complete_basis_leakage": 1e-14}, self.st, 1, ["CV", 4, 1, 0])
        self.assertTrue(guard["stage_physical_basis_pass"])

    def test_12_ck_requires_both_physical_and_charge_diagnostics(self):
        diag = {"outer_complete_basis_leakage": 0.0,
                "internal_volume_sector_leakage": 0.0,
                "complete_charge_basis_leakage": 1.0}
        guard = validate_physical_stage(diag, self.st, 2, ["CK", 2, 1, 1])
        self.assertEqual(guard["charge_leakage_max_unprojected"], 1.0)
        self.assertFalse(guard["charge_covariance_verified"])

    def test_13_ck_missing_internal_diagnostic_rejected(self):
        with self.assertRaisesRegex(ValueError, "MISSING_OPERATION_SPECIFIC_DIAGNOSTIC"):
            validate_physical_stage(
                {"outer_complete_basis_leakage": 0.0,
                 "complete_charge_basis_leakage": 0.0},
                self.st, 2, ["CK", 2, 1, 1])

    def test_14_charge_diagnostic_alone_not_cv_pass(self):
        with self.assertRaisesRegex(ValueError, "MISSING_OPERATION_SPECIFIC_DIAGNOSTIC"):
            validate_physical_stage(
                {"complete_charge_basis_leakage": 0.0}, self.st, 1, ["CV", 4, 1, 0])

    def test_15_high_physical_leak_rejected(self):
        with self.assertRaisesRegex(ValueError, "PHYSICAL_BASIS_LEAKAGE"):
            validate_physical_stage(
                {"outer_complete_basis_leakage": 1.0,
                 "internal_volume_sector_leakage": 0.0,
                 "complete_charge_basis_leakage": 0.0},
                self.st, 2, ["CK", 2, 1, 1])

    def test_16_empty_stage_is_not_pass(self):
        with self.assertRaisesRegex(ValueError, "EXACT_ZERO_STAGE_OUTPUT"):
            validate_physical_stage(
                {"CV_complete_basis_leakage": 0.0}, {}, 1, ["CV", 4, 1, 0])

    def test_17_zero_amplitude_support_is_not_pass(self):
        zero = {next(iter(self.st)): 0j}
        with self.assertRaisesRegex(ValueError, "EXACT_ZERO_STAGE_OUTPUT"):
            validate_physical_stage(
                {"CV_complete_basis_leakage": 0.0}, zero, 1, ["CV", 4, 1, 0])

    def test_18_nonfinite_diagnostic_rejected(self):
        with self.assertRaisesRegex(ValueError, "NONFINITE_OR_BAD_LEAKAGE"):
            validate_physical_stage(
                {"CV_complete_basis_leakage": float("nan")},
                self.st, 1, ["CV", 4, 1, 0])

    def test_19_complex_phase_matters(self):
        key = next(iter(self.st))
        self.assertGreater(compare_states({key: 1 + 0j}, {key: 0 + 1j})["relative_delta"], 1.4)

    def test_20_checkpoint_integrity(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "x.npz"
            digest = save_charged(path, self.st)
            self.assertEqual(load_charged(path, digest), self.st)
            with self.assertRaises(ValueError):
                load_charged(path, "0" * 64)

    def test_21_collector_rejects_stored_guard_without_diagnostics(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "stage"
            digest = save_charged(prefix.with_suffix(".npz"), self.st)
            prefix.with_suffix(".json").write_text(json.dumps(
                dict(self.one, passed=True, status="REAL_ORIGINAL_STAGE_ONLY",
                     npz_sha256=digest, output_state_sha256=digest,
                     state_support=1, stage=1,
                     guard={"stage_physical_basis_pass": True})))
            with self.assertRaisesRegex(ValueError, "MISSING_STAGE_DIAGNOSTICS"):
                read_stage(prefix)

    def test_22_collector_recomputes_guard_and_checks_output_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "stage"
            digest = save_charged(prefix.with_suffix(".npz"), self.st)
            diag = {"CV_complete_basis_leakage": 0.0}
            guard = validate_physical_stage(diag, self.st, 1, ["CV", 4, 1, 0])
            meta = dict(self.one, passed=True, status="REAL_ORIGINAL_STAGE_ONLY",
                        npz_sha256=digest, output_state_sha256=digest,
                        state_support=1, stage=1, diagnostics=diag,
                        operator_signature=["CV", 4, 1, 0], guard=guard)
            prefix.with_suffix(".json").write_text(json.dumps(meta))
            self.assertEqual(read_stage(prefix)[0]["npz_sha256"], digest)
            meta["guard"]["physical_leakage_max"] = 0.5
            prefix.with_suffix(".json").write_text(json.dumps(meta))
            with self.assertRaisesRegex(ValueError, "STORED_GUARD_RECOMPUTATION_MISMATCH"):
                read_stage(prefix)

    def test_23_read_original_v2(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "v2.npz"
            np.savez_compressed(path, spins=np.array([[1] * 32]),
                                Ks=np.array([[0] * 16]), amp=np.array([1 + 2j]))
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(list(read_v2({"npz": str(path), "sha256": digest}).values()), [1 + 2j])

    def test_24_reject_v2_corruption(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "v2.npz"
            np.savez_compressed(path, spins=np.empty((0, 32)),
                                Ks=np.empty((0, 16)), amp=np.empty((0,), dtype=complex))
            with self.assertRaisesRegex(ValueError, "INDEPENDENT_V2_SHA256_MISMATCH"):
                read_v2({"npz": str(path), "sha256": "0" * 64})


if __name__ == "__main__":
    unittest.main()
