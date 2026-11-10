#!/usr/bin/env python3
"""V21 structural regressions. Synthetic states test software only, never physics."""
import hashlib
import json
import math
from random import Random
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from v21_stage_logic import all_32_specs, first_leg_census, operations, predecessor_ok
from v21_checkpoint_guard import validate_physical_stage
from pl_v21_pair_collector import compare_states, read_stage, read_v2, assemble_pair_exact, sum_states_exact, sum_states_legacy_pruned, validate_path_metadata, validate_manifest_indices
from pl_v20_stage_worker import save_charged, load_charged
from v21_stage_attestation import attest_manifest
from pl_lorentzian_48_collect_v2 import validate_term_record
from v2_term_guards import validate_term_diagnostics


class V21C1Tests(unittest.TestCase):
    def setUp(self):
        self.st = {((1,) * 32, (0,) * 16, 0, 0, 1, 1): 1 + 2j}
        self.one = {
            "stage": 1, "passed": True, "status": "REAL_ORIGINAL_STAGE_ONLY", "state_support": 1, "source_commit": "a" * 40,
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
            with self.assertRaisesRegex(ValueError, "MISSING_EXTERNAL_STAGE_ATTESTATION"):
                read_stage(prefix, None)

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
            meta_path = prefix.with_suffix(".json")
            meta_path.write_text(json.dumps(meta))
            attestation = {
                "metadata_sha256": hashlib.sha256(meta_path.read_bytes()).hexdigest(),
                "npz_sha256": digest,
            }
            self.assertEqual(read_stage(prefix, attestation)[0]["npz_sha256"], digest)
            meta["guard"]["physical_leakage_max"] = 0.5
            meta_path.write_text(json.dumps(meta))
            attestation["metadata_sha256"] = hashlib.sha256(meta_path.read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError, "STORED_GUARD_RECOMPUTATION_MISMATCH"):
                read_stage(prefix, attestation)

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

    def test_25_attestation_builder_pins_both_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            item = {}
            for stage in (1, 2, 3):
                prefix = root / f"stage{stage}"
                digest = save_charged(prefix.with_suffix(".npz"), self.st)
                meta = {"stage": stage, "passed": True, "npz_sha256": digest,
                        "output_state_sha256": digest}
                prefix.with_suffix(".json").write_text(json.dumps(meta))
                item[f"stage{stage}"] = str(prefix)
            result = attest_manifest({"paths": [item]})
            pins = result["paths"][0]["stage_attestations"]
            self.assertEqual(set(pins), {"stage1", "stage2", "stage3"})
            self.assertEqual(pins["stage2"]["npz_sha256"],
                             hashlib.sha256(Path(item["stage2"]).with_suffix(".npz").read_bytes()).hexdigest())

    def test_26_attestation_builder_rejects_npz_tampering(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            item = {}
            for stage in (1, 2, 3):
                prefix = root / f"stage{stage}"
                digest = save_charged(prefix.with_suffix(".npz"), self.st)
                meta = {"stage": stage, "passed": True, "npz_sha256": digest,
                        "output_state_sha256": digest}
                prefix.with_suffix(".json").write_text(json.dumps(meta))
                item[f"stage{stage}"] = str(prefix)
            path = Path(item["stage2"]).with_suffix(".npz")
            path.write_bytes(path.read_bytes() + b"tamper")
            with self.assertRaisesRegex(ValueError, "NPZ_DIGEST_MISMATCH_STAGE2"):
                attest_manifest({"paths": [item]})

    def test_27_exact_pair_cancellation_is_distinguished(self):
        key = next(iter(self.st))
        pair, norm, status = assemble_pair_exact(
            {"fp": {key: 1 + 2j}, "fq": {key: 1 + 2j},
             "ap": {}, "aq": {}},
            {"fp": 1.0, "fq": -1.0, "ap": 0.0, "aq": 0.0})
        self.assertEqual(norm, 0.0)
        self.assertNotIn(key, pair)
        self.assertEqual(status, "VERIFIED_EXACT_ZERO_PAIR")

    def test_28_small_nonzero_pair_is_not_pruned(self):
        key = next(iter(self.st))
        pair, norm, status = assemble_pair_exact(
            {"fp": {key: 1e-14 + 0j}, "fq": {}, "ap": {}, "aq": {}},
            {"fp": 1.0, "fq": 0.0, "ap": 0.0, "aq": 0.0})
        self.assertGreater(norm, 0.0)
        self.assertEqual(pair[key], 1e-14 + 0j)
        self.assertEqual(status, "VERIFIED_NONZERO_PAIR")

    def test_29_reject_non_original_predecessor_status(self):
        with self.assertRaisesRegex(ValueError, "STAGE_NOT_PASSED"):
            predecessor_ok(dict(self.one, status="MOCK_OR_SYNTHETIC"),
                           self.two, 2, ["CV", 4, 1, 0])

    def test_30_auxiliary_sum_preserves_subthreshold_amplitudes(self):
        key = next(iter(self.st))
        total = sum_states_exact([{key: 1e-14 + 0j}, {}, {key: 2e-14 + 0j}])
        self.assertEqual(total[key], 3e-14 + 0j)

    def test_31_auxiliary_sum_removes_only_exact_cancellation(self):
        key = next(iter(self.st))
        total = sum_states_exact([{key: 1 + 2j}, {key: -1 - 2j}])
        self.assertNotIn(key, total)

    def test_32_stable_sum_recovers_small_residual_after_large_cancellation(self):
        key = next(iter(self.st))
        total = sum_states_exact([
            {key: 1e16 + 0j}, {key: 1.0 + 0j}, {key: -1e16 + 0j}
        ])
        self.assertEqual(total[key], 1.0 + 0j)

    def test_33_pair_assembly_is_order_stable_for_cancellation(self):
        key = next(iter(self.st))
        pair, norm, status = assemble_pair_exact(
            {"fp": {key: 1e16 + 0j}, "fq": {key: 1.0 + 0j},
             "ap": {key: -1e16 + 0j}, "aq": {}},
            {"fp": 1.0, "fq": 1.0, "ap": 1.0, "aq": 0.0})
        self.assertEqual(pair[key], 1.0 + 0j)
        self.assertEqual(norm, 1.0)
        self.assertEqual(status, "VERIFIED_NONZERO_PAIR")

    def test_34_legacy_sum_reproduces_ordered_pruning(self):
        key = next(iter(self.st))
        contributions = [{key: 0.75e-11 + 0j}, {key: 0.75e-11 + 0j}]
        legacy = sum_states_legacy_pruned(contributions, 1e-11)
        unpruned = sum_states_exact(contributions)
        # Historical add() drops each sub-threshold partial sum before the
        # next auxiliary path arrives; the unpruned sum retains their total.
        self.assertNotIn(key, legacy)
        self.assertEqual(unpruned[key], 1.5e-11 + 0j)

    def test_35_collector_ignores_attacker_controlled_guard_threshold(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "stage"
            digest = save_charged(prefix.with_suffix(".npz"), self.st)
            diagnostics = {"CV_complete_basis_leakage": 1e-6}
            permissive_guard = validate_physical_stage(
                diagnostics, self.st, 1, ["CV", 4, 1, 0], threshold=1.0)
            meta = dict(self.one, passed=True, status="REAL_ORIGINAL_STAGE_ONLY",
                        npz_sha256=digest, output_state_sha256=digest,
                        state_support=1, stage=1, diagnostics=diagnostics,
                        operator_signature=["CV", 4, 1, 0], guard=permissive_guard)
            meta_path = prefix.with_suffix(".json")
            meta_path.write_text(json.dumps(meta))
            attestation = {
                "metadata_sha256": hashlib.sha256(meta_path.read_bytes()).hexdigest(),
                "npz_sha256": digest,
            }
            with self.assertRaisesRegex(ValueError, "PHYSICAL_BASIS_LEAKAGE_EXCEEDS"):
                read_stage(prefix, attestation)

    def test_36_tiny_nonzero_stage_is_not_misclassified_as_zero(self):
        tiny = {next(iter(self.st)): 1e-200 + 0j}
        guard = validate_physical_stage(
            {"CV_complete_basis_leakage": 0.0}, tiny, 1, ["CV", 4, 1, 0])
        self.assertTrue(guard["stage_physical_basis_pass"])

    def test_37_tiny_nonzero_pair_is_not_misclassified_by_squared_norm_underflow(self):
        key = next(iter(self.st))
        pair, norm, status = assemble_pair_exact(
            {"fp": {key: 1e-200 + 0j}, "fq": {}, "ap": {}, "aq": {}},
            {"fp": 1.0, "fq": 0.0, "ap": 0.0, "aq": 0.0})
        self.assertEqual(status, "VERIFIED_NONZERO_PAIR")
        self.assertEqual(pair[key], 1e-200 + 0j)
        self.assertEqual(norm, 1e-200)

    def test_38_relative_comparison_does_not_hide_tiny_state_difference(self):
        key = next(iter(self.st))
        result = compare_states({key: 1e-200 + 0j}, {})
        self.assertEqual(result["relative_delta"], 1.0)
        self.assertEqual(result["absolute_delta"], 1e-200)

    def test_39_overflowed_relative_delta_is_not_accepted_as_nan(self):
        key = next(iter(self.st))
        result = compare_states({key: 1e308 + 0j}, {key: -1e308 + 0j})
        self.assertTrue(math.isinf(result["relative_delta"]))

    def test_40_v2_reference_must_match_attested_volume_backend(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            npz = root / "v2.npz"
            np.savez_compressed(npz, spins=np.array([[1] * 32]),
                                Ks=np.array([[0] * 16]), amp=np.array([1 + 0j]))
            metadata = {
                "passed": True, "mode": "forward", "index": 3,
                "operator_version": "tetrahedral-charged-volume-v2",
                "source_node": 0, "Jmax": 3.5, "checks": {"ok": True},
                "source_commit": "a" * 40,
                "volume_backend": "tetrahedral-volume-v1",
                "volume_backend_sha256": "c" * 64,
            }
            meta_path = root / "v2.json"
            meta_path.write_text(json.dumps(metadata))
            spec = {
                "npz": str(npz), "sha256": hashlib.sha256(npz.read_bytes()).hexdigest(),
                "metadata_json": str(meta_path),
                "metadata_sha256": hashlib.sha256(meta_path.read_bytes()).hexdigest(),
            }
            with self.assertRaisesRegex(ValueError, "V2_VOLUME_BACKEND_MISMATCH"):
                read_v2(spec, expected_source=0, expected_mode="forward",
                        expected_index=3, expected_source_commit="a" * 40,
                        expected_backend_sha="b" * 64)

    def test_41_attestation_builder_pins_four_v2_references(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = {
                "source_commit": "a" * 40,
                "volume_backend": "tetrahedral-volume-v1",
                "volume_backend_sha256": "b" * 64,
                "paths": [],
                "independent_v2": {},
            }
            for name in ("fp", "fq", "ap", "aq"):
                npz = root / f"{name}.npz"
                np.savez_compressed(npz, spins=np.array([[1] * 32]),
                                    Ks=np.array([[0] * 16]), amp=np.array([1 + 0j]))
                metadata = {
                    "passed": True, "checks": {
                "finite_covariant_norm": True, "finite_gauss_norm": True,
                "finite_covariant_amplitudes": True, "finite_gauss_amplitudes": True,
                "physical_basis_volume_leakage": True,
                "scalar_closure_fraction_or_exact_zero": True,
                "nonscalar_rejected_norm": True, "single_L_spin_wall": True,
                "PL_orientation_coefficient": True,
            },
                    "operator_version": "tetrahedral-charged-volume-v2",
                    "source_commit": "a" * 40,
                    "volume_backend": "tetrahedral-volume-v1",
                    "volume_backend_sha256": "b" * 64,
                }
                meta_path = root / f"{name}.json"
                meta_path.write_text(json.dumps(metadata))
                manifest["independent_v2"][name] = {
                    "npz": str(npz), "metadata_json": str(meta_path)
                }
            result = attest_manifest(manifest)
            for spec in result["independent_v2"].values():
                self.assertEqual(len(spec["sha256"]), 64)
                self.assertEqual(len(spec["metadata_sha256"]), 64)

    def test_42_attestation_builder_rejects_v2_backend_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = {
                "source_commit": "a" * 40,
                "volume_backend": "tetrahedral-volume-v1",
                "volume_backend_sha256": "b" * 64,
                "paths": [],
                "independent_v2": {},
            }
            for name in ("fp", "fq", "ap", "aq"):
                npz = root / f"{name}.npz"
                np.savez_compressed(npz, spins=np.array([[1] * 32]),
                                    Ks=np.array([[0] * 16]), amp=np.array([1 + 0j]))
                metadata = {
                    "passed": True, "checks": {"finite": True},
                    "operator_version": "tetrahedral-charged-volume-v2",
                    "source_commit": "a" * 40,
                    "volume_backend": ("legacy-volume" if name == "aq"
                                       else "tetrahedral-volume-v1"),
                    "volume_backend_sha256": "b" * 64,
                }
                meta_path = root / f"{name}.json"
                meta_path.write_text(json.dumps(metadata))
                manifest["independent_v2"][name] = {
                    "npz": str(npz), "metadata_json": str(meta_path)
                }
            with self.assertRaisesRegex(ValueError, "V2_REFERENCE_BACKEND_MISMATCH_aq"):
                attest_manifest(manifest)


    def test_43_v2_index_arrays_must_have_integer_dtype(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "v2-float-indices.npz"
            np.savez_compressed(path, spins=np.array([[1.5] * 32]),
                                Ks=np.array([[0] * 16]), amp=np.array([1 + 0j]))
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError, "V2_BAD_SHAPE_DTYPE_OR_NONFINITE"):
                read_v2({"npz": str(path), "sha256": digest})

    def test_44_v2_checks_must_be_boolean_true_not_truthy(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            npz = root / "v2.npz"
            np.savez_compressed(npz, spins=np.array([[1] * 32]),
                                Ks=np.array([[0] * 16]), amp=np.array([1 + 0j]))
            metadata = {
                "passed": True, "mode": "forward", "index": 3,
                "operator_version": "tetrahedral-charged-volume-v2",
                "source_node": 0, "Jmax": 3.5, "checks": {"ok": "yes"},
                "npz_sha256": hashlib.sha256(npz.read_bytes()).hexdigest(),
                "source_commit": "a" * 40,
                "volume_backend": "tetrahedral-volume-v1",
                "volume_backend_sha256": "c" * 64,
            }
            meta_path = root / "v2.json"
            meta_path.write_text(json.dumps(metadata))
            spec = {
                "npz": str(npz), "sha256": hashlib.sha256(npz.read_bytes()).hexdigest(),
                "metadata_json": str(meta_path),
                "metadata_sha256": hashlib.sha256(meta_path.read_bytes()).hexdigest(),
            }
            with self.assertRaisesRegex(ValueError, "V2_PRIMITIVE_FAILED_ACCEPTANCE_OR_INCOMPLETE_CHECKS"):
                read_v2(spec, expected_source=0, expected_mode="forward",
                        expected_index=3, expected_source_commit="a" * 40,
                        expected_backend_sha="c" * 64)

    def test_45_checkpoint_index_arrays_must_have_integer_dtype(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad-index.npz"
            np.savez_compressed(
                path, spins=np.array([[1.5] * 32]),
                K_other=np.array([[0] * 16]), charges=np.array([[0] * 4]),
                amp=np.array([1 + 0j]))
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError, "BAD_CHECKPOINT_SHAPE_OR_INDEX_DTYPE"):
                load_charged(path, digest)



    def test_46_stage_number_must_be_integer_not_boolean(self):
        with self.assertRaisesRegex(ValueError, "BAD_STAGE"):
            validate_physical_stage(
                {"CV_complete_basis_leakage": 0.0}, self.st, True,
                ["CV", 4, 1, 0])

    def test_47_operator_signature_indices_must_be_integer_not_boolean(self):
        with self.assertRaisesRegex(ValueError, "BAD_OPERATOR_SIGNATURE"):
            validate_physical_stage(
                {"CV_complete_basis_leakage": 0.0}, self.st, 1,
                ["CV", 4, True, 0])

    def test_48_v2_input_indices_must_be_integer_not_boolean(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            npz = root / "v2.npz"
            np.savez_compressed(npz, spins=np.array([[1] * 32]),
                                Ks=np.array([[0] * 16]), amp=np.array([1 + 0j]))
            metadata = {
                "passed": True, "mode": "forward", "index": False,
                "operator_version": "tetrahedral-charged-volume-v2",
                "source_node": 0, "Jmax": 3.5, "checks": {"ok": True},
                "source_commit": "a" * 40,
                "volume_backend": "tetrahedral-volume-v1",
                "volume_backend_sha256": "c" * 64,
            }
            meta_path = root / "v2.json"
            meta_path.write_text(json.dumps(metadata))
            spec = {
                "npz": str(npz), "sha256": hashlib.sha256(npz.read_bytes()).hexdigest(),
                "metadata_json": str(meta_path),
                "metadata_sha256": hashlib.sha256(meta_path.read_bytes()).hexdigest(),
            }
            with self.assertRaisesRegex(ValueError, "V2_PRIMITIVE_BAD_INTEGER_INDEX"):
                read_v2(spec, expected_source=0, expected_mode="forward",
                        expected_index=0, expected_source_commit="a" * 40,
                        expected_backend_sha="c" * 64)

    def test_49_attestation_stage_number_must_be_integer_not_boolean(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            item = {}
            for stage in (1, 2, 3):
                prefix = root / f"stage{stage}"
                digest = save_charged(prefix.with_suffix(".npz"), self.st)
                meta = {"stage": (True if stage == 1 else stage), "passed": True,
                        "npz_sha256": digest, "output_state_sha256": digest}
                prefix.with_suffix(".json").write_text(json.dumps(meta))
                item[f"stage{stage}"] = str(prefix)
            with self.assertRaisesRegex(ValueError, "STAGE_NUMBER_MISMATCH_STAGE1"):
                attest_manifest({"paths": [item]})


    def test_50_path_auxiliary_indices_reject_boolean_values(self):
        meta = {"mode": "forward", "order": "p", "aux": [False, 0, 1],
                "source_node": 0, "omit": 0, "cycle": 0, "Jmax2": 7}
        with self.assertRaisesRegex(ValueError, "INVALID_PATH_INDEX"):
            validate_path_metadata(meta)

    def test_51_path_provenance_indices_reject_boolean_values(self):
        meta = {"mode": "forward", "order": "p", "aux": [0, 0, 1],
                "source_node": True, "omit": 0, "cycle": 0, "Jmax2": 7}
        with self.assertRaisesRegex(ValueError, "BAD_STAGE_INTEGER_METADATA_SOURCE_NODE"):
            validate_path_metadata(meta)

    def test_52_path_metadata_requires_exact_integer_auxiliary_values(self):
        meta = {"mode": "forward", "order": "p", "aux": [0, 1.0, 1],
                "source_node": 0, "omit": 0, "cycle": 0, "Jmax2": 7}
        with self.assertRaisesRegex(ValueError, "INVALID_PATH_INDEX"):
            validate_path_metadata(meta)

    def test_53_valid_path_metadata_is_canonicalized(self):
        meta = {"mode": "adjoint", "order": "q", "aux": [1, 0, 1],
                "source_node": 0, "omit": 2, "cycle": 1, "Jmax2": 7}
        self.assertEqual(validate_path_metadata(meta), ("adjoint", "q", (1, 0, 1)))


    def test_54_manifest_source_node_rejects_boolean(self):
        manifest = {"source_node": True, "omit": 0, "cycle": 0, "Jmax2": 7}
        with self.assertRaisesRegex(ValueError, "BAD_MANIFEST_INTEGER_SOURCE_NODE"):
            validate_manifest_indices(manifest)

    def test_55_manifest_cycle_rejects_boolean(self):
        manifest = {"source_node": 0, "omit": 0, "cycle": False, "Jmax2": 7}
        with self.assertRaisesRegex(ValueError, "BAD_MANIFEST_INTEGER_CYCLE"):
            validate_manifest_indices(manifest)

    def test_56_manifest_indices_accept_canonical_integers(self):
        manifest = {"source_node": 0, "omit": 2, "cycle": 1, "Jmax2": 7}
        self.assertEqual(validate_manifest_indices(manifest), (0, 2, 1))


    def test_57_state_support_metadata_must_be_integer_not_boolean(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "stage"
            digest = save_charged(prefix.with_suffix(".npz"), self.st)
            diagnostics = {"CV_complete_basis_leakage": 0.0}
            guard = validate_physical_stage(
                diagnostics, self.st, 1, ["CV", 4, 1, 0])
            meta = dict(self.one, passed=True, status="REAL_ORIGINAL_STAGE_ONLY",
                        npz_sha256=digest, output_state_sha256=digest,
                        state_support=True, stage=1, diagnostics=diagnostics,
                        operator_signature=["CV", 4, 1, 0], guard=guard)
            meta_path = prefix.with_suffix(".json")
            meta_path.write_text(json.dumps(meta))
            attestation = {
                "metadata_sha256": hashlib.sha256(meta_path.read_bytes()).hexdigest(),
                "npz_sha256": digest,
            }
            with self.assertRaisesRegex(ValueError, "BAD_STATE_SUPPORT_METADATA"):
                read_stage(prefix, attestation)


    def test_58_predecessor_support_must_be_positive_nonboolean_integer(self):
        for invalid in (True, 1.0, 0, -1):
            with self.subTest(state_support=invalid):
                with self.assertRaisesRegex(ValueError, "BAD_PREDECESSOR_STATE_SUPPORT"):
                    predecessor_ok(dict(self.one, state_support=invalid),
                                   self.two, 2, ["CV", 4, 1, 0])


    def _v2_term_record(self, root, **overrides):
        path = Path(root) / "term_forward_0.json"
        npz = path.with_suffix(".npz")
        np.savez_compressed(npz, spins=np.array([[1] * 32]),
                            Ks=np.array([[0] * 16]), amp=np.array([1 + 0j]))
        metadata = {
            "passed": True, "mode": "forward", "index": 0, "source_node": 0,
            "PL_epsilon_coefficient": 1, "Jmax": 3.5, "gauss_support": 1,
            "checks": {
                "finite_covariant_norm": True, "finite_gauss_norm": True,
                "finite_covariant_amplitudes": True, "finite_gauss_amplitudes": True,
                "physical_basis_volume_leakage": True,
                "scalar_closure_fraction_or_exact_zero": True,
                "nonscalar_rejected_norm": True, "single_L_spin_wall": True,
                "PL_orientation_coefficient": True,
            },
            "npz_sha256": hashlib.sha256(npz.read_bytes()).hexdigest(),
        }
        metadata.update(overrides)
        return path, metadata

    def test_59_v2_term_passed_must_be_boolean_true(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._v2_term_record(directory, passed="yes")
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_NOT_PASSED_BOOLEAN_TRUE"):
                validate_term_record(path, meta)

    def test_60_v2_term_index_must_be_integer_not_boolean(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._v2_term_record(directory, index=True)
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_BAD_INTEGER_INDEX"):
                validate_term_record(path, meta)

    def test_61_v2_term_source_node_must_be_integer_not_boolean(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._v2_term_record(directory, source_node=False)
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_BAD_SOURCE_NODE"):
                validate_term_record(path, meta)

    def test_62_v2_orientation_coefficient_must_be_integer_not_boolean(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._v2_term_record(directory, PL_epsilon_coefficient=True)
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_BAD_ORIENTATION_COEFFICIENT"):
                validate_term_record(path, meta)

    def test_63_v2_term_npz_must_match_metadata_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._v2_term_record(directory, npz_sha256="0" * 64)
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_NPZ_DIGEST_MISMATCH"):
                validate_term_record(path, meta)

    def test_64_v2_term_record_accepts_exact_provenance_and_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._v2_term_record(directory)
            self.assertEqual(validate_term_record(path, meta), ("forward", 0, 0))


    def test_65_v2_metadata_must_bind_exact_npz_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            npz = root / "v2.npz"
            np.savez_compressed(npz, spins=np.array([[1] * 32]),
                                Ks=np.array([[0] * 16]), amp=np.array([1 + 0j]))
            digest = hashlib.sha256(npz.read_bytes()).hexdigest()
            metadata = {
                "passed": True, "mode": "forward", "index": 0,
                "operator_version": "tetrahedral-charged-volume-v2",
                "source_node": 0, "Jmax": 3.5, "checks": {"ok": True},
                "source_commit": "a" * 40,
                "volume_backend": "tetrahedral-volume-v1",
                "volume_backend_sha256": "c" * 64,
                "npz_sha256": "0" * 64,
            }
            meta_path = root / "v2.json"
            meta_path.write_text(json.dumps(metadata))
            spec = {
                "npz": str(npz), "sha256": digest,
                "metadata_json": str(meta_path),
                "metadata_sha256": hashlib.sha256(meta_path.read_bytes()).hexdigest(),
            }
            with self.assertRaisesRegex(ValueError, "V2_METADATA_NPZ_DIGEST_MISMATCH"):
                read_v2(spec, expected_source=0, expected_mode="forward",
                        expected_index=0, expected_source_commit="a" * 40,
                        expected_backend_sha="c" * 64)


    def test_66_v2_term_npz_shape_and_dtype_are_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._v2_term_record(directory)
            npz = path.with_suffix(".npz")
            np.savez_compressed(npz, spins=np.array([[1] * 31]),
                                Ks=np.array([[0] * 16]), amp=np.array([1 + 0j]))
            meta["npz_sha256"] = hashlib.sha256(npz.read_bytes()).hexdigest()
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_NPZ_BAD_SHAPE_DTYPE_OR_AMPLITUDE"):
                validate_term_record(path, meta)

    def test_67_v2_term_npz_duplicate_keys_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._v2_term_record(directory)
            npz = path.with_suffix(".npz")
            np.savez_compressed(npz, spins=np.array([[1] * 32, [1] * 32]),
                                Ks=np.array([[0] * 16, [0] * 16]),
                                amp=np.array([1 + 0j, 2 + 0j]))
            meta["gauss_support"] = 2
            meta["npz_sha256"] = hashlib.sha256(npz.read_bytes()).hexdigest()
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_NPZ_DUPLICATE_KEY"):
                validate_term_record(path, meta)

    def test_68_v2_term_support_must_be_integer_not_boolean(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._v2_term_record(directory, gauss_support=True)
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_NPZ_SUPPORT_MISMATCH"):
                validate_term_record(path, meta)


    def _attested_v2_read_spec(self, root, **overrides):
        root = Path(root)
        npz = root / "v2.npz"
        np.savez_compressed(npz, spins=np.array([[1] * 32]),
                            Ks=np.array([[0] * 16]), amp=np.array([1 + 0j]))
        digest = hashlib.sha256(npz.read_bytes()).hexdigest()
        metadata = {
            "passed": True, "mode": "forward", "index": 0,
            "operator_version": "tetrahedral-charged-volume-v2",
            "source_node": 0, "Jmax": 3.5, "gauss_support": 1,
            "checks": {
                "finite_covariant_norm": True, "finite_gauss_norm": True,
                "finite_covariant_amplitudes": True, "finite_gauss_amplitudes": True,
                "physical_basis_volume_leakage": True,
                "scalar_closure_fraction_or_exact_zero": True,
                "nonscalar_rejected_norm": True, "single_L_spin_wall": True,
                "PL_orientation_coefficient": True,
            }, "source_commit": "a" * 40,
            "volume_backend": "tetrahedral-volume-v1",
            "volume_backend_sha256": "c" * 64, "npz_sha256": digest,
        }
        metadata.update(overrides)
        meta_path = root / "v2.json"
        meta_path.write_text(json.dumps(metadata))
        spec = {
            "npz": str(npz), "sha256": digest,
            "metadata_json": str(meta_path),
            "metadata_sha256": hashlib.sha256(meta_path.read_bytes()).hexdigest(),
        }
        return spec

    def test_69_v2_cutoff_must_be_numeric_not_string(self):
        with tempfile.TemporaryDirectory() as directory:
            spec = self._attested_v2_read_spec(directory, Jmax="3.5")
            with self.assertRaisesRegex(ValueError, "V2_PRIMITIVE_WRONG_CUTOFF"):
                read_v2(spec, expected_source=0, expected_mode="forward",
                        expected_index=0, expected_source_commit="a" * 40,
                        expected_backend_sha="c" * 64)

    def test_70_v2_attested_support_must_be_exact_integer(self):
        with tempfile.TemporaryDirectory() as directory:
            spec = self._attested_v2_read_spec(directory, gauss_support=True)
            with self.assertRaisesRegex(ValueError, "V2_PRIMITIVE_SUPPORT_MISMATCH"):
                read_v2(spec, expected_source=0, expected_mode="forward",
                        expected_index=0, expected_source_commit="a" * 40,
                        expected_backend_sha="c" * 64)


    def test_71_v2_term_requires_complete_acceptance_check_inventory(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._v2_term_record(
                directory, checks={"finite_covariant_norm": True})
            with self.assertRaisesRegex(
                    RuntimeError, "V2_TERM_CHECKS_INCOMPLETE_OR_NOT_BOOLEAN_TRUE"):
                validate_term_record(path, meta)

    def test_72_v2_attested_reader_requires_complete_acceptance_check_inventory(self):
        with tempfile.TemporaryDirectory() as directory:
            spec = self._attested_v2_read_spec(
                directory, checks={"finite_covariant_norm": True})
            with self.assertRaisesRegex(
                    ValueError, "V2_PRIMITIVE_FAILED_ACCEPTANCE_OR_INCOMPLETE_CHECKS"):
                read_v2(spec, expected_source=0, expected_mode="forward",
                        expected_index=0, expected_source_commit="a" * 40,
                        expected_backend_sha="c" * 64)


    def test_73_v2_term_requires_all_physical_diagnostics(self):
        with self.assertRaisesRegex(ValueError, "MISSING_TERM_DIAGNOSTIC_CV_COMPLETE_BASIS_LEAKAGE"):
            validate_term_diagnostics({}, 1.0, 0.0)

    def test_74_v2_term_rejects_negative_or_nonfinite_leakage(self):
        good = {
            "CV_complete_basis_leakage": 0.0,
            "CK_outer_complete_basis_leakage": 0.0,
            "CK_internal_volume_sector_leakage": 0.0,
            "CK_complete_charge_basis_leakage": 0.0,
        }
        for key, bad in (("CV_complete_basis_leakage", -1e-9),
                         ("CK_complete_charge_basis_leakage", float("nan"))):
            with self.subTest(key=key, value=bad):
                diag = dict(good, **{key: bad})
                with self.assertRaisesRegex(ValueError, "INVALID_TERM_DIAGNOSTIC_"):
                    validate_term_diagnostics(diag, 1.0, 0.0)

    def test_75_v2_term_rejects_negative_or_nonfinite_projection_norms(self):
        good = {
            "CV_complete_basis_leakage": 0.0,
            "CK_outer_complete_basis_leakage": 0.0,
            "CK_internal_volume_sector_leakage": 0.0,
            "CK_complete_charge_basis_leakage": 0.0,
        }
        for accepted, rejected in ((-1.0, 0.0), (1.0, -1.0),
                                   (float("inf"), 0.0), (1.0, float("nan"))):
            with self.subTest(accepted=accepted, rejected=rejected):
                with self.assertRaisesRegex(ValueError, "INVALID_PROJECTION_NORM_"):
                    validate_term_diagnostics(good, accepted, rejected)


    def test_76_v2_diagnostic_rejects_boolean_numeric_value(self):
        diag = {
            "CV_complete_basis_leakage": True,
            "CK_outer_complete_basis_leakage": 0.0,
            "CK_internal_volume_sector_leakage": 0.0,
            "CK_complete_charge_basis_leakage": 0.0,
        }
        with self.assertRaisesRegex(ValueError, "INVALID_TERM_DIAGNOSTIC_CV_COMPLETE_BASIS_LEAKAGE"):
            validate_term_diagnostics(diag, 1.0, 0.0)

    def test_77_v2_diagnostic_rejects_numeric_string(self):
        diag = {
            "CV_complete_basis_leakage": "0.0",
            "CK_outer_complete_basis_leakage": 0.0,
            "CK_internal_volume_sector_leakage": 0.0,
            "CK_complete_charge_basis_leakage": 0.0,
        }
        with self.assertRaisesRegex(ValueError, "INVALID_TERM_DIAGNOSTIC_CV_COMPLETE_BASIS_LEAKAGE"):
            validate_term_diagnostics(diag, 1.0, 0.0)

    def test_78_v2_projection_norm_rejects_boolean(self):
        diag = {
            "CV_complete_basis_leakage": 0.0,
            "CK_outer_complete_basis_leakage": 0.0,
            "CK_internal_volume_sector_leakage": 0.0,
            "CK_complete_charge_basis_leakage": 0.0,
        }
        with self.assertRaisesRegex(ValueError, "INVALID_PROJECTION_NORM_ACCEPTED2"):
            validate_term_diagnostics(diag, True, 0.0)


    def test_79_exact_complex_aggregation_is_randomized_order_invariant(self):
        rng = Random(20261108)
        states = []
        for _ in range(80):
            state = {}
            for key_index in range(5):
                scale = 1e12 if rng.random() < 0.5 else 1.0
                sign = -1.0 if rng.random() < 0.5 else 1.0
                real = sign * scale + rng.choice((-1.0, 0.0, 1.0))
                imag = rng.uniform(-2.0, 2.0)
                state[((key_index,), (0,))] = complex(real, imag)
            states.append(state)
        expected = sum_states_exact(states)
        for _ in range(32):
            permuted = list(states)
            rng.shuffle(permuted)
            permuted = [dict(reversed(list(state.items()))) for state in permuted]
            self.assertEqual(sum_states_exact(permuted), expected)



if __name__ == "__main__":
    unittest.main()
