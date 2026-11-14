import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch, Mock
import subprocess

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_v10_w5_amp_tol_shadow_compare as shadow_compare
from bqg_v10_w5_amp_tol_shadow_compare import (
    closure_matrix_available, compare_vectors, compare_cache_runs, unpack_columns, run_one,
)


class TestV10ShadowVectorComparison(unittest.TestCase):
    def test_compares_complex_sparse_vectors_and_support(self):
        reference = {
            "input-A": {"out-1": 1 + 2j, "out-2": 0.5 - 0.25j},
        }
        candidate = {
            "input-A": {"out-1": 1 + 2j, "out-3": 0.5 - 0.25j},
        }
        result = compare_vectors(reference, candidate)
        self.assertEqual(result["support_symmetric_difference_total"], 2)
        self.assertAlmostEqual(result["max_coefficient_absolute_residual"], abs(0.5 - 0.25j))
        self.assertGreater(result["absolute_l2_residual"], 0.0)

    def test_full_shadow_matrix_is_comparable_but_not_certified(self):
        shadow = {
            "status": "INCOMPLETE_NO_PRUNING_COVERAGE",
            "pruning_mode": "shadow_no_prune",
            "algebraic_status": "PASS_ADJOINT_CONSISTENCY",
            "forward_B_columns_sparse": [],
        }
        self.assertTrue(closure_matrix_available(shadow, shadow=True))
        self.assertFalse(closure_matrix_available(shadow, shadow=False))

    def test_shadow_matrix_is_not_comparable_when_closure_is_capped(self):
        shadow = {
            "status": "INCOMPLETE_CLOSURE_LIMIT_REACHED",
            "pruning_mode": "shadow_no_prune",
            "algebraic_status": "PASS_ADJOINT_CONSISTENCY",
        }
        self.assertFalse(closure_matrix_available(shadow, shadow=True))

    def test_decodes_full_complex_columns(self):
        result = {
            "basis_keys_repr": ["input", "output"],
            "forward_B_columns_sparse": [
                [[1, 0.25, -0.5]],
                [[0, 1.0, 2.0]],
            ],
        }
        decoded = unpack_columns(result)
        self.assertEqual(decoded["input"]["output"], 0.25 - 0.5j)
        self.assertEqual(decoded["output"]["input"], 1.0 + 2.0j)


if __name__ == "__main__":
    unittest.main()


class TestV10CacheEquivalence(unittest.TestCase):
    @staticmethod
    def result_fixture():
        return {
            "status": "PASS_ADJOINT_CONSISTENCY",
            "basis_keys_repr": ["a", "b"],
            "forward_B_columns_sparse": [[[0, 1.0, 0.0]], [[1, 0.5, 0.25]]],
            "direct_Bdag_columns_sparse": [[[0, 1.0, 0.0]], [[1, 0.5, -0.25]]],
            "lower_level_pruning_ledger": {
                "event_count": 2,
                "events": [{"site_id": "x", "propagated_abs": 1e-13}],
                "coverage": {"complete": False, "static_threshold_inventory_complete": True, "fixed_operator_mask_error_bound_complete": True},
                "cache_policy": "cache-on",
            },
            "coefficient_loss_ledger": {"summary": {"amp_tol.forward": {"count": 1}}},
            "coherent_node_sum_ledger": [{"input_key": "a", "forward": {"post_sum_coherent_l1": 1.0}}],
            "diagnostics": {
                "spin_wall_omitted_branch_count": 0,
                "max_intermediate_spin_reached": 3.0,
                "closure_cap_counts": {"max_states": 0, "max_depth": 0},
            },
            "spectral_sensitivity": {"zero_factor": 1000.0, "zeroed_eigenvalue_count": 1},
        }

    def test_full_cache_matrix_and_every_ledger_compare_exactly(self):
        cache_on = self.result_fixture()
        cache_off = self.result_fixture()
        cache_off["lower_level_pruning_ledger"]["cache_policy"] = "cache-off"

        result = compare_cache_runs(cache_on, cache_off)
        self.assertEqual(result["status"], "PASS_CACHE_EQUIVALENCE")
        self.assertTrue(result["equivalence_passed"])
        self.assertTrue(result["static_coverage_complete"])
        self.assertTrue(result["full_ledger_exactly_equal"])
        self.assertTrue(result["observed_ledger_exactly_equal"])
        self.assertTrue(result["matrix_exactly_equal"])
        self.assertTrue(result["full_ledger_exactly_equal"])
        self.assertTrue(all(result["ledger_comparison"].values()))
        self.assertEqual(result["cache_on_event_count"], 2)
        self.assertEqual(result["cache_off_event_count"], 2)

    def test_cache_comparison_detects_spin_wall_and_spectral_differences(self):
        cache_on = self.result_fixture()
        cache_off = self.result_fixture()
        cache_off["diagnostics"]["spin_wall_omitted_branch_count"] = 1
        cache_off["spectral_sensitivity"]["zeroed_eigenvalue_count"] = 2

        result = compare_cache_runs(cache_on, cache_off)
        self.assertTrue(result["matrix_exactly_equal"])
        self.assertFalse(result["ledger_comparison"]["spin_wall_and_closure_diagnostics"])
        self.assertFalse(result["ledger_comparison"]["spectral_sensitivity_ledger"])
        self.assertFalse(result["full_ledger_exactly_equal"])
        self.assertEqual(result["status"], "FAIL_CACHE_EQUIVALENCE")
        self.assertFalse(result["equivalence_passed"])

    def test_cache_equivalence_uses_static_coverage_not_circular_runtime_flag(self):
        cache_on = self.result_fixture()
        cache_off = self.result_fixture()
        cache_off["lower_level_pruning_ledger"]["cache_policy"] = "cache-off"
        cache_on["lower_level_pruning_ledger"]["coverage"]["static_threshold_inventory_complete"] = False
        cache_off["lower_level_pruning_ledger"]["coverage"]["static_threshold_inventory_complete"] = False

        result = compare_cache_runs(cache_on, cache_off)
        self.assertTrue(result["matrix_exactly_equal"])
        self.assertTrue(result["observed_ledger_exactly_equal"])
        self.assertFalse(result["static_coverage_complete"])
        self.assertFalse(result["full_ledger_exactly_equal"])
        self.assertEqual(result["status"], "INCOMPLETE_CACHE_EQUIVALENCE_STATIC_COVERAGE")
        self.assertFalse(result["equivalence_passed"])

    def test_cache_comparison_fails_closed_when_required_ledger_is_missing(self):
        cache_on = self.result_fixture()
        cache_off = self.result_fixture()
        cache_off.pop("coherent_node_sum_ledger")

        result = compare_cache_runs(cache_on, cache_off)
        self.assertFalse(result["ledger_payloads_present"])
        self.assertFalse(result["observed_ledger_exactly_equal"])
        self.assertEqual(result["status"], "INCOMPLETE_CACHE_EQUIVALENCE_LEDGER_MISSING")
        self.assertFalse(result["equivalence_passed"])

    def test_missing_lower_level_ledger_fails_closed_without_exception(self):
        cache_on = self.result_fixture()
        cache_off = self.result_fixture()
        cache_off.pop("lower_level_pruning_ledger")

        result = compare_cache_runs(cache_on, cache_off)
        self.assertFalse(result["ledger_payloads_present"])
        self.assertFalse(result["static_coverage_complete"])
        self.assertEqual(result["status"], "INCOMPLETE_CACHE_EQUIVALENCE_LEDGER_MISSING")

    def test_cache_comparison_status_fails_on_matrix_mismatch(self):
        cache_on = self.result_fixture()
        cache_off = self.result_fixture()
        cache_off["forward_B_columns_sparse"][0][0][1] += 1e-14

        result = compare_cache_runs(cache_on, cache_off)
        self.assertFalse(result["matrix_exactly_equal"])
        self.assertTrue(result["full_ledger_exactly_equal"])
        self.assertEqual(result["status"], "FAIL_CACHE_EQUIVALENCE")
        self.assertFalse(result["equivalence_passed"])


class TestV10ShadowTimeoutHandling(unittest.TestCase):
    def test_unbounded_closure_flag_is_forwarded_to_both_processes(self):
        import tempfile

        args = SimpleNamespace(
            Jmax2=10, max_states=500, max_depth=8, nodes=[0],
            timeout_seconds=1.0, unbounded_closure=True,
            spectral_zero_factor=2000.0, cache_mode="off",
        )
        completed = Mock(returncode=0, stdout="", stderr="")
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "result.json"

            def fake_run(command, **kwargs):
                self.assertIn("--unbounded-closure", command)
                self.assertIn("--spectral-zero-factor", command)
                self.assertIn("2000.0", command)
                self.assertIn("--cache-mode", command)
                self.assertIn("off", command)
                output.write_text('{"status": "INCOMPLETE_CLOSURE_LIMIT_REACHED"}', encoding="utf-8")
                return completed

            with patch(
                "bqg_v10_w5_amp_tol_shadow_compare.subprocess.run",
                side_effect=fake_run,
            ):
                result = run_one(Path(temp) / "v10.py", output, 1e-12, args)
        self.assertIn("--unbounded-closure", result["_process"]["command"])

    def test_timeout_is_incomplete_and_partial_json_is_discarded(self):
        import tempfile

        args = SimpleNamespace(
            Jmax2=10, max_states=5, max_depth=0, nodes=[0],
            timeout_seconds=0.01,
        )
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "partial.json"
            output.write_text('{"partial": true}', encoding="utf-8")
            with patch(
                "bqg_v10_w5_amp_tol_shadow_compare.subprocess.run",
                side_effect=subprocess.TimeoutExpired(
                    cmd=["python"], timeout=0.01,
                    output=b"partial stdout", stderr=b"partial stderr",
                ),
            ):
                result = run_one(Path(temp) / "v10.py", output, 1e-12, args)

            self.assertEqual(result["status"], "INCOMPLETE_TIMEOUT")
            self.assertTrue(result["_process"]["partial_output_discarded"])
            self.assertFalse(output.exists())

    def test_secondary_audits_are_skipped_when_baseline_matrix_times_out(self):
        import contextlib
        import io
        import tempfile

        calls = []
        def incomplete_run_one(script, output, amp_tol, args, pruning_mode="production"):
            calls.append((amp_tol, pruning_mode, args.cache_mode))
            return {"status": "INCOMPLETE_TIMEOUT", "_process": {"timeout_seconds": 0.01}}

        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "comparison.json"
            argv = [
                "bqg_v10_w5_amp_tol_shadow_compare.py",
                "--compare-cache-modes",
                "--spectral-sensitivity-factors", "2000",
                "--timeout-seconds", "0.01",
                "--nodes", "0",
                "--output", str(output),
            ]
            with patch("sys.argv", argv), patch.object(
                shadow_compare, "run_one", side_effect=incomplete_run_one
            ), contextlib.redirect_stdout(io.StringIO()):
                exit_code = shadow_compare.main()

            result = __import__("json").loads(output.read_text(encoding="utf-8"))

        self.assertEqual(exit_code, 2)
        self.assertEqual(len(calls), 2)  # production and shadow only
        self.assertEqual(result["cache_equivalence"]["status"], "SKIPPED_BASELINE_MATRIX_INCOMPLETE")
        self.assertEqual(result["spectral_sensitivity"]["status"], "SKIPPED_BASELINE_MATRIX_INCOMPLETE")
        self.assertFalse(result["comparison_performed"])

    def test_completed_process_without_json_is_fail_closed(self):
        import tempfile

        args = SimpleNamespace(
            Jmax2=10, max_states=5, max_depth=0, nodes=[0],
            timeout_seconds=1.0,
        )
        completed = Mock(returncode=1, stdout="", stderr="crash")
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "missing.json"
            with patch(
                "bqg_v10_w5_amp_tol_shadow_compare.subprocess.run",
                return_value=completed,
            ):
                result = run_one(Path(temp) / "v10.py", output, 1e-12, args)
        self.assertEqual(result["status"], "INCOMPLETE_NO_RESULT_JSON")
        self.assertEqual(result["_process"]["returncode"], 1)


if __name__ == "__main__":
    unittest.main()
