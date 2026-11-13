import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch, Mock
import subprocess

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from bqg_v10_w5_amp_tol_shadow_compare import (
    closure_matrix_available, compare_vectors, unpack_columns, run_one,
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


class TestV10ShadowTimeoutHandling(unittest.TestCase):
    def test_unbounded_closure_flag_is_forwarded_to_both_processes(self):
        import tempfile

        args = SimpleNamespace(
            Jmax2=10, max_states=500, max_depth=8, nodes=[0],
            timeout_seconds=1.0, unbounded_closure=True,
        )
        completed = Mock(returncode=0, stdout="", stderr="")
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "result.json"

            def fake_run(command, **kwargs):
                self.assertIn("--unbounded-closure", command)
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
