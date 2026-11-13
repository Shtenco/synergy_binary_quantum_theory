import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_v10_w5_convergence_grid as grid


class TestW5ConvergenceGrid(unittest.TestCase):
    @staticmethod
    def complete_result(shadow=False):
        row = {
            "basis_keys_repr": ["basis-0", "basis-1"],
            "forward_B_columns_sparse": [
                [[0, 1.0, 0.0]],
                [[1, 0.5, 0.25]],
            ],
            "direct_Bdag_columns_sparse": [
                [[0, 1.0, 0.0]],
                [[1, 0.5, -0.25]],
            ],
            "diagnostics": {
                "spin_wall_omitted_branch_count": 3,
                "max_intermediate_spin_reached": 4.0,
            },
            "coefficient_loss_ledger": {
                "summary": {"closure_caps.forward": {"cap_reason_counts": {"max_depth": 1}}}
            },
            "closure_size": 2,
            "basis_id_sha256": "fixture-basis",
            "spectral_sensitivity": {"zero_factor": 1000.0},
        }
        if shadow:
            row.update({
                "status": "INCOMPLETE_NO_PRUNING_COVERAGE",
                "pruning_mode": "shadow_no_prune",
                "algebraic_status": "PASS_ADJOINT_CONSISTENCY",
            })
        else:
            row["status"] = "PASS_ADJOINT_CONSISTENCY"
        return row

    def test_matrix_digest_is_independent_of_dictionary_insertion_order(self):
        a = {"input-b": {"out-2": 2 + 1j}, "input-a": {"out-1": 1 - 3j}}
        b = {"input-a": {"out-1": 1 - 3j}, "input-b": {"out-2": 2 + 1j}}
        self.assertEqual(grid.matrix_digest(a), grid.matrix_digest(b))

    def test_grid_compares_full_complex_support_across_tol_and_spin_walls(self):
        args = SimpleNamespace(
            Jmax2_grid=[10, 12],
            amp_tol_grid=[1e-10, 1e-12],
            max_states=500,
            max_depth=8,
            unbounded_closure=False,
            nodes=[0],
            timeout_seconds=1.0,
            cache_mode="on",
            spectral_zero_factor=1000.0,
        )

        def fake_run_one(script, output, amp_tol, run_args, pruning_mode="production"):
            return self.complete_result(shadow=pruning_mode == "shadow_no_prune")

        with patch.object(grid, "run_one", side_effect=fake_run_one):
            result = grid.run_grid(args)

        self.assertEqual(result["status"], "COMPLETED_FINITE_GRID_SHADOW_REFERENCE_INCOMPLETE")
        self.assertEqual(result["expected_cells"], 4)
        self.assertEqual(result["completed_production_and_shadow_cells"], 4)
        self.assertEqual(len(result["amp_tol_convergence"]), 2)
        self.assertEqual(len(result["Jmax2_convergence"]), 2)
        self.assertTrue(result["all_production_matrices_complete"])
        self.assertTrue(result["all_shadow_matrices_complete"])
        self.assertFalse(result["zero_loss_proved"])
        for row in result["grid"]:
            self.assertEqual(row["production_vs_shadow_forward"]["support_symmetric_difference_total"], 0)
            self.assertEqual(row["production_vs_shadow_forward"]["max_coefficient_absolute_residual"], 0.0)
            self.assertEqual(row["production_vs_shadow_direct_adjoint"]["support_symmetric_difference_total"], 0)
            self.assertEqual(row["production_vs_shadow_direct_adjoint"]["max_coefficient_absolute_residual"], 0.0)
            self.assertEqual(row["production_diagnostics"]["spin_wall_omitted_branch_count"], 3)
            self.assertEqual(
                row["production_loss_summary"]["closure_caps.forward"]["cap_reason_counts"]["max_depth"],
                1,
            )
        self.assertEqual(
            result["amp_tol_convergence"][0]["direct_adjoint_comparison"]["support_symmetric_difference_total"],
            0,
        )
        self.assertEqual(
            result["Jmax2_convergence"][0]["direct_adjoint_comparison"]["max_coefficient_absolute_residual"],
            0.0,
        )


if __name__ == "__main__":
    unittest.main()
