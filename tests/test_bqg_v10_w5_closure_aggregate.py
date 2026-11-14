import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from bqg_v10_w5_closure_aggregate import aggregate_records
from bqg_v10_w5_closure_summary import matrix_digest


def record(j, node, k, digest="same"):
    return {
        "Jmax2": j, "source_nodes": [node], "seed_Ks": [k],
        "status": "PASS_ADJOINT_CONSISTENCY",
        "algebraic_status": "PASS_ADJOINT_CONSISTENCY",
        "matrices_present": True, "basis_id_sha256": "basis",
        "closure_size": 12,
        "forward_frobenius_norm": 1.0,
        "direct_adjoint_frobenius_norm": 1.0,
        "direct_adjoint_relative_residual": 1e-12,
        "max_component_leakage": 1e-14,
        "scalar_projection_rejected_terms": 0,
        "max_spin_reached": 3.5,
        "max_intermediate_spin_reached": 3.5,
        "spin_wall_omitted_branch_count": 0,
        "sine_engine_gate_passed": True,
        "forward_matrix_digest_10sig": digest,
        "direct_adjoint_matrix_digest_10sig": digest,
    }


class TestW5ClosureAggregate(unittest.TestCase):
    def test_matrix_digest_is_stable_to_sub_precision_noise(self):
        a = [[[0, 1.234567890123, -0.000000000123456789]]]
        b = [[[0, 1.234567890124, -0.00000000012345678901]]]
        self.assertEqual(matrix_digest(a), matrix_digest(b))

    def test_matrix_digest_changes_for_material_coefficient_change(self):
        a = [[[0, 1.234567890123, 0.0]]]
        b = [[[0, 1.234567891123, 0.0]]]
        self.assertNotEqual(matrix_digest(a), matrix_digest(b))

    def test_complete_matrix_and_cutoff_stability_pass(self):
        rows = [record(j, 0, 0) for j in (10, 12, 14)]
        result = aggregate_records(rows, nodes=(0,), seed_ks=(0,), jmax2_values=(10, 12, 14))
        self.assertTrue(result["passed"])
        self.assertEqual(result["expected_cell_count"], 3)

    def test_missing_cutoff_cell_fails_closed(self):
        rows = [record(j, 0, 0) for j in (10, 14)]
        result = aggregate_records(rows, nodes=(0,), seed_ks=(0,), jmax2_values=(10, 12, 14))
        self.assertFalse(result["passed"])
        self.assertEqual(result["missing_cells"], [[12, 0, 0]])

    def test_matrix_change_across_cutoff_fails_convergence(self):
        rows = [record(10, 0, 0), record(12, 0, 0, "changed"), record(14, 0, 0)]
        result = aggregate_records(rows, nodes=(0,), seed_ks=(0,), jmax2_values=(10, 12, 14))
        self.assertFalse(result["passed"])
        self.assertTrue(result["cutoff_convergence_failures"])

    def test_nonpass_algebraic_status_fails(self):
        rows = [record(j, 0, 0) for j in (10, 12, 14)]
        rows[1]["algebraic_status"] = "FAIL_ADJOINT_OR_LEAKAGE_GATE"
        result = aggregate_records(rows, nodes=(0,), seed_ks=(0,), jmax2_values=(10, 12, 14))
        self.assertFalse(result["passed"])
        self.assertTrue(result["cell_failures"])


if __name__ == "__main__":
    unittest.main()
