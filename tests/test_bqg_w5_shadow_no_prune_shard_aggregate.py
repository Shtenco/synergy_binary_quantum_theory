import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from bqg_w5_shadow_no_prune_shard_aggregate import aggregate


class TestW5ShadowShardAggregate(unittest.TestCase):
    def report(self, jmax2, node, mode, *, matrix_value=1.0, unbounded=True):
        return {
            "status": "PASS_ADJOINT_CONSISTENCY" if mode == "production" else "INCOMPLETE_NO_PRUNING_COVERAGE",
            "algebraic_status": "PASS_ADJOINT_CONSISTENCY",
            "lower_level_pruning_coverage_complete": False,
            "source_nodes": [node],
            "pruning_mode": mode,
            "cache_mode": "on",
            "Jmax2": jmax2,
            "seed_count": 5,
            "seed_Ks": [0, 2, 4, 6, 8],
            "max_states": None if unbounded else 500,
            "max_depth": None if unbounded else 8,
            "unbounded_closure": unbounded,
            "omitted_new_support_keys": 0,
            "closure_size": 2,
            "basis_id_sha256": hashlib.sha256(b"basis-0\nbasis-1").hexdigest(),
            "basis_keys_repr": ["basis-0", "basis-1"],
            "forward_B_columns_sparse": [[[0, matrix_value, 0.0]], [[1, 2.0, 0.0]]],
            "direct_Bdag_columns_sparse": [[[0, 3.0, 0.0]], [[1, 4.0, 0.0]]],
            "coefficient_loss_ledger": {"amp_tol": {"count": 1}},
            "coherent_node_sum_ledger": [{"input_key": "k", "forward": {"norm": 1.0}}],
            "lower_level_pruning_ledger": {"event_count": 1, "events": [], "by_site": {}, "by_site_and_kind": {}, "coverage": {}},
            "diagnostics": {"direct_adjoint_relative_residual": 0.0},
            "physical_status": {"physical_hermitian_ordering_selected": False, "full_physical_HL": False},
        }

    def write_pair(self, directory, jmax2, node, *, difference=False, unbounded=True):
        for mode in ("production", "shadow_no_prune"):
            value = 1.1 if difference and mode == "shadow_no_prune" else 1.0
            report = self.report(jmax2, node, mode, matrix_value=value, unbounded=unbounded)
            path = directory / f"W5_SHADOW_MATRIX_J{jmax2}_NODE{node}_MODE{mode}_RUN1.json"
            path.write_text(json.dumps(report), encoding="utf-8")

    def fill(self, directory, *, difference_pair=None, unbounded=True):
        for node in range(5):
            self.write_pair(directory, 10, node, difference=(node == difference_pair), unbounded=unbounded)

    def test_all_five_node_pairs_complete_but_claim_remains_open(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            result = aggregate(directory, (10,))
        self.assertTrue(result["finite_matrix_comparison_complete"], result["errors"])
        self.assertFalse(result["passed"])
        self.assertEqual(result["status"], "INCOMPLETE_NO_PRUNING_COVERAGE_ALL_SHARDS_COMPLETED")
        self.assertEqual(result["completed_pair_count"], 5)
        self.assertEqual(result["valid_mode_shard_count"], 10)

    def test_matrix_difference_is_reported_not_mislabelled_as_infrastructure_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory, difference_pair=0)
            result = aggregate(directory, (10,))
        self.assertTrue(result["finite_matrix_comparison_complete"])
        pair = next(x for x in result["pairs"] if x["node"] == 0)
        self.assertGreater(pair["matrix_comparisons"]["forward_B_columns_sparse"]["relative_l2_residual_on_common_domain"], 0.0)

    def test_missing_mode_shard_is_incomplete(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            (directory / "W5_SHADOW_MATRIX_J10_NODE4_MODEshadow_no_prune_RUN1.json").unlink()
            result = aggregate(directory, (10,))
        self.assertFalse(result["finite_matrix_comparison_complete"])
        self.assertEqual(len(result["missing_shards"]), 1)


    def test_incorrect_basis_digest_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            path = directory / "W5_SHADOW_MATRIX_J10_NODE0_MODEproduction_RUN1.json"
            report = json.loads(path.read_text())
            report["basis_id_sha256"] = "a" * 64
            path.write_text(json.dumps(report), encoding="utf-8")
            result = aggregate(directory, (10,))
        self.assertFalse(result["finite_matrix_comparison_complete"])
        self.assertTrue(any("failed checks: basis_hash" in error for error in result["errors"]))

    def test_bounded_closure_cannot_be_aggregated_as_unbounded(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory, unbounded=False)
            result = aggregate(directory, (10,))
        self.assertFalse(result["finite_matrix_comparison_complete"])
        self.assertTrue(any("failed checks: unbounded_closure" in error for error in result["errors"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
