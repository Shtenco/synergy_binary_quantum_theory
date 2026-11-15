import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from bqg_w5_shadow_no_prune_shard_aggregate import aggregate


def report(jmax2, node, mode, *, amplitude=1.0):
    status = "PASS_ADJOINT_CONSISTENCY" if mode == "production" else "INCOMPLETE_NO_PRUNING_COVERAGE"
    return {
        "status": status,
        "algebraic_status": "PASS_ADJOINT_CONSISTENCY",
        "source_nodes": [node],
        "seed_Ks": [0, 2, 4, 6, 8],
        "Jmax2": jmax2,
        "pruning_mode": mode,
        "cache_mode": "on",
        "max_states": None,
        "max_depth": None,
        "basis_keys_repr": ["basis-0"],
        "forward_B_columns_sparse": [[[0, amplitude, 0.0]]],
        "direct_Bdag_columns_sparse": [[[0, amplitude, 0.0]]],
        "coefficient_loss_ledger": {},
        "lower_level_pruning_ledger": {"event_count": 0},
        "coherent_node_sum_ledger": [],
    }


class TestW5ShadowPairedShardAggregate(unittest.TestCase):
    def write(self, directory, data, suffix="RUN1"):
        j, n, mode = data["Jmax2"], data["source_nodes"][0], data["pruning_mode"]
        path = directory / f"W5_SHADOW_MATRIX_J{j}_NODE{n}_MODE{mode}_{suffix}.json"
        path.write_text(json.dumps(data), encoding="utf-8")

    def fill(self, directory):
        for j in (10, 12, 14):
            for n in range(5):
                for mode in ("production", "shadow_no_prune"):
                    self.write(directory, report(j, n, mode))

    def test_all_fifteen_pairs_complete_but_never_claim_global_pass(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            result = aggregate(directory)
        self.assertTrue(result["finite_matrix_comparison_complete"], result["errors"])
        self.assertFalse(result["passed"])
        self.assertEqual(result["completed_pair_count"], 15)
        self.assertEqual(result["status"], "INCOMPLETE_NO_PRUNING_COVERAGE_ALL_SHARDS_COMPLETED")

    def test_missing_shadow_mode_is_incomplete(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            (directory / "W5_SHADOW_MATRIX_J14_NODE4_MODEshadow_no_prune_RUN1.json").unlink()
            result = aggregate(directory)
        self.assertFalse(result["finite_matrix_comparison_complete"])
        self.assertEqual(len(result["missing_mode_shards"]), 1)

    def test_non_finite_comparison_cannot_pass(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            bad = report(10, 0, "shadow_no_prune", amplitude=float("nan"))
            self.write(directory, bad)
            result = aggregate(directory)
        self.assertFalse(result["finite_matrix_comparison_complete"])
        self.assertTrue(any("duplicate shard identity" in e for e in result["errors"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
