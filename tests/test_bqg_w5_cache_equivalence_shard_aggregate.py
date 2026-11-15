import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from bqg_w5_cache_equivalence_shard_aggregate import aggregate


def report(jmax2, node, seed_k, mode, *, amplitude=1.0, status="PASS_ADJOINT_CONSISTENCY"):
    return {
        "status": status,
        "source_nodes": [node],
        "seed_Ks": [seed_k],
        "Jmax2": jmax2,
        "pruning_mode": "production",
        "cache_mode": mode,
        "max_states": None,
        "max_depth": None,
        "basis_keys_repr": ["basis-0"],
        "forward_B_columns_sparse": [[[0, amplitude, 0.0]]],
        "direct_Bdag_columns_sparse": [[[0, amplitude, 0.0]]],
        "lower_level_pruning_ledger": {
            "cache_policy": mode,
            "coverage": {
                "static_threshold_inventory_complete": True,
                "fixed_operator_mask_error_bound_complete": True,
            },
            "event_count": 0,
        },
        "coefficient_loss_ledger": {},
        "coherent_node_sum_ledger": [],
        "spectral_sensitivity": {},
        "diagnostics": {},
    }


class TestW5CachePairedShardAggregate(unittest.TestCase):
    def write(self, directory, data, suffix="RUN1"):
        j, n, k = data["Jmax2"], data["source_nodes"][0], data["seed_Ks"][0]
        mode = data["cache_mode"]
        path = directory / f"W5_CACHE_MATRIX_J{j}_NODE{n}_K{k}_MODE{mode}_{suffix}.json"
        path.write_text(json.dumps(data), encoding="utf-8")

    def fill(self, directory):
        for j in (10, 12, 14):
            for n in range(5):
                for k in (0, 2, 4, 6, 8):
                    for mode in ("on", "off"):
                        self.write(directory, report(j, n, k, mode))

    def test_all_75_pairs_pass_with_separate_mode_artifacts(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            result = aggregate(directory)
        self.assertTrue(result["passed"], result["errors"])
        self.assertEqual(result["completed_pair_count"], 75)
        self.assertEqual(result["valid_mode_shard_count"], 150)

    def test_missing_cache_off_shard_is_incomplete(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            (directory / "W5_CACHE_MATRIX_J14_NODE4_K8_MODEoff_RUN1.json").unlink()
            result = aggregate(directory)
        self.assertFalse(result["passed"])
        self.assertEqual(len(result["missing_mode_shards"]), 1)

    def test_matrix_mismatch_fails_pair(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            self.write(directory, report(10, 0, 0, "off", amplitude=2.0))
            result = aggregate(directory)
        self.assertFalse(result["passed"])
        self.assertTrue(any("duplicate shard identity" in e for e in result["errors"]))

    def test_incomplete_operator_gate_cannot_pass(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            self.write(directory, report(10, 0, 0, "on", status="INCOMPLETE_CLOSURE_LIMIT_REACHED"))
            result = aggregate(directory)
        self.assertFalse(result["passed"])
        self.assertTrue(any("failed checks" in e for e in result["errors"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
