import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from bqg_w5_cache_equivalence_shard_aggregate import aggregate


class TestW5CacheShardAggregate(unittest.TestCase):
    def report(self, jmax2, node, seed_k, mode, *, matrix_value=1.0, ledger_value=1):
        return {
            "status": "PASS_ADJOINT_CONSISTENCY",
            "algebraic_status": "PASS_ADJOINT_CONSISTENCY",
            "lower_level_pruning_coverage_complete": False,
            "source_nodes": [node],
            "pruning_mode": "production",
            "cache_mode": mode,
            "Jmax2": jmax2,
            "seed_count": 1,
            "seed_Ks": [seed_k],
            "max_states": None,
            "max_depth": None,
            "unbounded_closure": True,
            "omitted_new_support_keys": 0,
            "closure_size": 2,
            "basis_id_sha256": hashlib.sha256(b"basis-0\\nbasis-1").hexdigest(),
            "basis_keys_repr": ["basis-0", "basis-1"],
            "forward_B_columns_sparse": [[[0, matrix_value, 0.0]], [[1, 2.0, 0.0]]],
            "direct_Bdag_columns_sparse": [[[0, 3.0, 0.0]], [[1, 4.0, 0.0]]],
            "coefficient_loss_ledger": {"amp_tol": {"count": ledger_value}},
            "coherent_node_sum_ledger": [{"input_key": "k", "forward": {"norm": 1.0}}],
            "lower_level_pruning_ledger": {
                "event_count": 1,
                "events": [{"site": "site-a", "kind": "kept", "amplitude": 1.0}],
                "by_site": {"site-a": {"coefficient_count": 1}},
                "by_site_and_kind": {"site-a::kept": {"coefficient_count": 1}},
                "coverage": {"unclassified_sites": 0},
                "cache_policy": {"mode": mode},
            },
            "diagnostics": {"direct_adjoint_relative_residual": 0.0},
        }

    def write_pair(self, directory, jmax2, node, seed_k, *, bad_matrix=False, bad_ledger=False):
        for mode in ("on", "off"):
            report = self.report(
                jmax2, node, seed_k, mode,
                matrix_value=(1.1 if bad_matrix and mode == "off" else 1.0),
                ledger_value=(2 if bad_ledger and mode == "off" else 1),
            )
            path = directory / f"W5_CACHE_MATRIX_J{jmax2}_NODE{node}_K{seed_k}_MODE{mode}_RUN1.json"
            path.write_text(json.dumps(report), encoding="utf-8")

    def fill(self, directory, *, bad_pair=None, bad_matrix=False, bad_ledger=False):
        for node in range(5):
            for seed_k in (0, 2, 4, 6, 8):
                key = (10, node, seed_k)
                self.write_pair(
                    directory, *key,
                    bad_matrix=(key == bad_pair and bad_matrix),
                    bad_ledger=(key == bad_pair and bad_ledger),
                )

    def test_all_25_pairs_match(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            result = aggregate(directory, (10,))
        self.assertTrue(result["passed"], result["errors"])
        self.assertEqual(result["completed_pair_count"], 25)
        self.assertEqual(result["valid_mode_shard_count"], 50)

    def test_missing_mode_shard_is_incomplete(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            (directory / "W5_CACHE_MATRIX_J10_NODE4_K8_MODEoff_RUN1.json").unlink()
            result = aggregate(directory, (10,))
        self.assertFalse(result["passed"])
        self.assertEqual(len(result["missing_shards"]), 1)

    def test_matrix_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory, bad_pair=(10, 0, 0), bad_matrix=True)
            result = aggregate(directory, (10,))
        self.assertFalse(result["passed"])
        self.assertTrue(any("matrix residual exceeds" in error for error in result["errors"]))

    def test_ledger_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory, bad_pair=(10, 0, 0), bad_ledger=True)
            result = aggregate(directory, (10,))
        self.assertFalse(result["passed"])
        self.assertTrue(any("coefficient loss ledgers differ" in error for error in result["errors"]))


    def test_incorrect_basis_digest_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            path = directory / "W5_CACHE_MATRIX_J10_NODE0_K0_MODEon_RUN1.json"
            report = json.loads(path.read_text())
            report["basis_id_sha256"] = "a" * 64
            path.write_text(json.dumps(report), encoding="utf-8")
            result = aggregate(directory, (10,))
        self.assertFalse(result["passed"])
        self.assertTrue(any("failed checks: basis_hash" in error for error in result["errors"]))

    def test_matrix_column_count_must_match_basis_dimension(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            path = directory / "W5_CACHE_MATRIX_J10_NODE0_K0_MODEon_RUN1.json"
            report = json.loads(path.read_text())
            report["forward_B_columns_sparse"].pop()
            path.write_text(json.dumps(report), encoding="utf-8")
            result = aggregate(directory, (10,))
        self.assertFalse(result["passed"])
        self.assertTrue(any("failed checks: forward_B_columns_sparse" in error for error in result["errors"]))

    def test_bounded_closure_cannot_pass(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            path = directory / "W5_CACHE_MATRIX_J10_NODE0_K0_MODEon_RUN1.json"
            report = json.loads(path.read_text())
            report["unbounded_closure"] = False
            report["max_states"] = 500
            path.write_text(json.dumps(report), encoding="utf-8")
            result = aggregate(directory, (10,))
        self.assertFalse(result["passed"])
        self.assertTrue(any("failed checks: unbounded_closure" in error for error in result["errors"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
