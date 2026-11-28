#!/usr/bin/env python3
import json
import tempfile
import unittest
from pathlib import Path
from bqg_w5_closure_gram_gate import (
    allowed_k2, gram_certificate, load_keys, validate_basis_key,
)


class W5ClosureGramGateTests(unittest.TestCase):
    def seed_key(self, k2):
        return ((4,) * 10, (k2, 0, 0, 0, 0))

    def test_five_seed_keys_are_admissible_and_orthonormal(self):
        keys = [self.seed_key(k2) for k2 in (0, 2, 4, 6, 8)]
        result = gram_certificate(keys)
        self.assertTrue(result["passed"])
        self.assertEqual(result["basis_dimension"], 5)
        self.assertEqual(result["gram_rank"], 5)
        self.assertEqual(result["input_metric_Q"],
                         {"representation": "implicit_identity", "dimension": 5})
        self.assertEqual(result["gram_matrix"],
                         {"type": "identity", "dimension": 5})

    def test_duplicate_keys_fail_closed(self):
        key = self.seed_key(0)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            gram_certificate([key, key])

    def test_inadmissible_node_channel_fails_closed(self):
        spins = (4,) * 10
        # At node 0, local spins (4,4,4,4) allow K2=0; change one incident
        # edge so that pair (4,2) and pair (4,4) have no common K2=0 channel.
        changed = list(spins)
        from bqg_w5_closure_gram_gate import EDGE_INDEX
        changed[EDGE_INDEX[(0, 4)]] = 2
        key = (tuple(changed), (0, 0, 0, 0, 0))
        with self.assertRaisesRegex(ValueError, "inadmissible"):
            validate_basis_key(key)

    def test_allowed_channels_match_doubled_spin_triangle_rule(self):
        self.assertEqual(allowed_k2((4, 4, 4, 4)), (0, 2, 4, 6, 8))

    def test_complete_artifact_key_parser(self):
        keys = [self.seed_key(k2) for k2 in (0, 2, 4, 6, 8)]
        payload = {"basis_keys_repr": [repr(key) for key in keys]}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "complete.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            parsed, scope = load_keys(path)
        self.assertEqual(parsed, keys)
        self.assertEqual(scope, "complete_basis_list")

    def test_partial_artifact_parser_keeps_partial_scope(self):
        key = self.seed_key(0)
        payload = {
            "status": "INCOMPLETE_CLOSURE_LIMIT_REACHED",
            "columns_computed": 1,
            "partial_forward_columns_top_amplitudes": [
                {"input_key": repr(key), "support_size": 0}
            ]
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "partial.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            parsed, scope = load_keys(path)
        self.assertEqual(parsed, [key])
        self.assertEqual(scope, "computed_input_columns_only")

    def test_partial_status_cannot_masquerade_as_complete_basis(self):
        key = self.seed_key(0)
        payload = {
            "status": "INCOMPLETE_CLOSURE_LIMIT_REACHED",
            "basis_keys_repr": [repr(key)],
            "closure_size": 1,
            "forward_B_columns_sparse": [[]],
            "direct_Bdag_columns_sparse": [[]],
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "mislabelled.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "completed closure"):
                load_keys(path)

    def test_complete_artifact_requires_matching_matrix_dimensions(self):
        key = self.seed_key(0)
        payload = {
            "status": "PASS_ADJOINT_CONSISTENCY",
            "closure_size": 1,
            "basis_keys_repr": [repr(key)],
            "forward_B_columns_sparse": [],
            "direct_Bdag_columns_sparse": [[]],
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad_complete.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "forward-column count"):
                load_keys(path)

    def test_partial_artifact_requires_matching_column_count(self):
        key = self.seed_key(0)
        payload = {
            "status": "INCOMPLETE_CLOSURE_LIMIT_REACHED",
            "columns_computed": 2,
            "partial_forward_columns_top_amplitudes": [
                {"input_key": repr(key), "support_size": 0}
            ],
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad_partial.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "columns_computed"):
                load_keys(path)


if __name__ == "__main__":
    unittest.main()
