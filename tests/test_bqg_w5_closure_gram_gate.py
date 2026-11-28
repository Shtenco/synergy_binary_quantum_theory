#!/usr/bin/env python3
import unittest
import sympy as sp

from bqg_w5_closure_gram_gate import (
    allowed_k2, gram_certificate, validate_basis_key,
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
        self.assertEqual(result["input_metric_Q"], [
            ["1", "0", "0", "0", "0"],
            ["0", "1", "0", "0", "0"],
            ["0", "0", "1", "0", "0"],
            ["0", "0", "0", "1", "0"],
            ["0", "0", "0", "0", "1"],
        ])

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


if __name__ == "__main__":
    unittest.main()
