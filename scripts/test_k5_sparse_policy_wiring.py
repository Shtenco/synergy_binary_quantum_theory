#!/usr/bin/env python3
"""Focused regression tests for the K5 producer's explicit pruning modes."""
from __future__ import annotations

import inspect
import unittest

import k5_peter_weyl_safe_hda_column as k5


class K5SparsePolicyTests(unittest.TestCase):
    def test_add_dict_no_prune_retains_tiny_nonzero(self):
        dst = {}
        k5.add_dict(dst, {"tiny": 1e-14}, mode="no_prune")
        self.assertEqual(dst["tiny"], 1e-14)

    def test_legacy_add_dict_preserves_historical_cutoff(self):
        dst = {}
        k5.add_dict(dst, {"tiny": 1e-14}, mode="legacy")
        self.assertNotIn("tiny", dst)

    def test_no_prune_state_bypasses_state_threshold(self):
        state = {"tiny": 1e-14, "large": 1.0}
        self.assertEqual(k5.prune_state(state, 1e-8, "no_prune"), state)

    def test_legacy_state_pruning_remains_available(self):
        state = {"tiny": 1e-14, "large": 1.0}
        self.assertEqual(k5.prune_state(state, 1e-8, "legacy"), {"large": 1.0})

    def test_policy_is_part_of_t_cache_key(self):
        params = inspect.signature(k5.T_cached).parameters
        self.assertIn("mode", params)
        self.assertIsNotNone(k5.T_cached.cache_info)

    def test_raw_volume_hermiticity_is_recorded_before_symmetrization(self):
        k5.VOLUME_RAW_DEFECTS.clear()
        matrix = k5.volume123_matrix(1, 1, 1)
        self.assertEqual(matrix.shape, (8, 8))
        record = k5.VOLUME_RAW_DEFECTS[(1, 1, 1)]
        self.assertIn("relative_antihermitian_defect", record)
        self.assertGreaterEqual(record["relative_antihermitian_defect"], 0.0)
        self.assertLess(record["relative_antihermitian_defect"], 1e-12)


if __name__ == "__main__":
    unittest.main()
