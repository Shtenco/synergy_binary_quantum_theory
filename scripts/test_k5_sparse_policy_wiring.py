#!/usr/bin/env python3
"""Focused regression tests for the K5 producer's explicit pruning modes."""
from __future__ import annotations

import inspect
import unittest
from unittest.mock import patch

import numpy as np
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


    def test_primitive_no_prune_retains_subthreshold_local_channel(self):
        fake_branch = (((1,),), (np.array([1.0]),), 1.0 + 0.0j)
        patches = [
            patch.object(k5, "initial_factorized_oriented", return_value=fake_branch),
            patch.object(k5, "VERT", (0,)),
            patch.object(k5, "local_spins", return_value=()),
            patch.object(k5, "allowed_k2_t", return_value=(0,)),
            patch.object(k5, "oriented_intertwiner", return_value=np.array([1.0])),
            patch.object(k5.np, "vdot", return_value=1e-13 + 0.0j),
        ]
        with patches[0], patches[1], patches[2], patches[3], patches[4], patches[5]:
            strict = k5.apply_primitive_oriented(("input",), (), 5, "no_prune")
            legacy = k5.apply_primitive_oriented(("input",), (), 5, "legacy")
        self.assertEqual(len(strict), 1)
        self.assertAlmostEqual(next(iter(strict.values())).real, 1e-13)
        self.assertEqual(legacy, {})

    def test_no_prune_hh_composition_aggregates_before_threshold(self):
        state = {"a": 1.0 + 0.0j, "b": 1.0 + 0.0j}
        def fake_h(single_state, target_v, jmax2, mode="legacy"):
            key = next(iter(single_state))
            return {"coherent": 0.75e-8 if key == "a" else -0.25e-8}
        with patch.object(k5, "apply_H_cached_state", side_effect=fake_h):
            strict = k5.compose_on_sparse(state, 0, 5, "no_prune")
            legacy = k5.compose_on_sparse(state, 0, 5, "legacy")
        self.assertAlmostEqual(strict["coherent"].real, 0.5e-8)
        self.assertEqual(legacy, {})

    def test_cache_separates_legacy_and_no_prune_results(self):
        k5.T_cached.cache_clear()
        def fake_apply(key, v, a, b, c, jmax2, adj=False, mode="legacy"):
            return {"mode": 1.0 if mode == "legacy" else 2.0}
        with patch.object(k5, "apply_T_oriented_basis", side_effect=fake_apply):
            legacy = dict(k5.T_cached(("key",), 0, 1, 2, 3, 5, False, "legacy"))
            strict = dict(k5.T_cached(("key",), 0, 1, 2, 3, 5, False, "no_prune"))
        self.assertEqual(legacy, {"mode": 1.0})
        self.assertEqual(strict, {"mode": 2.0})
        k5.T_cached.cache_clear()

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
