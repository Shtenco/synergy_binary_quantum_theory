from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import collective_l1_amplitude_producer as producer


class L1AmplitudeProducerTests(unittest.TestCase):
    def test_intermediate_pruning_counterexample(self):
        result = producer.self_test()
        self.assertTrue(result["passed"], result)

    def test_stable_sum_is_shard_order_invariant(self):
        x = {b"a": 1.0 + 2.0j, b"b": -0.25j}
        y = {b"a": -0.5 - 1.0j, b"c": 3.0}
        a, _, _ = producer.stable_sum_states([x, y])
        b, _, _ = producer.stable_sum_states([y, x])
        self.assertEqual(a, b)

    def test_discarded_global_norm_is_reported(self):
        tiny = {b"a": 0.5e-10, b"b": 0.5e-10}
        kept, discarded, relative = producer.stable_sum_states([tiny])
        self.assertEqual(kept, {})
        self.assertGreater(discarded, 0.0)
        self.assertAlmostEqual(relative, 1.0)

    def test_duplicate_sparse_state_keys_fail_closed(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "duplicate.npz"
            np.savez_compressed(
                p,
                spins=np.zeros((2, 2), dtype=np.int16),
                Ks=np.zeros((2, 3), dtype=np.int16),
                amp=np.asarray([1.0 + 0j, 2.0 + 0j], dtype=np.complex128),
            )
            with self.assertRaisesRegex(ValueError, "duplicate state key"):
                producer.load_state(p, 2, 3)

    def test_nonfinite_amplitude_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "nonfinite.npz"
            np.savez_compressed(
                p,
                spins=np.zeros((1, 2), dtype=np.int16),
                Ks=np.zeros((1, 3), dtype=np.int16),
                amp=np.asarray([complex(float("nan"), 0.0)], dtype=np.complex128),
            )
            with self.assertRaisesRegex(ValueError, "non-finite"):
                producer.load_state(p, 2, 3)


if __name__ == "__main__":
    unittest.main()
