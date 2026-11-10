from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import pl_lorentzian_48_collect as COLLECT


class PL48FailClosedCollectorTests(unittest.TestCase):
    def test_accepts_well_formed_empty_sparse_state(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "empty.npz"
            np.savez_compressed(
                path,
                spins=np.zeros((0, 32), dtype=np.int16),
                Ks=np.zeros((0, 16), dtype=np.int16),
                amp=np.zeros((0,), dtype=np.complex128),
            )
            self.assertEqual(COLLECT.load_state(path), {})

    def test_rejects_nonfinite_amplitudes_before_pruning(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "nan.npz"
            np.savez_compressed(
                path,
                spins=np.zeros((1, 32), dtype=np.int16),
                Ks=np.zeros((1, 16), dtype=np.int16),
                amp=np.asarray([complex(float("nan"), 0.0)]),
            )
            with self.assertRaisesRegex(ValueError, "non-finite amplitudes"):
                COLLECT.load_state(path)

    def test_rejects_duplicate_sparse_keys(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.npz"
            np.savez_compressed(
                path,
                spins=np.zeros((2, 32), dtype=np.int16),
                Ks=np.zeros((2, 16), dtype=np.int16),
                amp=np.asarray([1.0 + 0.0j, 2.0 + 0.0j]),
            )
            with self.assertRaisesRegex(ValueError, "duplicate sparse-state key"):
                COLLECT.load_state(path)

    def test_pruning_refuses_nan_instead_of_silently_dropping_it(self):
        state = {((1,) * 32, (0,) * 16): complex(float("nan"), 0.0)}
        with self.assertRaisesRegex(ValueError, "refusing to prune"):
            COLLECT.prune(state)

    def test_accumulation_refuses_nonfinite_amplitudes(self):
        state = {((1,) * 32, (0,) * 16): complex(float("inf"), 0.0)}
        with self.assertRaisesRegex(ValueError, "non-finite input amplitude"):
            COLLECT.add({}, state)


if __name__ == "__main__":
    unittest.main()
