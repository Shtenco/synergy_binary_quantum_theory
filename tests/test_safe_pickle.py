import os
import pickle
import tempfile
import unittest
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "proof_tools"))
from safe_pickle import safe_load_bytes, safe_load_builtin_pickle, safe_load_numpy_pickle


MARKER_PATH = Path(tempfile.gettempdir()) / "bqg_pickle_should_not_execute"


class _MaliciousPayload:
    def __reduce__(self):
        return (os.system, (f"touch {MARKER_PATH}",))


class SafePickleTests(unittest.TestCase):
    def test_loads_plain_metadata_without_numpy_globals(self):
        payload = {
            "kind": "test",
            "blocks": {3: {"q_rows": {(0, 2, 4): 7}, "m": 2}},
            "shard": 4,
        }
        raw = pickle.dumps(payload, protocol=5)
        self.assertEqual(safe_load_bytes(raw), payload)

    def test_loads_numpy_numeric_arrays_protocol_four_and_five(self):
        payload = {
            (0, 2, 4): np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64),
            (2, 0, 2): np.array([[1, 2]], dtype=np.int64),
        }
        for protocol in (4, 5):
            with self.subTest(protocol=protocol):
                raw = pickle.dumps(payload, protocol=protocol)
                loaded = safe_load_bytes(raw, allow_numpy=True)
                self.assertEqual(set(loaded), set(payload))
                for key in payload:
                    np.testing.assert_array_equal(loaded[key], payload[key])

    def test_rejects_arbitrary_global_before_execution(self):
        marker = MARKER_PATH
        marker.unlink(missing_ok=True)
        raw = pickle.dumps(_MaliciousPayload(), protocol=5)
        with self.assertRaisesRegex(pickle.UnpicklingError, "forbidden pickle global"):
            safe_load_bytes(raw, allow_numpy=True)
        self.assertFalse(marker.exists())

    def test_rejects_object_dtype_arrays(self):
        raw = pickle.dumps({"bad": np.array([{"x": 1}], dtype=object)}, protocol=5)
        with self.assertRaises(pickle.UnpicklingError):
            safe_load_bytes(raw, allow_numpy=True)

    def test_builtin_loader_rejects_numpy_globals(self):
        raw = pickle.dumps({"array": np.array([1.0, 2.0])}, protocol=5)
        with self.assertRaises(pickle.UnpicklingError):
            safe_load_bytes(raw, allow_numpy=False)

    def test_file_helpers_use_same_restrictions(self):
        raw = pickle.dumps({"k": (1, 2), "v": [3, 4]}, protocol=5)
        with tempfile.TemporaryFile() as f:
            f.write(raw)
            f.seek(0)
            self.assertEqual(safe_load_builtin_pickle(f), {"k": (1, 2), "v": [3, 4]})
        raw = pickle.dumps(np.eye(2), protocol=5)
        with tempfile.TemporaryFile() as f:
            f.write(raw)
            f.seek(0)
            np.testing.assert_array_equal(safe_load_numpy_pickle(f), np.eye(2))


if __name__ == "__main__":
    unittest.main()
