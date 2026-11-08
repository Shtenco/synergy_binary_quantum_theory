import hashlib
import re
import tempfile
import unittest
from pathlib import Path

import numpy as np


def _sha(value):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def validate_term_record(path, metadata, expected_source=0):
    """Isolated copy of the V2 48-term collector's term-record guard."""
    if not isinstance(metadata, dict) or metadata.get("passed") is not True:
        raise RuntimeError("V2_TERM_NOT_PASSED_BOOLEAN_TRUE")
    mode = metadata.get("mode")
    index = metadata.get("index")
    source = metadata.get("source_node")
    coefficient = metadata.get("PL_epsilon_coefficient")
    if mode not in ("forward", "adjoint"):
        raise RuntimeError("V2_TERM_BAD_MODE")
    if type(index) is not int or index not in range(24):
        raise RuntimeError("V2_TERM_BAD_INTEGER_INDEX")
    if type(source) is not int or source != expected_source:
        raise RuntimeError("V2_TERM_BAD_SOURCE_NODE")
    if type(coefficient) is not int or coefficient not in (-1, 1):
        raise RuntimeError("V2_TERM_BAD_ORIENTATION_COEFFICIENT")
    jmax = metadata.get("Jmax")
    if (not isinstance(jmax, (int, float)) or isinstance(jmax, bool)
            or float(jmax) != 3.5):
        raise RuntimeError("V2_TERM_BAD_CUTOFF")
    checks = metadata.get("checks")
    if (not isinstance(checks, dict) or not checks
            or not all(value is True for value in checks.values())):
        raise RuntimeError("V2_TERM_CHECKS_NOT_BOOLEAN_TRUE")
    npz_path = Path(path).with_suffix(".npz")
    digest = metadata.get("npz_sha256")
    if not _sha(digest) or not npz_path.is_file():
        raise RuntimeError("V2_TERM_NPZ_DIGEST_MISSING")
    if hashlib.sha256(npz_path.read_bytes()).hexdigest() != digest:
        raise RuntimeError("V2_TERM_NPZ_DIGEST_MISMATCH")
    with np.load(npz_path, allow_pickle=False) as data:
        if set(data.files) != {"spins", "Ks", "amp"}:
            raise RuntimeError("V2_TERM_NPZ_BAD_KEYS")
        spins, Ks, amp = data["spins"], data["Ks"], data["amp"]
        if (spins.shape != (len(amp), 32) or Ks.shape != (len(amp), 16)
                or amp.ndim != 1 or spins.dtype.kind not in "iu"
                or Ks.dtype.kind not in "iu" or amp.dtype.kind != "c"
                or not np.isfinite(amp).all()):
            raise RuntimeError("V2_TERM_NPZ_BAD_SHAPE_DTYPE_OR_AMPLITUDE")
        support = metadata.get("gauss_support")
        if type(support) is not int or support != len(amp):
            raise RuntimeError("V2_TERM_NPZ_SUPPORT_MISMATCH")
        seen = set()
        for spin_row, K_row in zip(spins, Ks):
            key = tuple(map(int, spin_row)), tuple(map(int, K_row))
            if key in seen:
                raise RuntimeError("V2_TERM_NPZ_DUPLICATE_KEY")
            seen.add(key)
    return mode, index, source


class V2ArtifactGuardAudit(unittest.TestCase):
    def _record(self, root, **overrides):
        path = Path(root) / "term_forward_0.json"
        npz = path.with_suffix(".npz")
        np.savez_compressed(npz, spins=np.array([[1] * 32]),
                            Ks=np.array([[0] * 16]), amp=np.array([1 + 0j]))
        metadata = {
            "passed": True, "mode": "forward", "index": 0, "source_node": 0,
            "PL_epsilon_coefficient": 1, "Jmax": 3.5, "gauss_support": 1,
            "checks": {"finite": True},
            "npz_sha256": hashlib.sha256(npz.read_bytes()).hexdigest(),
        }
        metadata.update(overrides)
        return path, metadata

    def test_truthy_passed_string_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path, meta = self._record(d, passed="yes")
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_NOT_PASSED_BOOLEAN_TRUE"):
                validate_term_record(path, meta)

    def test_boolean_index_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path, meta = self._record(d, index=True)
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_BAD_INTEGER_INDEX"):
                validate_term_record(path, meta)

    def test_boolean_orientation_coefficient_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path, meta = self._record(d, PL_epsilon_coefficient=True)
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_BAD_ORIENTATION_COEFFICIENT"):
                validate_term_record(path, meta)

    def test_npz_digest_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path, meta = self._record(d, npz_sha256="0" * 64)
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_NPZ_DIGEST_MISMATCH"):
                validate_term_record(path, meta)

    def test_wrong_npz_shape_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path, meta = self._record(d)
            npz = path.with_suffix(".npz")
            np.savez_compressed(npz, spins=np.array([[1] * 31]),
                                Ks=np.array([[0] * 16]), amp=np.array([1 + 0j]))
            meta["npz_sha256"] = hashlib.sha256(npz.read_bytes()).hexdigest()
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_NPZ_BAD_SHAPE_DTYPE_OR_AMPLITUDE"):
                validate_term_record(path, meta)

    def test_duplicate_npz_keys_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path, meta = self._record(d)
            npz = path.with_suffix(".npz")
            np.savez_compressed(npz, spins=np.array([[1] * 32, [1] * 32]),
                                Ks=np.array([[0] * 16, [0] * 16]),
                                amp=np.array([1 + 0j, 2 + 0j]))
            meta["gauss_support"] = 2
            meta["npz_sha256"] = hashlib.sha256(npz.read_bytes()).hexdigest()
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_NPZ_DUPLICATE_KEY"):
                validate_term_record(path, meta)

    def test_boolean_support_count_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path, meta = self._record(d, gauss_support=True)
            with self.assertRaisesRegex(RuntimeError, "V2_TERM_NPZ_SUPPORT_MISMATCH"):
                validate_term_record(path, meta)

    def test_canonical_v2_artifact_is_accepted(self):
        with tempfile.TemporaryDirectory() as d:
            path, meta = self._record(d)
            self.assertEqual(validate_term_record(path, meta), ("forward", 0, 0))


if __name__ == "__main__":
    unittest.main(verbosity=2)
