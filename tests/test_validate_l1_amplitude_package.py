from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import validate_l1_amplitude_package as validator


class L1AmplitudePackageValidatorTests(unittest.TestCase):
    def _fixture(self, root: Path):
        arrays = {}
        for e in validator.EDGES:
            for prefix in ("u", "v"):
                stem = f"{prefix}_{e}"
                arrays[f"{stem}_spins"] = np.zeros((0, 768), np.int16)
                arrays[f"{stem}_Ks"] = np.zeros((0, 384), np.int16)
                arrays[f"{stem}_amp"] = np.zeros((0,), np.complex128)
            for n in validator.NODES:
                for prefix in ("u_shard", "v_shard"):
                    stem = f"{prefix}_{e}_{n}"
                    arrays[f"{stem}_spins"] = np.zeros((0, 768), np.int16)
                    arrays[f"{stem}_Ks"] = np.zeros((0, 384), np.int16)
                    arrays[f"{stem}_amp"] = np.zeros((0,), np.complex128)
        archive = root / "L1_AMPLITUDE_PACKAGE.npz"
        np.savez_compressed(archive, **arrays)
        shard_root = root / "shards"
        shard_root.mkdir()
        shard_hashes = {}
        for e in validator.EDGES:
            for n in validator.NODES:
                for name in (f"edge_{e}_node_{n}.json",
                             f"u_{e}_node_{n}.npz",
                             f"v_{e}_node_{n}.npz"):
                    item = shard_root / name
                    item.write_bytes((name + "\\n").encode())
                    shard_hashes[name] = validator.sha256(item)
        package = {
            "schema": "L1_AMPLITUDE_PACKAGE",
            "schema_version": 1,
            "package_status": "PASS_E_ONLY_REGRESSION_PRECURSOR",
            "science_status": "L1_EUCLIDEAN_E_ONLY_PRECURSOR_NOT_FULL_PHYSICS",
            "complete_physical_package": False,
            "missing_required_operators": ["S", "R"],
            "source_commit": "a" * 40,
            "source_commit_attestation": "caller-supplied; must equal git rev-parse HEAD at production time; not independently queried by this script",
            "amplitude_archive_sha256": validator.sha256(archive),
            "provenance": {
                "shard_file_sha256": shard_hashes,
                "shard_manifest_sha256": validator.canonical_sha(shard_hashes),
                "source_worker_blob_sha": "9b1845aa263bacf23a0f79d5140bb82ceb56c8a1",
                "source_collector_blob_sha": "3ba2514ecd553fc54a59d3b56c1c43ce24aae09f",
                "producer_source_sha256": validator.sha256(ROOT / "scripts/collective_l1_amplitude_producer.py"),
                "producer_blob_sha": validator.git_blob_sha(ROOT / "scripts/collective_l1_amplitude_producer.py"),
                "worker_identity_limitation": "worker/collector source blobs are pinned by this schema, but legacy shard metadata does not self-attest its generating source blob",
            },
            "operator_registry": {
                "E": {"source_blob_sha": "9b1845aa263bacf23a0f79d5140bb82ceb56c8a1",
                      "backend": "LocalPLPeterWeylEuclidean exact active-cone"},
                "S": {"status": "MISSING"}, "R": {"status": "MISSING"},
            },
            "carrier": {
                "nodes": 384, "dual_links": 768, "parent_coarse_tetra": 0,
                "parent_fine_nodes": 24, "edge_representatives": [0, 1, 5],
                "Jmax": 2.5, "basis_status": "NOT_YET_JOINT_E_S_R_W",
            },
            "checks": {
                "all_72_shards_present_and_passed": True,
                "source_states_identical": True,
                "finite_amplitudes_and_moments": True,
                "legacy_collector_state_regression": True,
                "legacy_collector_moment_regression": True,
                "S4_three_orbit_consistency": True,
                "nonnegative_depth2_variance_with_frozen_tolerance": True,
            },
        }
        package["package_sha256"] = validator.canonical_sha(package)
        manifest = root / "L1_AMPLITUDE_PACKAGE.json"
        manifest.write_text(json.dumps(package, indent=2), encoding="utf-8")
        return manifest, archive, shard_root

    def test_valid_integrity_fixture(self):
        with tempfile.TemporaryDirectory() as td:
            manifest, archive, shard_root = self._fixture(Path(td))
            result = validator.validate(manifest, archive, ROOT / "scripts/collective_l1_amplitude_producer.py", shard_root)
            self.assertTrue(result["passed"])
            self.assertFalse(result["physical_completeness"])
            self.assertEqual(result["checked_sparse_state_groups"], 150)

    def test_archive_tampering_is_detected(self):
        with tempfile.TemporaryDirectory() as td:
            manifest, archive, shard_root = self._fixture(Path(td))
            with archive.open("ab") as f:
                f.write(b"tamper")
            with self.assertRaisesRegex(ValueError, "archive SHA256 mismatch"):
                validator.validate(manifest, archive, ROOT / "scripts/collective_l1_amplitude_producer.py", shard_root)

    def test_source_shard_tampering_is_detected(self):
        with tempfile.TemporaryDirectory() as td:
            manifest, archive, shard_root = self._fixture(Path(td))
            target = shard_root / "edge_0_node_0.json"
            target.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "modified source shard"):
                validator.validate(manifest, archive, ROOT / "scripts/collective_l1_amplitude_producer.py", shard_root)

    def test_missing_required_check_is_detected(self):
        with tempfile.TemporaryDirectory() as td:
            manifest, archive, shard_root = self._fixture(Path(td))
            package = json.loads(manifest.read_text())
            package["checks"].pop("S4_three_orbit_consistency")
            body = dict(package)
            body.pop("package_sha256")
            package["package_sha256"] = validator.canonical_sha(body)
            manifest.write_text(json.dumps(package))
            with self.assertRaisesRegex(ValueError, "missing, unexpected, or failed"):
                validator.validate(manifest, archive, ROOT / "scripts/collective_l1_amplitude_producer.py", shard_root)


if __name__ == "__main__":
    unittest.main()
