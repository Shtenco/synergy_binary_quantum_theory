#!/usr/bin/env python3
"""Synthetic Git-repository tests for the deduplicated LZMA branch snapshot."""
from __future__ import annotations

import hashlib
import json
import subprocess
import tarfile
import tempfile
import unittest
from pathlib import Path

import bqg_branch_snapshot_lzma as snapshot


def git(cwd: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", *args], cwd=cwd, check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    return result.stdout.strip()


class BranchSnapshotLzmaTests(unittest.TestCase):
    def make_repo(self, root: Path) -> Path:
        repo = root / "repo"
        repo.mkdir()
        git(repo, "init", "-b", "main")
        git(repo, "config", "user.email", "audit@example.invalid")
        git(repo, "config", "user.name", "Audit Fixture")
        (repo / "README.md").write_text("first\n", encoding="utf-8")
        (repo / "docs").mkdir()
        (repo / "docs" / "β.md").write_text("unicode path\n", encoding="utf-8")
        git(repo, "add", ".")
        git(repo, "commit", "-m", "first")
        first = git(repo, "rev-parse", "HEAD")
        git(repo, "update-ref", "refs/remotes/origin/main", first)
        git(repo, "update-ref", "refs/remotes/origin/alias", first)

        (repo / "README.md").write_text("second\n", encoding="utf-8")
        git(repo, "add", "README.md")
        git(repo, "commit", "-m", "second")
        second = git(repo, "rev-parse", "HEAD")
        git(repo, "update-ref", "refs/remotes/origin/feature", second)
        return repo

    def test_archive_deduplicates_blob_content_and_preserves_branch_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            repo = self.make_repo(root)
            archive_path = root / "snapshot.tar.xz"
            manifest_path = root / "manifest.json"
            result = snapshot.build_snapshot(repo, archive_path, manifest_path, batch_size=1)

            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(result["branch_count"], 3)
            self.assertEqual(result["unique_head_commit_count"], 2)
            self.assertEqual(manifest["unique_blob_count"], 3)
            self.assertEqual(
                manifest["commit_graph"][manifest["branch_heads"]["feature"]]["parents"],
                [manifest["branch_heads"]["main"]],
            )
            self.assertEqual(
                manifest["branch_heads"]["main"],
                manifest["branch_heads"]["alias"],
            )
            paths = {
                entry["path"]
                for entries in manifest["snapshots"].values()
                for entry in entries
            }
            self.assertIn("docs/β.md", paths)

            with tarfile.open(archive_path, mode="r:xz") as tar:
                members = tar.getnames()
                blobs = [name for name in members if name.startswith("blobs/")]
                self.assertEqual(len(blobs), manifest["unique_blob_count"])
                embedded = json.load(tar.extractfile("manifest.json"))
                self.assertEqual(embedded["branch_count"], 3)
                for name in blobs:
                    payload = tar.extractfile(name).read()
                    git_sha = name.split("/", 1)[1]
                    actual = hashlib.sha1(
                        b"blob " + str(len(payload)).encode("ascii") + b"\0" + payload
                    ).hexdigest()
                    self.assertEqual(actual, git_sha)

            self.assertEqual(
                result["archive_sha256"],
                hashlib.sha256(archive_path.read_bytes()).hexdigest(),
            )

    def test_invalid_controls_fail_before_git_access(self):
        with self.assertRaisesRegex(ValueError, "batch_size"):
            snapshot.build_snapshot(".", "x.tar.xz", "x.json", batch_size=0)
        with self.assertRaisesRegex(ValueError, "preset"):
            snapshot.build_snapshot(".", "x.tar.xz", "x.json", preset=10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
