#!/usr/bin/env python3
"""Create a deduplicated LZMA snapshot of every fetched Git branch head.

The archive contains one copy of each unique Git blob plus a manifest mapping
all branch heads and commit snapshots to their original paths and blob SHAs.
Run inside a checkout after fetching refs/remotes/origin/*.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import io
import json
import subprocess
import tarfile
from pathlib import Path


def git_output(arguments, input_bytes=None, cwd=None) -> bytes:
    result = subprocess.run(
        ["git", *arguments],
        cwd=cwd,
        input=input_bytes,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode:
        message = result.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"git {' '.join(arguments)} failed: {message}")
    return result.stdout


def branch_heads(cwd=None) -> dict[str, str]:
    raw = git_output([
        "for-each-ref", "--format=%(refname)%09%(objectname)",
        "refs/remotes/origin/",
    ], cwd=cwd)
    branches = {}
    for line in raw.decode("utf-8", errors="replace").splitlines():
        if not line:
            continue
        ref, sha = line.split("\t", 1)
        if ref.endswith("/HEAD"):
            continue
        branches[ref.removeprefix("refs/remotes/origin/")] = sha
    if "main" not in branches:
        raise RuntimeError("refs/remotes/origin/main is missing; fetch all branch heads first")
    return dict(sorted(branches.items()))


def git_commit_metadata(commit: str, cwd=None) -> dict:
    raw = git_output(["cat-file", "-p", commit], cwd=cwd)
    header = raw.split(b"\n\n", 1)[0].decode("utf-8", errors="replace")
    metadata = {"parents": []}
    for line in header.splitlines():
        if line.startswith("tree "):
            metadata["tree_sha"] = line[5:]
        elif line.startswith("parent "):
            metadata["parents"].append(line[7:])
    return metadata


def git_tree(commit: str, cwd=None) -> list[dict]:
    raw = git_output(["ls-tree", "-r", "-z", commit], cwd=cwd)
    entries = []
    for record in raw.split(b"\0"):
        if not record:
            continue
        metadata, path_bytes = record.split(b"\t", 1)
        mode, kind, sha = metadata.decode("ascii").split()
        item = {
            "mode": mode,
            "type": kind,
            "path": path_bytes.decode("utf-8", errors="replace"),
            "path_bytes_base64": base64.b64encode(path_bytes).decode("ascii"),
        }
        if kind == "blob":
            item["blob_sha"] = sha
        elif kind == "commit":
            item["gitlink_commit_sha"] = sha
        entries.append(item)
    return entries


def git_blob_batch(shas: list[str], cwd=None) -> dict[str, bytes]:
    if not shas:
        return {}
    request = "".join(f"{sha}\n" for sha in shas).encode("ascii")
    raw = git_output(["cat-file", "--batch"], request, cwd=cwd)
    blobs = {}
    offset = 0
    for sha in shas:
        line_end = raw.find(b"\n", offset)
        if line_end < 0:
            raise RuntimeError(f"truncated git cat-file header for {sha}")
        header = raw[offset:line_end].split()
        if len(header) != 3 or header[1] != b"blob":
            raise RuntimeError(f"git object {sha} is missing or not a blob")
        size = int(header[2])
        start = line_end + 1
        finish = start + size
        payload = raw[start:finish]
        if len(payload) != size or raw[finish:finish + 1] != b"\n":
            raise RuntimeError(f"truncated git blob payload for {sha}")
        actual = hashlib.sha1(
            b"blob " + str(size).encode("ascii") + b"\0" + payload
        ).hexdigest()
        if actual != sha:
            raise RuntimeError(f"Git blob SHA mismatch for {sha}: got {actual}")
        blobs[sha] = payload
        offset = finish + 1
    if offset != len(raw):
        raise RuntimeError("unexpected trailing bytes from git cat-file --batch")
    return blobs


def add_bytes(archive: tarfile.TarFile, name: str, payload: bytes) -> None:
    info = tarfile.TarInfo(name)
    info.size = len(payload)
    info.mode = 0o644
    info.mtime = 0
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    archive.addfile(info, io.BytesIO(payload))


def build_snapshot(
    cwd: str | Path,
    archive_path: str | Path,
    manifest_path: str | Path,
    *,
    batch_size: int = 64,
    preset: int = 6,
) -> dict:
    if not isinstance(batch_size, int) or isinstance(batch_size, bool) or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")
    if not isinstance(preset, int) or isinstance(preset, bool) or not 0 <= preset <= 9:
        raise ValueError("LZMA preset must be an integer from 0 through 9")

    branches = branch_heads(cwd)
    commits = sorted(set(branches.values()))
    snapshots = {}
    commit_graph = {}
    unique_blobs = set()
    total_path_entries = 0
    gitlinks = 0
    for commit in commits:
        entries = git_tree(commit, cwd)
        snapshots[commit] = entries
        commit_graph[commit] = git_commit_metadata(commit, cwd)
        total_path_entries += len(entries)
        for entry in entries:
            if entry["type"] == "blob":
                unique_blobs.add(entry["blob_sha"])
            elif entry["type"] == "commit":
                gitlinks += 1

    manifest = {
        "schema": "bqg-deduplicated-branch-snapshot/v1",
        "compression": "tar.xz (LZMA2)",
        "lzma_preset": preset,
        "branch_count": len(branches),
        "unique_head_commit_count": len(commits),
        "main_head": branches["main"],
        "branch_heads": branches,
        "commit_graph": commit_graph,
        "snapshot_commit_count": len(snapshots),
        "snapshot_path_entry_count": total_path_entries,
        "gitlink_entry_count": gitlinks,
        "unique_blob_count": len(unique_blobs),
        "blob_storage": "blobs/<git-blob-sha>; each unique blob is stored once",
        "submodule_note": "gitlink commit IDs are recorded but submodule trees are not dereferenced",
        "git_lfs_note": "tracked Git LFS pointer blobs are included; external LFS payloads are not fetched",
        "snapshots": snapshots,
    }
    manifest_bytes = (
        json.dumps(manifest, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        + "\n"
    ).encode("utf-8")

    archive_path = Path(archive_path)
    manifest_path = Path(manifest_path)
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_bytes(manifest_bytes)

    blob_list = sorted(unique_blobs)
    with tarfile.open(archive_path, mode="w:xz", preset=preset) as archive:
        add_bytes(archive, "manifest.json", manifest_bytes)
        for start in range(0, len(blob_list), batch_size):
            batch = blob_list[start:start + batch_size]
            payloads = git_blob_batch(batch, cwd)
            for sha in batch:
                add_bytes(archive, f"blobs/{sha}", payloads[sha])

    digest = hashlib.sha256()
    with archive_path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {
        "archive": str(archive_path),
        "manifest": str(manifest_path),
        "archive_bytes": archive_path.stat().st_size,
        "manifest_bytes": len(manifest_bytes),
        "archive_sha256": digest.hexdigest(),
        "branch_count": len(branches),
        "unique_head_commit_count": len(commits),
        "snapshot_path_entry_count": total_path_entries,
        "unique_blob_count": len(unique_blobs),
        "lzma_preset": preset,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        default=Path("verification_results/BQG_BRANCH_SNAPSHOT_2026-11-11.tar.xz"))
    parser.add_argument("--manifest-out", type=Path,
                        default=Path("verification_results/BQG_BRANCH_SNAPSHOT_2026-11-11.manifest.json"))
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--preset", type=int, default=6)
    args = parser.parse_args()
    result = build_snapshot(".", args.out, args.manifest_out,
                            batch_size=args.batch_size, preset=args.preset)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
