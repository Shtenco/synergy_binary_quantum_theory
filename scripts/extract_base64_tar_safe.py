#!/usr/bin/env python3
"""Decode and safely extract a base64-encoded gzip tar proof bundle."""
from __future__ import annotations

import argparse
import base64
import io
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROOF_TOOLS = ROOT / "proof_tools"
import sys
if str(PROOF_TOOLS) not in sys.path:
    sys.path.insert(0, str(PROOF_TOOLS))

from safe_archive import safe_extract_tar


def extract_base64_tar(
    archive_path: str | Path,
    destination: str | Path,
    *,
    max_archive_bytes: int = 64 * 1024 * 1024,
    max_uncompressed_bytes: int = 256 * 1024 * 1024,
) -> Path:
    if max_archive_bytes <= 0 or max_uncompressed_bytes <= 0:
        raise ValueError("archive size limits must be positive")
    encoded = "".join(Path(archive_path).read_text(encoding="ascii").split())
    raw = base64.b64decode(encoded, validate=True)
    if len(raw) > max_archive_bytes:
        raise ValueError(
            f"compressed archive exceeds limit: {len(raw)} > {max_archive_bytes}"
        )
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as archive:
        total = sum(member.size for member in archive.getmembers() if member.isfile())
        if total > max_uncompressed_bytes:
            raise ValueError(
                f"uncompressed archive exceeds limit: {total} > {max_uncompressed_bytes}"
            )
        return safe_extract_tar(archive, destination)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--max-archive-bytes", type=int, default=64 * 1024 * 1024)
    parser.add_argument("--max-uncompressed-bytes", type=int, default=256 * 1024 * 1024)
    args = parser.parse_args()
    result = extract_base64_tar(
        args.archive,
        args.destination,
        max_archive_bytes=args.max_archive_bytes,
        max_uncompressed_bytes=args.max_uncompressed_bytes,
    )
    print(f"SAFE_ARCHIVE_EXTRACTED {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
