#!/usr/bin/env python3
"""Download all matching GitHub Actions artifacts from a workflow run, with pagination.

Designed for large sharded BQG proof runs where artifact counts exceed 100.
Fails closed on duplicate names, incomplete expected counts, HTTP errors, or
malformed ZIP archives.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROOF_TOOLS = ROOT / "proof_tools"
if str(PROOF_TOOLS) not in sys.path:
    sys.path.insert(0, str(PROOF_TOOLS))

from safe_archive import safe_artifact_directory, safe_extract_zip

API = "https://api.github.com"


class SameOriginAuthRedirectHandler(urllib.request.HTTPRedirectHandler):
    """Never forward a GitHub bearer token to a different HTTPS origin."""

    @staticmethod
    def origin(url: str) -> tuple[str, str, int]:
        parsed = urllib.parse.urlsplit(url)
        port = parsed.port or (443 if parsed.scheme.lower() == "https" else 80)
        return parsed.scheme.lower(), (parsed.hostname or "").lower(), port

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        parsed = urllib.parse.urlsplit(newurl)
        if parsed.scheme.lower() != "https":
            raise urllib.error.HTTPError(
                req.full_url, code, "refusing non-HTTPS artifact redirect", headers, fp
            )
        redirected = super().redirect_request(req, fp, code, msg, headers, newurl)
        if redirected is not None and self.origin(req.full_url) != self.origin(newurl):
            redirected.remove_header("Authorization")
        return redirected


SAFE_OPENER = urllib.request.build_opener(SameOriginAuthRedirectHandler())


def request(url: str, token: str) -> bytes:
    req = urllib.request.Request(
        url,
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "BQG-artifact-downloader",
        },
    )
    with SAFE_OPENER.open(req, timeout=120) as r:
        return r.read()

def download_to_file(url: str, token: str, destination: Path) -> tuple[int, str]:
    req = urllib.request.Request(
        url,
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "BQG-artifact-downloader",
        },
    )
    digest = hashlib.sha256()
    size = 0
    with SAFE_OPENER.open(req, timeout=120) as response, destination.open("xb") as output:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            output.write(chunk)
            digest.update(chunk)
            size += len(chunk)
    return size, digest.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="owner/name")
    ap.add_argument("--run-id", type=int, required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--expected", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if not token:
        raise SystemExit("GH_TOKEN/GITHUB_TOKEN is required")

    artifacts = []
    page = 1
    while True:
        url = (
            f"{API}/repos/{args.repo}/actions/runs/{args.run_id}/artifacts"
            f"?per_page=100&page={page}"
        )
        payload = json.loads(request(url, token))
        batch = payload.get("artifacts", [])
        artifacts.extend(
            a for a in batch
            if not a.get("expired", False)
            and str(a.get("name", "")).startswith(args.prefix)
        )
        if len(batch) < 100:
            break
        page += 1

    by_name = {}
    for a in artifacts:
        name = str(a["name"])
        if name in by_name:
            raise SystemExit(f"duplicate artifact name: {name}")
        by_name[name] = a

    if len(by_name) != args.expected:
        raise SystemExit(
            f"artifact coverage mismatch for prefix={args.prefix!r}: "
            f"expected {args.expected}, got {len(by_name)}"
        )

    args.out.mkdir(parents=True, exist_ok=True)
    digest_verified = 0
    for n, name in enumerate(sorted(by_name), 1):
        a = by_name[name]
        expected_digest = str(a.get("digest") or "")
        if not expected_digest.startswith("sha256:"):
            raise SystemExit(f"artifact {name!r} has no SHA-256 digest; refusing unverified input")
        expected_sha256 = expected_digest.split(":", 1)[1].lower()
        if len(expected_sha256) != 64 or any(ch not in "0123456789abcdef" for ch in expected_sha256):
            raise SystemExit(f"artifact {name!r} has malformed SHA-256 digest")

        dest = safe_artifact_directory(args.out, name)
        dest.mkdir(parents=True, exist_ok=True)
        zip_path = args.out / f".{name}.download.zip"
        size, actual_sha256 = download_to_file(
            f"{API}/repos/{args.repo}/actions/artifacts/{int(a['id'])}/zip",
            token,
            zip_path,
        )
        try:
            expected_size = int(a.get("size_in_bytes", size))
            if size != expected_size:
                raise SystemExit(
                    f"artifact {name!r} size mismatch: expected {expected_size}, got {size}"
                )
            if actual_sha256 != expected_sha256:
                raise SystemExit(
                    f"artifact {name!r} SHA-256 mismatch: expected {expected_sha256}, got {actual_sha256}"
                )
            with zipfile.ZipFile(zip_path) as zf:
                safe_extract_zip(zf, dest)
            digest_verified += 1
        finally:
            zip_path.unlink(missing_ok=True)

        if n % 20 == 0 or n == args.expected:
            print(f"downloaded and digest-verified {n}/{args.expected}", flush=True)

    print(json.dumps({
        "repo": args.repo,
        "run_id": args.run_id,
        "prefix": args.prefix,
        "expected": args.expected,
        "downloaded": len(by_name),
        "digest_verified": digest_verified,
        "output": str(args.out),
    }, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
