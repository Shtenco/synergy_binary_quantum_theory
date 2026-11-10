#!/usr/bin/env python3
"""Branch-wide static audit of every Python source version in a GitHub repository.

Enumerates all branch heads, deduplicates commit snapshots and blob contents,
parses each distinct Python source blob, and records risky API usage. It does
not execute repository code.
"""
from __future__ import annotations

import argparse
import ast
import base64
import json
import hashlib
import io
import stat
import tarfile
import zipfile
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

API = "https://api.github.com"
PY_EXTS = {".py", ".pyi", ".pyx", ".pxd"}
IPYNB_EXT = ".ipynb"
MAX_WORKERS = 8
ARCHIVE_SUFFIXES = (
    ".tar.gz.b64", ".tgz.b64", ".zip.b64",
    ".tar.gz", ".tgz", ".tar.bz2", ".tar.xz", ".tar",
    ".zip", ".whl", ".egg",
)
MAX_EMBEDDED_ARCHIVE_BYTES = 64 * 1024 * 1024
MAX_EMBEDDED_MEMBER_BYTES = 16 * 1024 * 1024
MAX_EMBEDDED_TOTAL_SOURCE_BYTES = 128 * 1024 * 1024
MAX_EMBEDDED_MEMBERS = 20000


def api_json(url: str, token: str):
    req = urllib.request.Request(
        url,
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "BQG-branch-wide-python-audit",
        },
    )
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, timeout=90) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            if exc.code not in {403, 429, 500, 502, 503, 504} or attempt == 5:
                raise
            retry_after = exc.headers.get("Retry-After")
            remaining = exc.headers.get("X-RateLimit-Remaining")
            reset = exc.headers.get("X-RateLimit-Reset")
            if retry_after:
                delay = max(1.0, float(retry_after))
            elif remaining == "0" and reset:
                delay = max(1.0, float(reset) - time.time() + 1.0)
                if delay > 120:
                    raise RuntimeError(
                        f"GitHub API rate limit exhausted; reset in {delay:.0f}s"
                    ) from exc
            else:
                delay = min(30.0, 2.0 ** attempt * 2.0)
            time.sleep(delay)


def blob_bytes(repo: str, sha: str, commit: str, path: str) -> bytes:
    # Raw content is not subject to the authenticated REST API quota. Verify
    # the Git object SHA-1 so fetched bytes remain pinned to the tree entry.
    quoted_path = urllib.parse.quote(path, safe="/")
    url = f"https://raw.githubusercontent.com/{repo}/{commit}/{quoted_path}"
    req = urllib.request.Request(url, headers={"User-Agent": "BQG-branch-wide-python-audit"})
    with urllib.request.urlopen(req, timeout=90) as response:
        raw = response.read()
    header = b"blob " + str(len(raw)).encode("ascii") + b"\0"
    actual_sha = hashlib.sha1(header + raw).hexdigest()
    if actual_sha != sha:
        raise RuntimeError(
            f"raw blob SHA mismatch for {commit}:{path}: expected {sha}, got {actual_sha}"
        )
    return raw


def blob_text(repo: str, sha: str, commit: str, path: str) -> str:
    return blob_bytes(repo, sha, commit, path).decode("utf-8", errors="replace")


def dotted(node):
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        left = dotted(node.value)
        return f"{left}.{node.attr}" if left else node.attr
    if isinstance(node, ast.Call):
        return dotted(node.func)
    return ""


def record(findings, category, path, blob, node, detail, locations=None):
    item = {
        "category": category,
        "path": path,
        "blob_sha": blob,
        "line": getattr(node, "lineno", None),
        "detail": detail,
    }
    if locations:
        item["locations"] = locations
    findings.append(item)


def scan_source(source: str, path: str, blob: str, findings: list, totals: Counter, locations=None):
    try:
        tree = ast.parse(source, filename=path)
    except SyntaxError as exc:
        record(findings, "syntax_error", path, blob, exc, f"{exc.msg}; offset={exc.offset}", locations)
        totals["syntax_errors"] += 1
        return

    totals["parsed_sources"] += 1
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                totals[f"import:{alias.name.split('.')[0]}"] += 1
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                totals[f"import:{node.module.split('.')[0]}"] += 1
        elif isinstance(node, ast.Assert):
            totals["assert_statements"] += 1
        elif isinstance(node, ast.ExceptHandler):
            if node.type is None:
                record(findings, "bare_except", path, blob, node, "bare except catches BaseException", locations)
            elif isinstance(node.type, ast.Name) and node.type.id in {"Exception", "BaseException"}:
                record(findings, "broad_except", path, blob, node, dotted(node.type), locations)
        elif isinstance(node, ast.Call):
            name = dotted(node.func)
            if name in {"eval", "exec"} or name in {"builtins.eval", "builtins.exec"}:
                record(findings, "dynamic_code_execution", path, blob, node, name, locations)
            if name == "os.system" or name.startswith("subprocess."):
                record(findings, "process_execution", path, blob, node, name, locations)
            if name in {"pickle.load", "pickle.loads", "_pickle.load", "_pickle.loads"}:
                record(findings, "pickle_deserialization", path, blob, node, name, locations)
            if name in {"torch.load", "torch.jit.load"}:
                record(findings, "torch_deserialization", path, blob, node, name, locations)
            if name in {"tarfile.TarFile.extractall", "TarFile.extractall"} or name.endswith(".extractall"):
                record(findings, "archive_extraction", path, blob, node, name, locations)
            if name in {"yaml.load", "yaml.unsafe_load"}:
                record(findings, "yaml_deserialization", path, blob, node, name, locations)
            if name in {"urllib.request.urlopen", "requests.get", "requests.post",
                        "httpx.get", "httpx.post", "aiohttp.ClientSession"}:
                record(findings, "network_access", path, blob, node, name, locations)
            if name in {"numpy.load", "np.load"}:
                if any(k.arg == "allow_pickle" and isinstance(k.value, ast.Constant) and k.value.value is True
                       for k in node.keywords):
                    record(findings, "numpy_pickle_load", path, blob, node, name, locations)
            if name in {"random.seed", "numpy.random.seed", "np.random.seed"}:
                totals["explicit_random_seeds"] += 1

    for node in tree.body:
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            record(findings, "module_level_call", path, blob, node, dotted(node.value.func), locations)
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call):
            record(findings, "module_level_assignment_call", path, blob, node, dotted(node.value.func), locations)


def scan_notebook(source: str, path: str, blob: str, findings: list, totals: Counter, locations=None):
    try:
        nb = json.loads(source)
    except Exception as exc:
        record(findings, "notebook_json_error", path, blob, None, str(exc), locations)
        totals["notebook_json_errors"] += 1
        return
    for i, cell in enumerate(nb.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        code = "".join(cell.get("source", []))
        scan_source(code, f"{path}#cell-{i}", blob, findings, totals, locations)



def safe_archive_member(name: str) -> bool:
    normalized = name.replace("\\", "/")
    parts = Path(normalized).parts
    return (
        bool(normalized)
        and "\x00" not in normalized
        and not normalized.startswith("/")
        and not re.match(r"^[A-Za-z]:", normalized)
        and ".." not in parts
    )


def scan_embedded_source(payload: bytes, member_path: str, findings: list,
                         totals: Counter, locations=None):
    source_sha = hashlib.sha256(payload).hexdigest()
    try:
        source = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        record(findings, "embedded_source_decode_error", member_path,
               f"sha256:{source_sha}", None, str(exc), locations)
        totals["embedded_source_decode_errors"] += 1
        return
    if Path(member_path.split("!/", 1)[-1]).suffix.lower() == IPYNB_EXT:
        scan_notebook(source, member_path, f"sha256:{source_sha}",
                      findings, totals, locations)
    else:
        scan_source(source, member_path, f"sha256:{source_sha}",
                    findings, totals, locations)
    totals["embedded_python_sources"] += 1
    totals["embedded_python_source_bytes"] += len(payload)


def scan_embedded_archive(payload: bytes, archive_path: str, archive_sha: str,
                          findings: list, totals: Counter, locations=None):
    suffix = archive_path.lower()
    try:
        if suffix.endswith(".b64"):
            compact = b"".join(payload.split()).decode("ascii")
            raw = base64.b64decode(compact, validate=True)
        else:
            raw = payload
    except Exception as exc:
        record(findings, "embedded_archive_decode_error", archive_path,
               archive_sha, None, str(exc), locations)
        totals["embedded_archive_errors"] += 1
        return
    if len(raw) > MAX_EMBEDDED_ARCHIVE_BYTES:
        record(findings, "embedded_archive_size_limit", archive_path,
               archive_sha, None, f"archive bytes {len(raw)} exceed limit", locations)
        totals["embedded_archive_errors"] += 1
        return

    suffix = archive_path.lower()
    member_count = 0
    source_bytes = 0
    try:
        if suffix.endswith((".tar.gz.b64", ".tgz.b64", ".tar.gz", ".tgz",
                            ".tar.bz2", ".tar.xz", ".tar")):
            with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as archive:
                for member in archive:
                    member_count += 1
                    if member_count > MAX_EMBEDDED_MEMBERS:
                        raise RuntimeError("embedded archive member-count limit exceeded")
                    if member.issym() or member.islnk() or member.isdev() or member.isfifo():
                        record(findings, "embedded_archive_link_or_device",
                               f"{archive_path}!/{member.name}", archive_sha,
                               None, "link/device member is not followed", locations)
                        totals["embedded_archive_unsafe_members"] += 1
                        continue
                    if not member.isfile():
                        continue
                    if not safe_archive_member(member.name):
                        record(findings, "embedded_archive_unsafe_path",
                               f"{archive_path}!/{member.name}", archive_sha,
                               None, "absolute or traversal member path", locations)
                        totals["embedded_archive_unsafe_members"] += 1
                        continue
                    if Path(member.name).suffix.lower() not in PY_EXTS | {IPYNB_EXT}:
                        continue
                    if member.size > MAX_EMBEDDED_MEMBER_BYTES:
                        record(findings, "embedded_archive_member_size_limit",
                               f"{archive_path}!/{member.name}", archive_sha,
                               None, f"member bytes {member.size} exceed limit", locations)
                        totals["embedded_archive_errors"] += 1
                        continue
                    stream = archive.extractfile(member)
                    if stream is None:
                        raise RuntimeError(f"could not read regular member {member.name}")
                    payload = stream.read(MAX_EMBEDDED_MEMBER_BYTES + 1)
                    if len(payload) != member.size or len(payload) > MAX_EMBEDDED_MEMBER_BYTES:
                        raise RuntimeError(f"member size mismatch or limit exceeded: {member.name}")
                    source_bytes += len(payload)
                    if source_bytes > MAX_EMBEDDED_TOTAL_SOURCE_BYTES:
                        raise RuntimeError("embedded archive total source-byte limit exceeded")
                    scan_embedded_source(payload, f"{archive_path}!/{member.name}",
                                         findings, totals, locations)
        elif suffix.endswith((".zip.b64", ".zip", ".whl", ".egg")):
            with zipfile.ZipFile(io.BytesIO(raw)) as archive:
                for member in archive.infolist():
                    member_count += 1
                    if member_count > MAX_EMBEDDED_MEMBERS:
                        raise RuntimeError("embedded archive member-count limit exceeded")
                    mode = member.external_attr >> 16
                    if (stat.S_ISLNK(mode) or stat.S_ISCHR(mode)
                            or stat.S_ISBLK(mode) or stat.S_ISFIFO(mode)):
                        record(findings, "embedded_archive_link_or_device",
                               f"{archive_path}!/{member.filename}", archive_sha,
                               None, "symlink member is not followed", locations)
                        totals["embedded_archive_unsafe_members"] += 1
                        continue
                    if member.is_dir():
                        continue
                    if not safe_archive_member(member.filename):
                        record(findings, "embedded_archive_unsafe_path",
                               f"{archive_path}!/{member.filename}", archive_sha,
                               None, "absolute or traversal member path", locations)
                        totals["embedded_archive_unsafe_members"] += 1
                        continue
                    if Path(member.filename).suffix.lower() not in PY_EXTS | {IPYNB_EXT}:
                        continue
                    if member.file_size > MAX_EMBEDDED_MEMBER_BYTES:
                        record(findings, "embedded_archive_member_size_limit",
                               f"{archive_path}!/{member.filename}", archive_sha,
                               None, f"member bytes {member.file_size} exceed limit", locations)
                        totals["embedded_archive_errors"] += 1
                        continue
                    payload = archive.read(member)
                    if len(payload) != member.file_size:
                        raise RuntimeError(f"member size mismatch: {member.filename}")
                    source_bytes += len(payload)
                    if source_bytes > MAX_EMBEDDED_TOTAL_SOURCE_BYTES:
                        raise RuntimeError("embedded archive total source-byte limit exceeded")
                    scan_embedded_source(payload, f"{archive_path}!/{member.filename}",
                                         findings, totals, locations)
        else:
            raise RuntimeError("unsupported embedded archive suffix")
    except Exception as exc:
        record(findings, "embedded_archive_error", archive_path,
               archive_sha, None, f"{type(exc).__name__}: {exc}", locations)
        totals["embedded_archive_errors"] += 1
        return

    totals["embedded_archives_scanned"] += 1
    totals["embedded_archive_members"] += member_count
    totals["embedded_archive_source_bytes"] += source_bytes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY", "Shtenco/synergy_binary_quantum_theory"))
    ap.add_argument("--out", type=Path, default=Path("python-branch-audit.json"))
    args = ap.parse_args()
    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if not token:
        raise SystemExit("GITHUB_TOKEN or GH_TOKEN is required")

    t0 = time.time()
    refs = api_json(f"{API}/repos/{args.repo}/git/matching-refs/heads/", token)
    branches = {x["ref"].removeprefix("refs/heads/"): x["object"]["sha"] for x in refs}
    if "main" not in branches:
        raise SystemExit("branch-wide audit did not receive a main branch head")
    main_head = branches["main"]
    commits = sorted(set(branches.values()))
    snapshots = []
    occurrences = []
    archive_occurrences = []
    errors = []

    def get_tree(commit):
        obj = api_json(f"{API}/repos/{args.repo}/git/trees/{commit}?recursive=1", token)
        if obj.get("truncated"):
            raise RuntimeError(f"truncated tree for commit {commit}")
        entries = obj.get("tree", [])
        files = [e for e in entries if e.get("type") == "blob"
                 and Path(e["path"]).suffix.lower() in PY_EXTS | {IPYNB_EXT}]
        archives = [e for e in entries if e.get("type") == "blob"
                    and e["path"].lower().endswith(ARCHIVE_SUFFIXES)]
        return commit, len(entries), files, archives

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = [pool.submit(get_tree, c) for c in commits]
        for future in as_completed(futures):
            try:
                commit, tree_count, files, archives = future.result()
                snapshots.append({"commit": commit, "tree_entries": tree_count,
                                  "python_files": len(files), "embedded_archives": len(archives)})
                for f in files:
                    occurrences.append({
                        "commit": commit, "path": f["path"], "blob_sha": f["sha"],
                        "size": f.get("size", 0), "extension": Path(f["path"]).suffix.lower(),
                    })
                for f in archives:
                    archive_occurrences.append({
                        "commit": commit, "path": f["path"], "blob_sha": f["sha"],
                        "size": f.get("size", 0),
                    })
            except Exception as exc:
                errors.append(str(exc))

    if errors:
        raise SystemExit("branch tree enumeration incomplete: " + " | ".join(errors[:10]))

    blob_locations = defaultdict(list)
    for item in occurrences:
        blob_locations[item["blob_sha"]].append({
            "commit": item["commit"],
            "path": item["path"],
            "extension": item["extension"],
        })
    for locations in blob_locations.values():
        locations.sort(key=lambda x: (x["path"], x["commit"], x["extension"]))
    unique_blobs = sorted(blob_locations)
    findings = []
    totals = Counter()
    fetched = 0

    def get_blob(sha):
        location = blob_locations[sha][0]
        return sha, blob_text(
            args.repo, sha, location["commit"], location["path"]
        )

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = [pool.submit(get_blob, sha) for sha in unique_blobs]
        for future in as_completed(futures):
            try:
                sha, source = future.result()
                locations = blob_locations[sha]
                representative = locations[0]["path"]
                extensions = {x["extension"] for x in locations}
                fetched += 1
                if IPYNB_EXT in extensions:
                    scan_notebook(source, representative, sha, findings, totals, locations)
                if extensions & PY_EXTS:
                    scan_source(source, representative, sha, findings, totals, locations)
            except Exception as exc:
                errors.append(str(exc))

    if errors:
        raise SystemExit("blob audit incomplete: " + " | ".join(errors[:10]))

    archive_locations = defaultdict(list)
    for item in archive_occurrences:
        archive_locations[item["blob_sha"]].append({
            "commit": item["commit"], "path": item["path"]
        })
    for locations in archive_locations.values():
        locations.sort(key=lambda x: (x["path"], x["commit"]))
    for archive_sha in sorted(archive_locations):
        locations = archive_locations[archive_sha]
        representative = locations[0]["path"]
        try:
            location = locations[0]
            payload = blob_bytes(
                args.repo, archive_sha, location["commit"], location["path"]
            )
            scan_embedded_archive(payload, representative, archive_sha,
                                  findings, totals, locations)
        except Exception as exc:
            record(findings, "embedded_archive_fetch_error", representative,
                   archive_sha, None, f"{type(exc).__name__}: {exc}", locations)
            totals["embedded_archive_errors"] += 1

    path_versions = {(x["path"], x["blob_sha"]) for x in occurrences}
    unique_paths = {x["path"] for x in occurrences}
    main_locations = [x for x in occurrences if x["commit"] == main_head]
    main_blob_shas = {x["blob_sha"] for x in main_locations}
    main_syntax_errors = sum(
        1 for item in findings
        if item["category"] == "syntax_error" and item["blob_sha"] in main_blob_shas
    )
    imports = {k.removeprefix("import:"): v for k, v in totals.items() if k.startswith("import:")}
    category_counts = Counter(x["category"] for x in findings)
    result = {
        "repo": args.repo,
        "branch_count": len(branches),
        "branch_heads": dict(sorted(branches.items())),
        "main_head_sha": main_head,
        "main_python_file_occurrences": len(main_locations),
        "main_unique_python_blobs": len(main_blob_shas),
        "main_syntax_errors": main_syntax_errors,
        "main_ast_status": "PASS" if main_syntax_errors == 0 else "FAIL",
        "unique_branch_head_commits": len(commits),
        "branch_snapshot_python_file_occurrences": len(occurrences),
        "unique_python_paths_across_branches": len(unique_paths),
        "unique_path_and_content_versions": len(path_versions),
        "unique_blob_contents_audited": fetched,
        "embedded_archive_file_occurrences": len(archive_occurrences),
        "unique_embedded_archive_blobs": len(archive_locations),
        "embedded_archives_scanned": totals["embedded_archives_scanned"],
        "embedded_archive_members": totals["embedded_archive_members"],
        "embedded_python_sources": totals["embedded_python_sources"],
        "embedded_python_source_bytes": totals["embedded_python_source_bytes"],
        "embedded_archive_unsafe_members": totals["embedded_archive_unsafe_members"],
        "embedded_archive_errors": totals["embedded_archive_errors"],
        "tree_snapshots": len(snapshots),
        "tree_entries_across_unique_heads": sum(x["tree_entries"] for x in snapshots),
        "parsed_sources": totals["parsed_sources"],
        "syntax_errors": totals["syntax_errors"],
        "assert_statements": totals["assert_statements"],
        "explicit_random_seeds": totals["explicit_random_seeds"],
        "imports": dict(sorted(imports.items())),
        "finding_counts": dict(sorted(category_counts.items())),
        "findings": sorted(findings, key=lambda x: (x["category"], x["path"], x["line"] or 0)),
        "snapshots": sorted(snapshots, key=lambda x: x["commit"]),
        "elapsed_seconds": round(time.time() - t0, 2),
        "scope_note": "Static inspection only; repository code was not executed. Every distinct Python blob reachable from a branch head was parsed once.",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in (
        "branch_count", "unique_branch_head_commits", "main_head_sha",
        "main_python_file_occurrences", "main_syntax_errors", "main_ast_status",
        "branch_snapshot_python_file_occurrences", "unique_python_paths_across_branches",
        "unique_path_and_content_versions", "unique_blob_contents_audited",
        "parsed_sources", "syntax_errors", "embedded_archive_file_occurrences",
        "unique_embedded_archive_blobs", "embedded_archives_scanned",
        "embedded_archive_members", "embedded_python_sources",
        "embedded_archive_unsafe_members", "embedded_archive_errors",
        "finding_counts", "elapsed_seconds",
    )}, indent=2))
    return 1 if main_syntax_errors or totals['embedded_archive_errors'] else 0


if __name__ == "__main__":
    raise SystemExit(main())
