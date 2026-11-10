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
import os
import re
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

API = "https://api.github.com"
PY_EXTS = {".py", ".pyi", ".pyx", ".pxd"}
IPYNB_EXT = ".ipynb"
MAX_WORKERS = 8


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
    with urllib.request.urlopen(req, timeout=90) as response:
        return json.loads(response.read().decode("utf-8"))


def blob_text(repo: str, sha: str, token: str):
    obj = api_json(f"{API}/repos/{repo}/git/blobs/{sha}", token)
    if obj.get("encoding") != "base64":
        raise RuntimeError(f"unexpected blob encoding for {sha}: {obj.get('encoding')}")
    return base64.b64decode(obj["content"]).decode("utf-8", errors="replace")


def dotted(node):
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        left = dotted(node.value)
        return f"{left}.{node.attr}" if left else node.attr
    if isinstance(node, ast.Call):
        return dotted(node.func)
    return ""


def record(findings, category, path, blob, node, detail):
    findings.append({
        "category": category,
        "path": path,
        "blob_sha": blob,
        "line": getattr(node, "lineno", None),
        "detail": detail,
    })


def scan_source(source: str, path: str, blob: str, findings: list, totals: Counter):
    try:
        tree = ast.parse(source, filename=path)
    except SyntaxError as exc:
        record(findings, "syntax_error", path, blob, exc, f"{exc.msg}; offset={exc.offset}")
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
                record(findings, "bare_except", path, blob, node, "bare except catches BaseException")
            elif isinstance(node.type, ast.Name) and node.type.id in {"Exception", "BaseException"}:
                record(findings, "broad_except", path, blob, node, dotted(node.type))
        elif isinstance(node, ast.Call):
            name = dotted(node.func)
            if name in {"eval", "exec"} or name.endswith(".eval") or name.endswith(".exec"):
                record(findings, "dynamic_code_execution", path, blob, node, name)
            if name == "os.system" or name.startswith("subprocess."):
                record(findings, "process_execution", path, blob, node, name)
            if name in {"pickle.load", "pickle.loads", "_pickle.load", "_pickle.loads"}:
                record(findings, "pickle_deserialization", path, blob, node, name)
            if name in {"torch.load", "torch.jit.load"}:
                record(findings, "torch_deserialization", path, blob, node, name)
            if name in {"tarfile.TarFile.extractall", "TarFile.extractall"} or name.endswith(".extractall"):
                record(findings, "archive_extraction", path, blob, node, name)
            if name in {"yaml.load", "yaml.unsafe_load"}:
                record(findings, "yaml_deserialization", path, blob, node, name)
            if name in {"urllib.request.urlopen", "requests.get", "requests.post",
                        "httpx.get", "httpx.post", "aiohttp.ClientSession"}:
                record(findings, "network_access", path, blob, node, name)
            if name in {"numpy.load", "np.load"}:
                if any(k.arg == "allow_pickle" and isinstance(k.value, ast.Constant) and k.value.value is True
                       for k in node.keywords):
                    record(findings, "numpy_pickle_load", path, blob, node, name)
            if name in {"random.seed", "numpy.random.seed", "np.random.seed"}:
                totals["explicit_random_seeds"] += 1

    for node in tree.body:
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            record(findings, "module_level_call", path, blob, node, dotted(node.value.func))
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call):
            record(findings, "module_level_assignment_call", path, blob, node, dotted(node.value.func))


def scan_notebook(source: str, path: str, blob: str, findings: list, totals: Counter):
    try:
        nb = json.loads(source)
    except Exception as exc:
        record(findings, "notebook_json_error", path, blob, None, str(exc))
        totals["notebook_json_errors"] += 1
        return
    for i, cell in enumerate(nb.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        code = "".join(cell.get("source", []))
        scan_source(code, f"{path}#cell-{i}", blob, findings, totals)


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
    commits = sorted(set(branches.values()))
    snapshots = []
    occurrences = []
    errors = []

    def get_tree(commit):
        obj = api_json(f"{API}/repos/{args.repo}/git/trees/{commit}?recursive=1", token)
        if obj.get("truncated"):
            raise RuntimeError(f"truncated tree for commit {commit}")
        entries = obj.get("tree", [])
        files = [e for e in entries if e.get("type") == "blob"
                 and Path(e["path"]).suffix.lower() in PY_EXTS | {IPYNB_EXT}]
        return commit, len(entries), files

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = [pool.submit(get_tree, c) for c in commits]
        for future in as_completed(futures):
            try:
                commit, tree_count, files = future.result()
                snapshots.append({"commit": commit, "tree_entries": tree_count, "python_files": len(files)})
                for f in files:
                    occurrences.append({
                        "commit": commit, "path": f["path"], "blob_sha": f["sha"],
                        "size": f.get("size", 0), "extension": Path(f["path"]).suffix.lower(),
                    })
            except Exception as exc:
                errors.append(str(exc))

    if errors:
        raise SystemExit("branch tree enumeration incomplete: " + " | ".join(errors[:10]))

    unique_blob_meta = {}
    for item in occurrences:
        unique_blob_meta.setdefault(item["blob_sha"], item)
    unique_blobs = list(unique_blob_meta.values())
    findings = []
    totals = Counter()
    fetched = 0

    def get_blob(item):
        return item["blob_sha"], blob_text(args.repo, item["blob_sha"], token)

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = [pool.submit(get_blob, item) for item in unique_blobs]
        for future in as_completed(futures):
            try:
                sha, source = future.result()
                item = unique_blob_meta[sha]
                fetched += 1
                if item["extension"] == IPYNB_EXT:
                    scan_notebook(source, item["path"], sha, findings, totals)
                else:
                    scan_source(source, item["path"], sha, findings, totals)
            except Exception as exc:
                errors.append(str(exc))

    if errors:
        raise SystemExit("blob audit incomplete: " + " | ".join(errors[:10]))

    path_versions = {(x["path"], x["blob_sha"]) for x in occurrences}
    unique_paths = {x["path"] for x in occurrences}
    imports = {k.removeprefix("import:"): v for k, v in totals.items() if k.startswith("import:")}
    category_counts = Counter(x["category"] for x in findings)
    result = {
        "repo": args.repo,
        "branch_count": len(branches),
        "unique_branch_head_commits": len(commits),
        "branch_snapshot_python_file_occurrences": len(occurrences),
        "unique_python_paths_across_branches": len(unique_paths),
        "unique_path_and_content_versions": len(path_versions),
        "unique_blob_contents_audited": fetched,
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
        "branch_count", "unique_branch_head_commits", "branch_snapshot_python_file_occurrences",
        "unique_python_paths_across_branches", "unique_path_and_content_versions",
        "unique_blob_contents_audited", "parsed_sources", "syntax_errors", "finding_counts",
        "elapsed_seconds",
    )}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
