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
import subprocess
import time
import textwrap
from collections import Counter, defaultdict
from pathlib import Path

API = "https://api.github.com"
PY_EXTS = {".py", ".pyi", ".pyx", ".pxd"}
IPYNB_EXT = ".ipynb"
INLINE_EXTS = {".md", ".rst", ".yml", ".yaml", ".sh"}
ARCHIVE_SUFFIXES = (
    ".tar.gz.b64", ".tgz.b64", ".zip.b64",
    ".tar.gz", ".tgz", ".tar.bz2", ".tar.xz", ".tar",
    ".zip", ".whl", ".egg",
)
MAX_EMBEDDED_ARCHIVE_BYTES = 64 * 1024 * 1024
MAX_EMBEDDED_MEMBER_BYTES = 16 * 1024 * 1024
MAX_EMBEDDED_TOTAL_SOURCE_BYTES = 128 * 1024 * 1024
MAX_EMBEDDED_MEMBERS = 20000


def git_output(arguments, input_bytes=None) -> bytes:
    result = subprocess.run(
        ["git", *arguments],
        input=input_bytes,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode:
        message = result.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"git {' '.join(arguments)} failed: {message}")
    return result.stdout


def branch_heads():
    raw = git_output([
        "for-each-ref", "--format=%(refname)%09%(objectname)",
        "refs/remotes/origin/",
    ])
    branches = {}
    for line in raw.decode("utf-8", errors="replace").splitlines():
        if not line:
            continue
        ref, sha = line.split("\t", 1)
        if ref.endswith("/HEAD"):
            continue
        name = ref.removeprefix("refs/remotes/origin/")
        branches[name] = sha
    return branches


def git_tree(commit: str):
    raw = git_output(["ls-tree", "-r", "-t", "-z", commit])
    entries = []
    for record_bytes in raw.split(b"\0"):
        if not record_bytes:
            continue
        metadata, path_bytes = record_bytes.split(b"\t", 1)
        mode, kind, sha = metadata.decode("ascii").split()
        path = path_bytes.decode("utf-8", errors="replace")
        entries.append({"mode": mode, "type": kind, "sha": sha, "path": path})
    return entries


def git_blob_batch(shas: list[str]) -> dict[str, bytes]:
    if not shas:
        return {}
    request = "".join(f"{sha}\n" for sha in shas).encode("ascii")
    raw = git_output(["cat-file", "--batch"], request)
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
        git_header = b"blob " + str(size).encode("ascii") + b"\0"
        actual_sha = hashlib.sha1(git_header + payload).hexdigest()
        if actual_sha != sha:
            raise RuntimeError(f"git cat-file SHA mismatch for {sha}: got {actual_sha}")
        blobs[sha] = payload
        offset = finish + 1
    if offset != len(raw):
        raise RuntimeError("unexpected trailing bytes from git cat-file --batch")
    return blobs


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
            if not (path.startswith("tests/") or "/tests/" in path or path.endswith("_test.py")):
                try:
                    condition = ast.unparse(node.test)
                except Exception:
                    condition = "assert condition"
                record(
                    findings, "assert_statement", path, blob, node,
                    f"optimization-sensitive guard: {condition}", locations
                )
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
                shell_true = any(
                    keyword.arg == "shell"
                    and isinstance(keyword.value, ast.Constant)
                    and keyword.value.value is True
                    for keyword in node.keywords
                )
                if name in {"os.system", "os.popen", "subprocess.getoutput", "subprocess.getstatusoutput"} or shell_true:
                    record(findings, "shell_execution", path, blob, node, name, locations)
            if name in {"pickle.load", "pickle.loads", "_pickle.load", "_pickle.loads"}:
                record(findings, "pickle_deserialization", path, blob, node, name, locations)
            if name in {"torch.load", "torch.jit.load"}:
                record(findings, "torch_deserialization", path, blob, node, name, locations)
            if name in {"tarfile.TarFile.extractall", "TarFile.extractall"} or name.endswith(".extractall"):
                category = (
                    "guarded_archive_extraction"
                    if path.replace("\\", "/").endswith("proof_tools/safe_archive.py")
                    else "archive_extraction"
                )
                record(findings, category, path, blob, node, name, locations)
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



def iter_inline_python_blocks(source: str, path: str):
    """Yield (kind, line, code) for explicitly marked Python snippets.

    Markdown/RST Python fences and Python heredocs in YAML/shell files are
    scanned. Arbitrary prose and shell fragments are deliberately not guessed
    to be Python.
    """
    suffix = Path(path).suffix.lower()
    lines = source.splitlines()

    if suffix in {".md", ".rst"}:
        i = 0
        fence_open = re.compile(r"""^[ \t]{0,3}([\x60]{3,}|~{3,})[ \t]*([A-Za-z0-9_+.-]+)?(?:[ \t].*)?$""")
        while i < len(lines):
            match = fence_open.match(lines[i])
            if not match:
                i += 1
                continue
            fence = match.group(1)
            language = (match.group(2) or "").lower()
            if language not in {"python", "python3", "py", "ipython"}:
                i += 1
                continue
            marker = fence[0]
            minimum = len(fence)
            close = re.compile(
                r"^[ \t]{0,3}" + re.escape(marker) + "{" + str(minimum) + r",}[ \t]*$"
            )
            start = i + 1
            j = start
            while j < len(lines) and not close.match(lines[j]):
                j += 1
            if j < len(lines):
                yield (f"fence-{language}", i + 1, textwrap.dedent("\n".join(lines[start:j])))
                i = j + 1
            else:
                i += 1

    if suffix in {".yml", ".yaml", ".sh"}:
        heredoc = re.compile(
            r"""^[ \t]*python(?:3)?(?:[ \t]+-[A-Za-z]+)*[ \t]+-[ \t]*<<['"]?([A-Za-z_][A-Za-z0-9_]*)['"]?[ \t]*$"""
        )
        i = 0
        while i < len(lines):
            match = heredoc.match(lines[i])
            if not match:
                i += 1
                continue
            delimiter = match.group(1)
            start = i + 1
            j = start
            while j < len(lines) and lines[j].strip() != delimiter:
                j += 1
            if j < len(lines):
                yield (f"heredoc-{delimiter}", i + 1, textwrap.dedent("\n".join(lines[start:j])))
                i = j + 1
            else:
                i += 1


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

    t0 = time.time()
    branches = branch_heads()
    if "main" not in branches:
        raise SystemExit("local shallow fetch does not contain refs/remotes/origin/main")
    main_head = branches["main"]
    commits = sorted(set(branches.values()))
    snapshots = []
    occurrences = []
    archive_occurrences = []
    inline_occurrences = []

    for commit in commits:
        entries = git_tree(commit)
        files = [
            entry for entry in entries
            if entry["type"] == "blob"
            and Path(entry["path"]).suffix.lower() in PY_EXTS | {IPYNB_EXT}
        ]
        archives = [
            entry for entry in entries
            if entry["type"] == "blob"
            and entry["path"].lower().endswith(ARCHIVE_SUFFIXES)
        ]
        inline_files = [
            entry for entry in entries
            if entry["type"] == "blob"
            and Path(entry["path"]).suffix.lower() in INLINE_EXTS
        ]
        snapshots.append({
            "commit": commit,
            "tree_entries": len(entries),
            "python_files": len(files),
            "embedded_archives": len(archives),
            "inline_source_files": len(inline_files),
        })
        for item in files:
            occurrences.append({
                "commit": commit,
                "path": item["path"],
                "blob_sha": item["sha"],
                "extension": Path(item["path"]).suffix.lower(),
            })
        for item in archives:
            archive_occurrences.append({
                "commit": commit,
                "path": item["path"],
                "blob_sha": item["sha"],
            })
        for item in inline_files:
            inline_occurrences.append({
                "commit": commit,
                "path": item["path"],
                "blob_sha": item["sha"],
                "extension": Path(item["path"]).suffix.lower(),
            })

    blob_locations = defaultdict(list)
    for item in occurrences:
        blob_locations[item["blob_sha"]].append({
            "commit": item["commit"],
            "path": item["path"],
            "extension": item["extension"],
        })
    for locations in blob_locations.values():
        locations.sort(key=lambda x: (x["path"], x["commit"], x["extension"]))

    archive_locations = defaultdict(list)
    for item in archive_occurrences:
        archive_locations[item["blob_sha"]].append({
            "commit": item["commit"], "path": item["path"]
        })
    for locations in archive_locations.values():
        locations.sort(key=lambda x: (x["path"], x["commit"]))

    inline_locations = defaultdict(list)
    for item in inline_occurrences:
        inline_locations[item["blob_sha"]].append({
            "commit": item["commit"],
            "path": item["path"],
            "extension": item["extension"],
        })
    for locations in inline_locations.values():
        locations.sort(key=lambda x: (x["path"], x["commit"], x["extension"]))

    all_blob_shas = sorted(set(blob_locations) | set(archive_locations) | set(inline_locations))
    blob_contents = git_blob_batch(all_blob_shas)

    findings = []
    totals = Counter()
    for sha in sorted(blob_locations):
        locations = blob_locations[sha]
        representative = locations[0]["path"]
        extensions = {x["extension"] for x in locations}
        source = blob_contents[sha].decode("utf-8", errors="replace")
        if IPYNB_EXT in extensions:
            scan_notebook(source, representative, sha, findings, totals, locations)
        if extensions & PY_EXTS:
            scan_source(source, representative, sha, findings, totals, locations)

    for archive_sha in sorted(archive_locations):
        locations = archive_locations[archive_sha]
        representative = locations[0]["path"]
        scan_embedded_archive(
            blob_contents[archive_sha], representative, archive_sha,
            findings, totals, locations
        )

    for inline_sha in sorted(inline_locations):
        locations = inline_locations[inline_sha]
        representative = locations[0]["path"]
        source = blob_contents[inline_sha].decode("utf-8", errors="replace")
        for index, (kind, line, code) in enumerate(
            iter_inline_python_blocks(source, representative), 1
        ):
            block_path = f"{representative}#{kind}-{index}"
            scan_source(code, block_path, inline_sha, findings, totals, locations)
            totals["inline_python_blocks"] += 1

    path_versions = {(x["path"], x["blob_sha"]) for x in occurrences}
    unique_paths = {x["path"] for x in occurrences}
    main_locations = [x for x in occurrences if x["commit"] == main_head]
    main_blob_shas = {x["blob_sha"] for x in main_locations}
    main_embedded_archive_locations = [
        x for x in archive_occurrences if x["commit"] == main_head
    ]
    main_embedded_syntax_errors = sum(
        1 for item in findings
        if item["category"] == "syntax_error"
        and item["blob_sha"].startswith("sha256:")
        and any(loc["commit"] == main_head for loc in item.get("locations", []))
    )
    main_syntax_errors = sum(
        1 for item in findings
        if item["category"] == "syntax_error" and item["blob_sha"] in main_blob_shas
    ) + main_embedded_syntax_errors

    imports = {
        k.removeprefix("import:"): v
        for k, v in totals.items() if k.startswith("import:")
    }
    category_counts = Counter(x["category"] for x in findings)
    result = {
        "repo": args.repo,
        "branch_count": len(branches),
        "branch_heads": dict(sorted(branches.items())),
        "main_head_sha": main_head,
        "main_python_file_occurrences": len(main_locations),
        "main_unique_python_blobs": len(main_blob_shas),
        "main_embedded_archive_file_occurrences": len(main_embedded_archive_locations),
        "main_inline_source_file_occurrences": sum(
            1 for x in inline_occurrences if x["commit"] == main_head
        ),
        "branch_snapshot_inline_source_file_occurrences": len(inline_occurrences),
        "unique_inline_source_blobs_audited": len(inline_locations),
        "inline_python_blocks": totals["inline_python_blocks"],
        "main_syntax_errors": main_syntax_errors,
        "main_ast_status": "PASS" if main_syntax_errors == 0 else "FAIL",
        "unique_branch_head_commits": len(commits),
        "branch_snapshot_python_file_occurrences": len(occurrences),
        "unique_python_paths_across_branches": len(unique_paths),
        "unique_path_and_content_versions": len(path_versions),
        "unique_blob_contents_audited": len(blob_locations),
        "unique_git_blobs_loaded": len(blob_contents),
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
        "scope_note": (
            "Static inspection only; repository code was not executed. Every distinct "
            "tracked Python blob, every Python source member in supported tar/ZIP "
            "archives, and explicitly marked Python fence/heredoc block at branch heads "
            "was parsed. Nested archives, unmarked inline fragments, unsupported archive "
            "formats, tags, and non-head historical commits are outside this scan."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in (
        "branch_count", "unique_branch_head_commits", "main_head_sha",
        "main_python_file_occurrences", "main_embedded_archive_file_occurrences",
        "main_syntax_errors", "main_ast_status",
        "branch_snapshot_python_file_occurrences", "unique_python_paths_across_branches",
        "unique_path_and_content_versions", "unique_blob_contents_audited",
        "unique_git_blobs_loaded", "parsed_sources", "syntax_errors",
        "embedded_archive_file_occurrences", "unique_embedded_archive_blobs",
        "embedded_archives_scanned", "embedded_archive_members",
        "embedded_python_sources", "embedded_archive_unsafe_members",
        "embedded_archive_errors", "finding_counts", "elapsed_seconds",
    )}, indent=2))
    return 1 if main_syntax_errors or totals["embedded_archive_errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
