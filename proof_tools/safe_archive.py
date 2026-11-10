"""Fail-closed extraction helpers for repository proof artifacts.

Only regular files and directories are accepted. Links, devices, absolute paths,
parent traversal, duplicate members, and non-empty extraction destinations are
rejected. Tar extraction additionally requires Python's data extraction filter.
"""
from __future__ import annotations

import shutil
import stat
import tarfile
import zipfile
from pathlib import Path, PurePosixPath


def _empty_root(destination: str | Path) -> Path:
    root = Path(destination).resolve()
    if root.exists():
        if not root.is_dir():
            raise ValueError(f"extraction destination is not a directory: {root}")
        if any(root.iterdir()):
            raise ValueError(f"extraction destination is not empty: {root}")
    else:
        root.mkdir(parents=True, exist_ok=True)
    return root


def _member_target(root: Path, name: str) -> Path:
    if not isinstance(name, str) or not name or "\x00" in name:
        raise ValueError(f"invalid archive member name: {name!r}")
    if "\\" in name:
        raise ValueError(f"backslash in archive member path is forbidden: {name!r}")
    posix = PurePosixPath(name)
    if posix.is_absolute() or any(part in {".."} for part in posix.parts):
        raise ValueError(f"archive member escapes extraction root: {name!r}")
    if posix.parts and len(posix.parts[0]) >= 2 and posix.parts[0][1] == ":":
        raise ValueError(f"drive-qualified archive member is forbidden: {name!r}")
    target = (root / Path(*posix.parts)).resolve()
    try:
        target.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"archive member escapes extraction root: {name!r}") from exc
    return target


def safe_extract_tar(archive: tarfile.TarFile, destination: str | Path) -> Path:
    """Extract a tar archive containing only unique regular files/directories."""
    root = _empty_root(destination)
    members = archive.getmembers()
    seen: set[Path] = set()
    for member in members:
        target = _member_target(root, member.name)
        if target in seen:
            raise ValueError(f"duplicate tar member target: {member.name!r}")
        seen.add(target)
        if member.issym() or member.islnk() or member.isdev() or member.isfifo():
            raise ValueError(f"non-regular tar member is forbidden: {member.name!r}")
        if not (member.isdir() or member.isfile()):
            raise ValueError(f"unsupported tar member type: {member.name!r}")
    try:
        archive.extractall(path=root, members=members, filter="data")
    except TypeError as exc:
        raise RuntimeError(
            "safe tar extraction requires Python tarfile filter='data'"
        ) from exc
    return root


def safe_extract_zip(archive: zipfile.ZipFile, destination: str | Path) -> Path:
    """Extract a ZIP archive without links, traversal, duplicates, or overwrites."""
    root = _empty_root(destination)
    members = archive.infolist()
    seen: set[Path] = set()
    for info in members:
        target = _member_target(root, info.filename)
        if target in seen:
            raise ValueError(f"duplicate ZIP member target: {info.filename!r}")
        seen.add(target)
        mode = (info.external_attr >> 16) & 0xFFFF
        if stat.S_ISLNK(mode) or stat.S_ISCHR(mode) or stat.S_ISBLK(mode) or stat.S_ISFIFO(mode):
            raise ValueError(f"non-regular ZIP member is forbidden: {info.filename!r}")
        if info.file_size < 0:
            raise ValueError(f"invalid ZIP member size: {info.filename!r}")

    for info in members:
        target = _member_target(root, info.filename)
        if info.is_dir():
            target.mkdir(parents=True, exist_ok=True)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        with archive.open(info, "r") as source, target.open("xb") as output:
            shutil.copyfileobj(source, output)
    return root


def safe_artifact_directory(root: str | Path, artifact_name: str) -> Path:
    """Resolve a GitHub artifact name as one safe directory component."""
    if (
        not isinstance(artifact_name, str)
        or not artifact_name
        or artifact_name in {".", ".."}
        or "/" in artifact_name
        or "\\" in artifact_name
        or "\x00" in artifact_name
    ):
        raise ValueError(f"unsafe artifact name: {artifact_name!r}")
    base = Path(root).resolve()
    base.mkdir(parents=True, exist_ok=True)
    destination = (base / artifact_name).resolve()
    try:
        destination.relative_to(base)
    except ValueError as exc:
        raise ValueError(f"artifact directory escapes output root: {artifact_name!r}") from exc
    return destination
