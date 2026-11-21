#!/usr/bin/env python3
"""Validate GitHub Markdown/LaTeX boundaries.

By default every repository Markdown file is scanned. Optional positional paths
allow CI to enforce a strict public-surface gate while historical research notes
are cleaned separately.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEX_COMMAND = re.compile(
    r"\\(?:begin|end|frac|sqrt|sum|prod|int|mathcal|mathfrak|mathrm|text|"
    r"operatorname|alpha|beta|gamma|delta|lambda|mu|pi|rho|sigma|theta|"
    r"Delta|Lambda|Phi|Psi|Theta|Omega|cdot|times|approx|leq|geq|infty)\b"
)
ENVIRONMENT = re.compile(r"\\(begin|end)\{([^{}]+)\}")


@dataclass(frozen=True)
class Problem:
    path: Path
    line: int
    message: str

    def __str__(self) -> str:
        try:
            display = self.path.relative_to(ROOT)
        except ValueError:
            display = self.path
        return f"{display}:{self.line}: {self.message}"


def validate_expression(source: str, path: Path, line: int) -> list[Problem]:
    problems: list[Problem] = []
    depth = 0
    escaped = False
    for char in source:
        if escaped:
            escaped = False
        elif char == "\\":
            escaped = True
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth < 0:
                problems.append(Problem(path, line, "unmatched '}' in math"))
                depth = 0
    if depth:
        problems.append(Problem(path, line, f"{depth} unclosed '{{' group(s) in math"))

    stack: list[str] = []
    for kind, name in ENVIRONMENT.findall(source):
        if kind == "begin":
            stack.append(name)
        elif not stack or stack.pop() != name:
            problems.append(Problem(path, line, f"unmatched \\end{{{name}}}"))
    for name in stack:
        problems.append(Problem(path, line, f"unclosed \\begin{{{name}}}"))
    return problems


def validate_markdown(path: Path) -> tuple[int, list[Problem]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    problems: list[Problem] = []
    expressions = 0
    in_fence = False
    display_source: list[str] | None = None
    display_line = 0

    for number, line in enumerate(lines, 1):
        # Normalize the standard LaTeX delimiters supported by GitHub.
        # Without this, valid \\( ... \\) and \\[ ... \\] math is
        # incorrectly scanned as plain Markdown prose.
        line = line.replace(r"\[", lambda _: "$" + "$").replace(r"\]", lambda _: "$" + "$")
        line = line.replace(r"\(", "$").replace(r"\)", "$")
        if re.match(r"^\s*(```|~~~)", line):
            in_fence = not in_fence
            continue
        if in_fence:
            continue

        outside: list[str] = []
        position = 0
        while position < len(line):
            if display_source is not None:
                closing = line.find("$$", position)
                if closing < 0:
                    display_source.append(line[position:])
                    position = len(line)
                    continue
                display_source.append(line[position:closing])
                source = "\n".join(display_source)
                problems.extend(validate_expression(source, path, display_line))
                expressions += 1
                display_source = None
                position = closing + 2
                continue

            opening = line.find("$", position)
            if opening < 0:
                outside.append(line[position:])
                break
            if opening and line[opening - 1] == "\\":
                outside.append(line[position : opening + 1])
                position = opening + 1
                continue
            outside.append(line[position:opening])
            if line.startswith("$$", opening):
                display_source = []
                display_line = number
                position = opening + 2
                continue

            closing = opening + 1
            while True:
                closing = line.find("$", closing)
                if closing < 0:
                    problems.append(Problem(path, number, "unclosed inline '$' delimiter"))
                    position = len(line)
                    break
                if line[closing - 1] != "\\":
                    source = line[opening + 1 : closing]
                    problems.extend(validate_expression(source, path, number))
                    if line.startswith("|") and "|" in source:
                        problems.append(Problem(path, number, "unescaped '|' inside table math"))
                    expressions += 1
                    position = closing + 1
                    break
                closing += 1

        plain = re.sub(r"`[^`]*`", "", " ".join(outside))
        if TEX_COMMAND.search(plain):
            problems.append(Problem(path, number, "TeX command outside a math delimiter"))

    if in_fence:
        problems.append(Problem(path, len(lines), "unclosed fenced code block"))
    if display_source is not None:
        problems.append(Problem(path, display_line, "unclosed '$$' delimiter"))
    return expressions, problems


def resolve_paths(values: list[str]) -> tuple[list[Path], list[str]]:
    errors: list[str] = []
    if not values:
        return (
            sorted(path for path in ROOT.rglob("*.md") if ".git" not in path.parts),
            errors,
        )

    paths: list[Path] = []
    for value in values:
        candidate = (ROOT / value).resolve()
        try:
            candidate.relative_to(ROOT)
        except ValueError:
            errors.append(f"unsafe path outside repository: {value}")
            continue
        if not candidate.is_file():
            errors.append(f"Markdown file not found: {value}")
            continue
        if candidate.suffix.lower() != ".md":
            errors.append(f"not a Markdown file: {value}")
            continue
        paths.append(candidate)
    return sorted(set(paths)), errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "paths",
        nargs="*",
        help="optional repository-relative Markdown files; default scans all .md files",
    )
    args = parser.parse_args()

    paths, path_errors = resolve_paths(args.paths)
    total = 0
    problems: list[Problem] = []
    for path in paths:
        count, found = validate_markdown(path)
        total += count
        problems.extend(found)

    for message in path_errors:
        print(message, file=sys.stderr)
    for problem in problems:
        print(problem, file=sys.stderr)

    if path_errors or problems:
        print(
            f"FAIL: {len(path_errors) + len(problems)} problem(s) in "
            f"{total} math expression(s)",
            file=sys.stderr,
        )
        return 1

    print(f"PASS: {total} GitHub LaTeX expression(s) in {len(paths)} Markdown file(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
