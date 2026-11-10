#!/usr/bin/env python3
"""Fail CI when retired side sectors leak back into the active GR/QM core.

The repository has one canonical programme: discrete quantum microstructure ->
quantum geometry -> coarse-grained smooth geometry -> GR/HDA continuum limit.

Historical archives and the explicit repository-audit ledger are evidence about
past branches, so they are allowed to mention retired vocabulary. They are not
active theory surfaces and must not make the canonical scope gate fail merely
for preserving scientific history.
"""
from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELF = Path(__file__).resolve()
TEXT_SUFFIXES = {".py", ".md", ".yml", ".yaml", ".json"}

# Historical/evidentiary surfaces are quarantined from the active-claim scan.
# They remain version controlled and readable, but are not canonical theory.
EXCLUDED_FILES = {
    Path("REPOSITORY_AUDIT_2026-08-29.md"),
    Path("AUDIT_REPORT_2026-10-22.md"),
    Path("docs/README_TECHNICAL_ARCHIVE_20261109.md"),
    Path("docs/README_STORY_ARCHIVE_20261109.md"),
    Path("docs/FIGURE_ATLAS.md"),
}
EXCLUDED_PREFIXES = {
    ("docs", "archive"),
    ("docs", "information_infoton_review"),
}

# README and status may describe the explicitly open information-mode proposal.
# This is permission to discuss it, not to import retired physics modules.
DISCUSSION_TERMS = {
    # The README uses "mirror" only in archived figure labels/links; this is
    # permission to reference historical material, not to restore the sector.
    Path("README.md"): {"info" + "ton", "mir" + "ror"},
    Path("THEORY_STATUS.md"): {"info" + "ton"},
}

# Construct retired vocabulary without embedding the retired tokens literally in
# this audit source, so repository-wide textual searches remain meaningful.
RETIRED_TERMS = {
    "mir" + "ror",
    "anti" + "grav",
    "anti" + "matter",
    "info" + "ton",
    "gold" + "stone",
    "fifth" + " force",
}

# Every module removed by the scope reset is pinned here. AST parsing of every
# active Python file then prevents stale imports from surviving even if the
# retired concept is not mentioned in prose.
RETIRED_MODULES = {
    "bcqg_unified_verification",
    "bcqg_critical_phase_demo",
    "action_tunneling_blitz",
    "eml_blitz",
    "fermion_zero_blitz",
    "grav" + "iton_in" + "foton_foam_gate",
    "orientation_odd_hda_gate",
    "mir" + "ror_16cell_orientation_eta_gate",
    "mir" + "ror_chirality_gravity_gate",
    "mir" + "ror_force_normalization_gate",
    "mir" + "ror_goldstone_source_gate",
    "mir" + "ror_heisenberg_parent_gate",
    "mir" + "ror_hodge_stiffness_gate",
    "mir" + "ror_master_criterion_gate",
    "mir" + "ror_matter_matrix_element_gate",
    "mir" + "ror_order_16cell_gate",
    "mir" + "ror_order_recursive_pl_gate",
    "mir" + "ror_sigma_range_gate",
    "mir" + "ror_wilson_matter_gate",
}


def is_excluded(rel: Path) -> bool:
    if rel in EXCLUDED_FILES:
        return True
    parts = rel.parts
    return any(parts[: len(prefix)] == prefix for prefix in EXCLUDED_PREFIXES)


def candidate_files():
    for path in ROOT.rglob("*"):
        if not path.is_file() or path == SELF or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        if ".git" in path.parts or "verification_results" in path.parts:
            continue
        rel = path.relative_to(ROOT)
        if is_excluded(rel):
            continue
        yield path


def scan_text(path: Path):
    text = path.read_text(encoding="utf-8", errors="replace").lower()
    allowed = DISCUSSION_TERMS.get(path.relative_to(ROOT), set())
    return sorted(term for term in RETIRED_TERMS - allowed if term in text)


def local_imports(path: Path):
    if path.suffix.lower() != ".py":
        return set()
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module.split(".")[0])
    return names


def main() -> int:
    failures = []
    scanned = 0
    python_files = 0

    for path in candidate_files():
        scanned += 1
        rel = path.relative_to(ROOT)
        hits = scan_text(path)
        if hits:
            failures.append(f"{rel}: retired vocabulary: {', '.join(hits)}")

        if path.suffix.lower() == ".py":
            python_files += 1
            imports = local_imports(path)
            stale = sorted(imports & RETIRED_MODULES)
            if stale:
                failures.append(f"{rel}: imports retired module(s): {', '.join(stale)}")

    if failures:
        print("CORE SCOPE AUDIT: FAIL")
        for item in failures:
            print(f" - {item}")
        return 1

    print(
        "CORE SCOPE AUDIT: PASS "
        f"({scanned} active text/code files; {python_files} active Python files; "
        "historical archive quarantined; no retired vocabulary or retired-module imports)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
