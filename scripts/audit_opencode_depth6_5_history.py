#!/usr/bin/env python3
"""Read-only literal audit for historical depth-6 [5] claims in OpenCode.

This scans JSON string values in the local OpenCode SQLite message store. It
prints message identifiers and matched terms only; it never prints message
contents. SQLite is opened immutable/read-only so the audit cannot alter the
history database.
"""
from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path

TERMS = ("nonvacuum_rank", "A_[5]", "B_[5]")


def strings(value: object):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for item in value.values():
            yield from strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from strings(item)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db",
        type=Path,
        default=Path.home() / ".local" / "share" / "opencode" / "opencode.db",
        help="OpenCode SQLite database (default: ~/.local/share/opencode/opencode.db)",
    )
    args = parser.parse_args()
    db = args.db.expanduser().resolve()
    if not db.is_file():
        parser.error(f"OpenCode database not found: {db}")

    uri = db.as_uri() + "?mode=ro&immutable=1"
    connection = sqlite3.connect(uri, uri=True)
    matches = {term: [] for term in TERMS}
    try:
        placeholders = " OR ".join("instr(data, ?) > 0" for _ in TERMS)
        rows = connection.execute(
            f"SELECT session_id, type, seq, data FROM session_message WHERE {placeholders}",
            TERMS,
        )
        for session_id, message_type, seq, data in rows:
            try:
                payload = json.loads(data or "null")
            except (TypeError, json.JSONDecodeError):
                payload = data
            values = list(strings(payload))
            found = [term for term in TERMS if any(term.casefold() in value.casefold() for value in values)]
            for term in found:
                matches[term].append(
                    {"session_id": session_id, "type": message_type, "seq": seq}
                )
    finally:
        connection.close()

    print(
        json.dumps(
            {
                "status": "READ_ONLY_LITERAL_SCAN_COMPLETE",
                "database_source": "local OpenCode SQLite",
                "terms": {
                    term: {"message_count": len(records), "messages": records}
                    for term, records in matches.items()
                },
                "scope_note": "Exact literal search across JSON string values only; absence does not rule out differently named or external evidence.",
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
