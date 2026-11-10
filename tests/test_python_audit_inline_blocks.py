from __future__ import annotations

import ast
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_python_all_branch_heads import iter_inline_python_blocks


class InlinePythonAuditTests(unittest.TestCase):
    def test_markdown_python_fences_are_extracted_and_dedented(self):
        source = """# Example
\x60\x60\x60python title="safe"
    import math
    value = math.sqrt(4)
\x60\x60\x60

\x60\x60\x60bash
python -c 'print(1)'
\x60\x60\x60
"""
        blocks = list(iter_inline_python_blocks(source, "README.md"))
        self.assertEqual(len(blocks), 1)
        self.assertEqual(blocks[0][0], "fence-python")
        tree = ast.parse(blocks[0][2])
        self.assertTrue(tree.body)

    def test_yaml_python_heredocs_are_extracted_and_dedented(self):
        source = """steps:
  - run: |
      python - <<'PY'
      import json
      value = {"ok": True}
      print(json.dumps(value))
      PY
"""
        blocks = list(iter_inline_python_blocks(source, ".github/workflows/example.yml"))
        self.assertEqual(len(blocks), 1)
        self.assertEqual(blocks[0][0], "heredoc-PY")
        tree = ast.parse(blocks[0][2])
        self.assertTrue(tree.body)

    def test_multiple_python_fences_are_kept_separate(self):
        source = """\x60\x60\x60py
x = 1
\x60\x60\x60
\x60\x60\x60python3
y = 2
\x60\x60\x60
"""
        blocks = list(iter_inline_python_blocks(source, "notes.md"))
        self.assertEqual([x[0] for x in blocks], ["fence-py", "fence-python3"])
        self.assertEqual([ast.parse(x[2]).body[0].targets[0].id for x in blocks], ["x", "y"])


if __name__ == "__main__":
    unittest.main()
