#!/usr/bin/env python3
"""Regression tests for GitHub Markdown LaTeX delimiter normalization."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from validate_github_latex import validate_markdown


class GithubLatexValidatorTests(unittest.TestCase):
    def validate_text(self, source: str):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "case.md"
            path.write_text(source, encoding="utf-8")
            return validate_markdown(path)

    def test_multiline_bracket_display_math_is_one_expression(self):
        count, problems = self.validate_text(
            "Before\n\\[\n\\frac{a}{b}+c\n\\]\nAfter\n"
        )
        self.assertEqual(1, count)
        self.assertEqual([], problems)

    def test_inline_parenthesis_math_is_normalized(self):
        count, problems = self.validate_text(
            "The value is \\(\\frac{1}{2}\\).\n"
        )
        self.assertEqual(1, count)
        self.assertEqual([], problems)

    def test_unclosed_bracket_display_math_fails_closed(self):
        count, problems = self.validate_text(
            "Start\n\\[\n\\frac{a}{b}\n"
        )
        self.assertEqual(0, count)
        self.assertTrue(any("unclosed '$$' delimiter" in p.message for p in problems))

    def test_tex_command_outside_math_is_rejected(self):
        count, problems = self.validate_text("Unwrapped \\frac{a}{b}.\n")
        self.assertEqual(0, count)
        self.assertTrue(any("TeX command outside a math delimiter" in p.message for p in problems))


if __name__ == "__main__":
    unittest.main()
