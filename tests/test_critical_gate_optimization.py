from __future__ import annotations

import ast
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "bqg_canonical_refinement_residual_gate.py"

sys.path.insert(0, str(ROOT / "scripts"))
from bqg_canonical_refinement_residual_gate import require


class CriticalGateOptimizationTests(unittest.TestCase):
    def test_require_raises_even_when_python_assertions_are_optimized_out(self):
        with self.assertRaisesRegex(AssertionError, "deliberate failure"):
            require(False, "deliberate failure")

    def test_canonical_gate_contains_no_optimization_sensitive_asserts(self):
        tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
        self.assertFalse(
            any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "canonical refinement gate must use explicit runtime checks",
        )

    def test_canonical_gate_passes_under_python_optimization(self):
        result = subprocess.run(
            [sys.executable, "-O", str(SCRIPT)],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
            timeout=60,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("RESIDUAL DECOMPOSITION PASS", result.stdout)
        self.assertIn("ORTHOGONAL DEFECT SPLIT PASS", result.stdout)
        self.assertIn("PASS", result.stdout)

    def test_volume_block_gate_contains_no_optimization_sensitive_asserts(self):
        volume_script = ROOT / "docs" / "numerical_continuation" / "verify_volume_blocks.py"
        tree = ast.parse(volume_script.read_text(encoding="utf-8"))
        self.assertFalse(
            any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "exact volume-block gate must use explicit runtime checks",
        )

    def test_volume_block_gate_passes_under_python_optimization(self):
        volume_script = ROOT / "docs" / "numerical_continuation" / "verify_volume_blocks.py"
        result = subprocess.run(
            [sys.executable, "-O", str(volume_script)],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
            timeout=120,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("exact_local_blocks", result.stdout)
        self.assertIn('"physicalization_complete": false', result.stdout)

    def test_assert_based_scientific_gates_refuse_optimized_execution(self):
        paths = (
            "scripts/bqg_three_node_reduced_chain_gate.py",
            "scripts/bqg_local_information_resistance_gate.py",
            "scripts/bqg_two_node_reduced_gluing_gate.py",
            "scripts/bqg_projector_graph_transition_gate.py",
            "scripts/bqg_minimal_entanglement_volume_gate.py",
        )
        for relative in paths:
            with self.subTest(path=relative):
                tree = ast.parse((ROOT / relative).read_text(encoding="utf-8"))
                guards = [
                    node for node in tree.body
                    if isinstance(node, ast.If)
                    and isinstance(node.test, ast.UnaryOp)
                    and isinstance(node.test.op, ast.Not)
                    and isinstance(node.test.operand, ast.Name)
                    and node.test.operand.id == "__debug__"
                ]
                self.assertTrue(
                    guards,
                    f"{relative} must refuse python -O while it contains assert gates",
                )


if __name__ == "__main__":
    unittest.main()
