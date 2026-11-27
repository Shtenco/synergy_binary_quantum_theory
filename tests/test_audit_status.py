"""Regression checks for negative controls and completion consistency."""
import contextlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_validator():
    spec = importlib.util.spec_from_file_location(
        "theory_gate_ledger", ROOT / "scripts" / "verify_theory_gates.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AuditStatusTests(unittest.TestCase):
    def run_ledger(self, module, data):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.json"
            module.LEDGER = path
            path.write_text(json.dumps(data), encoding="utf-8")
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                code = module.main()
            return code, json.loads(output.getvalue())

    def test_pure_ancilla_no_go_is_explicit_nonblocking_negative_control(self):
        module = load_validator()
        data = json.loads((ROOT / "theory_gates.json").read_text(encoding="utf-8"))
        gate = next(g for g in data["gates"] if g["id"] == "PURE_ANCILLA_NODE_REFINEMENT")
        self.assertEqual(gate["status"], "no_go")
        self.assertEqual(gate["closure_role"], "negative_control")

        code, result = self.run_ledger(module, data)
        self.assertEqual(code, 0)
        self.assertTrue(result["valid"])
        self.assertTrue(result["structural_candidate_closed"])
        self.assertEqual(result["negative_control_count"], 1)
        self.assertEqual(result["negative_controls"][0]["id"], "PURE_ANCILLA_NODE_REFINEMENT")
        self.assertFalse(result["negative_controls"][0]["blocking"])

    def test_no_go_cannot_be_misclassified_as_positive_core(self):
        module = load_validator()
        data = json.loads((ROOT / "theory_gates.json").read_text(encoding="utf-8"))
        gate = next(g for g in data["gates"] if g["id"] == "PURE_ANCILLA_NODE_REFINEMENT")
        gate["closure_role"] = "core"

        code, result = self.run_ledger(module, data)
        self.assertEqual(code, 1)
        self.assertFalse(result["valid"])
        self.assertTrue(any("status 'no_go' is not legal for closure_role 'core'" in e for e in result["errors"]))

    def test_negative_control_cannot_claim_a_positive_status(self):
        module = load_validator()
        data = json.loads((ROOT / "theory_gates.json").read_text(encoding="utf-8"))
        gate = next(g for g in data["gates"] if g["id"] == "PURE_ANCILLA_NODE_REFINEMENT")
        gate["status"] = "proved"

        code, result = self.run_ledger(module, data)
        self.assertEqual(code, 1)
        self.assertFalse(result["valid"])
        self.assertTrue(any("status 'proved' is not legal for closure_role 'negative_control'" in e for e in result["errors"]))

    def test_all_nine_physical_conditions_required(self):
        import ast

        def value(file, name):
            tree = ast.parse((ROOT / file).read_text(encoding="utf-8"))
            return next(
                ast.literal_eval(node.value)
                for node in tree.body
                if isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == name for target in node.targets)
            )

        actual = value("bcqg_bit_to_gravity_final.py", "REQUIRED_PHYSICAL")
        expected = value("scripts/verify_physicalization_gates.py", "EXPECTED_REQUIRED_PHYSICAL")
        self.assertEqual(actual, expected)
        self.assertEqual(len(actual), 9)


if __name__ == "__main__":
    unittest.main()
