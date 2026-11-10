import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github/workflows/bqg-depth6-mixed-stage-b-final.yml"
LEDGER = ROOT / "depth6_frontier.json"


class Depth6StageBWorkflowLedgerTests(unittest.TestCase):
    def test_stage_b_workflow_does_not_reopen_closed_s4_sector(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("s4['sector_status']=='CLOSED_FINITE_NUMERICAL'", workflow)
        self.assertIn("s4['obsolete_gate']['status']=='SUPERSEDED_DO_NOT_USE_AS_BLOCKER'", workflow)
        self.assertNotIn("s4['sector_status']=='ACTIVE_NOT_CLOSED'", workflow)
        self.assertNotIn("s4['giant_component']['square_minor_plan']", workflow)
        self.assertNotIn("giant_sparse_minor_shape", workflow)

    def test_stage_b_output_preserves_full_theorem_boundary(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("'finite_depth6_theorem_status':'NOT_YET_PROVED'", workflow)
        self.assertIn("'status':'NUMERICAL_WITNESSES_PENDING'", workflow)
        self.assertIn("continuum limit", workflow)

    def test_canonical_ledger_still_marks_s4_closed_and_three_witnesses_pending(self):
        import json
        ledger = json.loads(LEDGER.read_text(encoding="utf-8"))
        s4 = ledger["s4_sign_frontier"]
        self.assertEqual(s4["sector_status"], "CLOSED_FINITE_NUMERICAL")
        self.assertEqual(s4["obsolete_gate"]["status"], "SUPERSEDED_DO_NOT_USE_AS_BLOCKER")
        self.assertEqual(ledger["remaining_numerical_irreps"], ["[3,2]", "[3,1,1]", "[2,2,1]"])
        self.assertEqual(ledger["finite_depth6_theorem_status"], "NOT_YET_PROVED")


if __name__ == "__main__":
    unittest.main()
