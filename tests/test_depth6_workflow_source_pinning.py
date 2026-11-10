from __future__ import annotations

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_SHA = "f6bfcd775b8c3cd33c6c235d48293e40657baf9d"
BUNDLE_SHA = "cbbaa80f6b6b61d353458db1fa1f03f824730532"


class Depth6WorkflowSourcePinningTests(unittest.TestCase):
    def test_mixed_stage_workflows_pin_source_commit_and_bundle(self):
        paths = (
            ROOT / ".github/workflows/bqg-depth6-mixed-stage-a.yml",
            ROOT / ".github/workflows/bqg-depth6-mixed-stage-b-final.yml",
            ROOT / ".github/workflows/bqg-depth6-32-keymeta-salvage.yml",
        )
        for path in paths:
            with self.subTest(path=path.name):
                source = path.read_text(encoding="utf-8")
                self.assertNotIn("ref: research/depth6-s4sign-proof", source)
                self.assertIn(SOURCE_SHA, source)
                self.assertIn(BUNDLE_SHA, source)

    def test_s4sign_aggregate_pins_the_triggering_producer_commit(self):
        path = ROOT / ".github/workflows/bqg-depth6-s4sign-aggregate.yml"
        source = path.read_text(encoding="utf-8")
        self.assertIn("ref: ${{ github.event.workflow_run.head_sha }}", source)
        self.assertIn("test \"${actual}\" = \"${EXPECTED_HEAD}\"", source)


if __name__ == "__main__":
    unittest.main()
