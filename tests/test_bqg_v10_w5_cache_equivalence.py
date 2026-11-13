import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from bqg_v10_w5_cache_equivalence import summarize_run


class TestW5CacheEquivalenceRunner(unittest.TestCase):
    def test_incomplete_timeout_without_ledgers_is_summarized_safely(self):
        summary = summarize_run({
            "status": "INCOMPLETE_TIMEOUT",
            "_process": {"timeout_seconds": 900},
        })
        self.assertEqual(summary["status"], "INCOMPLETE_TIMEOUT")
        self.assertIsNone(summary["cache_policy"])
        self.assertIsNone(summary["event_count"])
        self.assertEqual(summary["process"]["timeout_seconds"], 900)

    def test_completed_run_includes_cache_policy_and_event_count(self):
        summary = summarize_run({
            "status": "PASS_ADJOINT_CONSISTENCY",
            "closure_size": 12,
            "basis_id_sha256": "abc",
            "lower_level_pruning_ledger": {
                "cache_policy": {"mode": "on"},
                "event_count": 4,
            },
        })
        self.assertEqual(summary["closure_size"], 12)
        self.assertEqual(summary["cache_policy"], {"mode": "on"})
        self.assertEqual(summary["event_count"], 4)


if __name__ == "__main__":
    unittest.main(verbosity=2)
