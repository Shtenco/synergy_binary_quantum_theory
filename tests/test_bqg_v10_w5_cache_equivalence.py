import argparse
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from bqg_v10_w5_cache_equivalence import summarize_run
from bqg_v10_w5_amp_tol_shadow_compare import run_one


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

    def test_run_one_forwards_selected_seed_to_fresh_process(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "closure.json"
            args = argparse.Namespace(
                Jmax2=10, max_states=500, max_depth=8, spectral_zero_factor=1000.0,
                cache_mode="on", pruning_mode="production", nodes=[0],
                seed_Ks=[0], timeout_seconds=30, unbounded_closure=False,
            )
            captured = {}

            def fake_run(command, **kwargs):
                captured["command"] = command
                output.write_text(json.dumps({"status": "INCOMPLETE_CLOSURE_LIMIT_REACHED"}))
                return type("Completed", (), {"returncode": 2, "stderr": "", "stdout": ""})()

            with patch("bqg_v10_w5_amp_tol_shadow_compare.subprocess.run", side_effect=fake_run):
                run_one(Path("closure.py"), output, 1e-12, args)
            command = captured["command"]
            seed_index = command.index("--seed-K")
            self.assertEqual(command[seed_index + 1], "0")
            self.assertIn("--nodes", command)

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
