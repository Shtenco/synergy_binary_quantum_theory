import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import peter_weyl_lorentzian_sine_ordered_triple_gate as sine_gate


class TestLorentzianSineExecutorDispatch(unittest.TestCase):
    def test_explicit_executor_is_invoked_instead_of_raw_run(self):
        sentinel = {"passed": True, "source": "executor"}
        callback = lambda: sentinel
        with patch.object(sine_gate.RAW, "run", side_effect=AssertionError("RAW.run must not execute")):
            self.assertIs(sine_gate._dispatch_raw(0, None, None, None, executor=callback), sentinel)

    def test_legacy_call_without_executor_preserves_raw_run(self):
        sentinel = {"passed": True, "source": "legacy"}
        with patch.object(sine_gate.RAW, "run", return_value=sentinel) as raw_run:
            self.assertIs(sine_gate._dispatch_raw(3, 1, 2, 4), sentinel)
            raw_run.assert_called_once_with(3, 1, 2, 4)

    def test_executor_result_must_be_diagnostics_mapping(self):
        with patch.object(sine_gate.RAW, "run", return_value=None):
            # Dispatch itself preserves the return value; run() validates the
            # mapping before applying the wrapper diagnostics.
            self.assertIsNone(sine_gate._dispatch_raw(0, None, None, None))


if __name__ == "__main__":
    unittest.main()
