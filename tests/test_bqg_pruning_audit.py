import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_pruning_audit as audit


class TestPruningAuditCollector(unittest.TestCase):
    def test_incremental_merge_records_removed_partial_sum(self):
        with audit.capture() as collector:
            audit.record_merge(
                "test.merge", ("state", 1), 0.75e-11, -0.25e-11,
                0.5e-11, 1e-11,
            )
        result = collector.snapshot()
        self.assertEqual(result["event_count"], 1)
        event = result["events"][0]
        self.assertEqual(event["kind"], "incremental_merge_result")
        self.assertAlmostEqual(event["real"], 0.5e-11)
        self.assertAlmostEqual(event["prior_real"], 0.75e-11)
        self.assertAlmostEqual(event["incoming_real"], -0.25e-11)
        self.assertEqual(result["by_site"]["test.merge"]["coefficient_count"], 1)

    def test_exact_cancellation_is_retained_as_zero_mass_event(self):
        with audit.capture() as collector:
            audit.record_merge("test.cancel", "k", 1 + 2j, -1 - 2j, 0j, 1e-11)
        result = collector.snapshot()
        self.assertEqual(result["event_count"], 1)
        self.assertEqual(result["events"][0]["abs"], 0.0)
        self.assertEqual(result["events"][0]["kind"], "exact_cancellation")
        self.assertEqual(result["by_site"]["test.cancel"]["pathwise_l1"], 0.0)

    def test_threshold_records_only_discarded_values(self):
        with audit.capture() as collector:
            audit.record_threshold("test.threshold", "small", 0.5e-12, 1e-12)
            audit.record_threshold("test.threshold", "large", 2e-12, 1e-12)
        result = collector.snapshot()
        self.assertEqual(result["event_count"], 1)
        self.assertEqual(result["events"][0]["key_repr"], "'small'")

    def test_structural_zero_is_not_counted_as_amplitude_loss(self):
        with audit.capture() as collector:
            audit.record_threshold("test.zero", "zero", 0j, 1e-12)
        self.assertEqual(collector.snapshot()["event_count"], 0)

    def test_shadow_bypass_keeps_nonzero_subthreshold_values(self):
        with audit.capture() as collector:
            with audit.bypass():
                self.assertTrue(audit.keep(0.5e-12, 1e-12))
                self.assertFalse(audit.keep(0j, 1e-12))
                audit.record_threshold("test.shadow", "small", 0.5e-12, 1e-12)
        event = collector.snapshot()["events"][0]
        self.assertEqual(event["kind"], "would_prune_bypassed")

    def test_snapshot_separates_production_shadow_and_cancellation_kinds(self):
        with audit.capture() as collector:
            audit.record_threshold("site", "prod", 0.5e-12, 1e-12)
            with audit.bypass():
                audit.record_threshold("site", "shadow", 0.5e-12, 1e-12)
            audit.record_merge("site", "cancel", 1 + 0j, -1 + 0j, 0j, 1e-12)
        by_kind = collector.snapshot()["by_site_and_kind"]
        self.assertIn("site::coefficient_threshold", by_kind)
        self.assertIn("site::would_prune_bypassed", by_kind)
        self.assertIn("site::exact_cancellation", by_kind)
        self.assertEqual(by_kind["site::exact_cancellation"]["pathwise_l1"], 0.0)

    def test_collector_is_disabled_outside_capture(self):
        audit.record_merge("test.inactive", "k", 1, -1, 0, 1e-11)
        with audit.capture() as collector:
            pass
        self.assertEqual(collector.snapshot()["event_count"], 0)

    def test_nested_capture_restores_outer_collector(self):
        with audit.capture() as outer:
            audit.record_threshold("outer.before", "a", 1e-13, 1e-12)
            with audit.capture() as inner:
                audit.record_threshold("inner", "b", 1e-13, 1e-12)
            audit.record_threshold("outer.after", "c", 1e-13, 1e-12)
        self.assertEqual(outer.snapshot()["event_count"], 2)
        self.assertEqual(inner.snapshot()["event_count"], 1)


if __name__ == "__main__":
    unittest.main()
