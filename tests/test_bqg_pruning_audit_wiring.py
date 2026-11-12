import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_pruning_audit as audit
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_covariant_composition_gate as COMP
import peter_weyl_lorentzian_ordered_triple_gate as RAW


class TestPruningAuditOperatorWiring(unittest.TestCase):
    def test_prune_state_and_complete_to_covariant_projection_are_observed(self):
        spins = (4,) * 10
        with audit.capture() as collector:
            pruned = PW.prune_state({"small": 0.5e-8}, 1e-8)
            source_label = ("C", 0, 0, 0, 0)
            labels = (source_label,) + tuple(("G", 0) for _ in range(4))
            projected = COMP.complete_to_covariant(
                {(spins, labels): 0.5e-11}, 0, tol=1e-11
            )

        self.assertEqual(pruned, {})
        self.assertEqual(projected, {})
        sites = collector.snapshot()["by_site"]
        self.assertIn("k5.prune_state", sites)
        self.assertIn("covariant_composition.complete_to_covariant", sites)

    def test_shadow_mode_preserves_subthreshold_coefficients_at_wired_sites(self):
        with audit.capture() as collector:
            with audit.bypass():
                a = {}
                PW.add_dict(a, {"k": 0.5e-13})
                b = {}
                CK.add_state(b, {"k": 0.5e-11})
                c = {}
                COMP.add(c, {"k": 0.5e-11})
                d = {}
                RAW.add(d, {"k": 0.5e-11})
                e = PW.prune_state({"k": 0.5e-8}, 1e-8)

        self.assertEqual(a, {"k": 0.5e-13})
        self.assertEqual(b, {"k": 0.5e-11})
        self.assertEqual(c, {"k": 0.5e-11})
        self.assertEqual(d, {"k": 0.5e-11})
        self.assertEqual(e, {"k": 0.5e-8})
        events = collector.snapshot()["events"]
        self.assertEqual(len(events), 5)
        self.assertTrue(all(x["kind"] == "would_prune_bypassed" for x in events))

    def test_core_operator_merge_sites_are_observed_without_semantic_change(self):
        with audit.capture() as collector:
            a = {}
            PW.add_dict(a, {"k": 0.5e-13})
            b = {}
            CK.add_state(b, {"k": 0.5e-11})
            c = {}
            COMP.add(c, {"k": 0.5e-11})
            d = {}
            RAW.add(d, {"k": 0.5e-11})

        self.assertEqual(a, {})
        self.assertEqual(b, {})
        self.assertEqual(c, {})
        self.assertEqual(d, {})
        by_site = collector.snapshot()["by_site"]
        self.assertEqual(
            set(by_site),
            {
                "k5.add_dict.1e-13",
                "covariant_K_leg.add_state",
                "covariant_composition.add",
                "lorentzian_ordered_triple.add",
            },
        )
        self.assertTrue(all(row["coefficient_count"] == 1 for row in by_site.values()))


if __name__ == "__main__":
    unittest.main()
