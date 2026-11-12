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
