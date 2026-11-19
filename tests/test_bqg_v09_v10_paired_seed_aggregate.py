import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from bqg_v09_v10_paired_seed_aggregate import aggregate, CUTOFFS, MODES, SEEDS


def make_report(cutoff, mode, seed, perturb=0.0):
    value = complex(seed + 1 + perturb, 0.25 * seed)
    row = [repr(("state", seed)), value.real, value.imag]
    return {
        "Jmax2": cutoff,
        "mode": mode,
        "seed_Ks": [seed],
        "acceptance": {"paired_seed": True},
        "V09_vs_V10_forward_columns": [{
            "seed_K": seed,
            "V09_raw_sparse": [row],
            "V10_seed_sparse": [row],
        }],
        "V10": {
            "basis_id_sha256": "same-basis",
            "closure_size": 1,
            "forward_matrix_digest": "same-forward",
            "direct_adjoint_matrix_digest": "same-adjoint",
        },
    }


class TestSeedShardedPairedAggregate(unittest.TestCase):
    def fixture(self):
        return [make_report(j, mode, seed)
                for j in CUTOFFS for mode in MODES for seed in SEEDS]

    def test_complete_seed_sharded_grid_passes(self):
        result = aggregate(self.fixture())
        self.assertTrue(result["passed"])
        self.assertEqual(result["expected_report_count"], 30)
        self.assertEqual(len(result["production_vs_no_pruning"]), 15)
        self.assertEqual(len(result["cutoff_convergence"]), 20)

    def test_production_shadow_difference_fails(self):
        rows = self.fixture()
        idx = next(i for i, r in enumerate(rows) if r["Jmax2"] == 10 and r["mode"] == "production" and r["seed_Ks"] == [0])
        rows[idx] = make_report(10, "production", 0, perturb=1e-4)
        result = aggregate(rows)
        self.assertFalse(result["passed"])
        self.assertTrue(result["pruning_failures"])

    def test_cutoff_difference_fails(self):
        rows = self.fixture()
        idx = next(i for i, r in enumerate(rows)
                   if r["Jmax2"] == 12 and r["mode"] == "shadow_no_prune" and r["seed_Ks"] == [0])
        rows[idx] = make_report(12, "shadow_no_prune", 0, perturb=1e-4)
        result = aggregate(rows)
        self.assertFalse(result["passed"])
        self.assertTrue(result["cutoff_convergence_failures"])

    def test_missing_seed_cutoff_mode_cell_fails_closed(self):
        result = aggregate(self.fixture()[:-1])
        self.assertFalse(result["passed"])
        self.assertEqual(len(result["missing_cells"]), 1)


if __name__ == "__main__":
    unittest.main()
