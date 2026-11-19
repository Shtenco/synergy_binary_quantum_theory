import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from bqg_v09_v10_paired_aggregate import aggregate, CUTOFFS, MODES, SEEDS


def make_report(cutoff, mode, perturb=0.0):
    columns = []
    for seed in SEEDS:
        value = complex(seed + 1 + perturb, 0.25 * seed)
        row = [repr(("state", seed)), value.real, value.imag]
        columns.append({
            "seed_K": seed,
            "V09_raw_sparse": [row],
            "V10_seed_sparse": [row],
        })
    return {
        "Jmax2": cutoff,
        "mode": mode,
        "acceptance": {"paired": True},
        "V09_columns": [],
        "V09_vs_V10_forward_columns": columns,
        "V10": {
            "basis_id_sha256": "same-basis",
            "closure_size": 5,
            "forward_matrix_digest": "same-forward",
            "direct_adjoint_matrix_digest": "same-adjoint",
        },
    }


class TestPairedAggregate(unittest.TestCase):
    def fixture(self):
        return [make_report(c, m) for c in CUTOFFS for m in MODES]

    def test_complete_identical_matrix_passes(self):
        report = aggregate(self.fixture())
        self.assertTrue(report["passed"])
        self.assertEqual(report["expected_report_count"], 6)

    def test_production_shadow_difference_fails(self):
        rows = self.fixture()
        rows[1] = make_report(10, "production", perturb=1e-4)
        report = aggregate(rows)
        self.assertFalse(report["passed"])
        self.assertTrue(report["pruning_failures"])

    def test_cutoff_difference_fails(self):
        rows = self.fixture()
        rows[2] = make_report(12, "shadow_no_prune", perturb=1e-4)
        report = aggregate(rows)
        self.assertFalse(report["passed"])
        self.assertTrue(report["cutoff_convergence_failures"])

    def test_missing_cell_fails_closed(self):
        report = aggregate(self.fixture()[:-1])
        self.assertFalse(report["passed"])
        self.assertEqual(len(report["missing_cells"]), 1)


if __name__ == "__main__":
    unittest.main()
