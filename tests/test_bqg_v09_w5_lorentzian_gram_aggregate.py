import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from bqg_v09_w5_lorentzian_gram_aggregate import assemble_w5_gram

class W5GramAggregationTests(unittest.TestCase):
    def test_cross_gram_preserves_interference_and_shard_coverage(self):
        records = []
        for node in range(5):
            for col, K in enumerate((0, 2, 4, 6, 8)):
                key = f"state-{col}"
                records.append({
                    "source_vertex": node,
                    "input_K": K,
                    "Jmax2": 10,
                    "physical_status": {"hermitian_completion_of_HL": False},
                    "A_Euclidean_sine_Gauss": [[key, 1.0, 0.0]],
                    "B_Lorentzian_raw_KKV_Gauss": [[key, -1.0, 0.0]],
                    "B_scalar_projection_rejections": [],
                    "term_diagnostics": [{"term_index": i} for i in range(24)],
                })
        result = assemble_w5_gram(records)
        GEE = result["gram_matrices"]["G_EE"]
        GEL = result["gram_matrices"]["G_EL_raw"]
        GLL = result["gram_matrices"]["G_LL_raw"]
        self.assertAlmostEqual(GEE[0][0][0], 5.0)
        self.assertAlmostEqual(GEL[0][0][0], -5.0)
        self.assertAlmostEqual(GLL[0][0][0], 5.0)

    def test_rejects_missing_epsilon_terms(self):
        records = []
        for node in range(5):
            for col, K in enumerate((0, 2, 4, 6, 8)):
                records.append({
                    "source_vertex": node,
                    "input_K": K,
                    "Jmax2": 10,
                    "physical_status": {"hermitian_completion_of_HL": False},
                    "A_Euclidean_sine_Gauss": [[f"state-{col}", 1.0, 0.0]],
                    "B_Lorentzian_raw_KKV_Gauss": [],
                    "B_scalar_projection_rejections": [],
                    "term_diagnostics": [{"term_index": i} for i in range(23)],
                })
        with self.assertRaises(ValueError):
            assemble_w5_gram(records)

if __name__ == "__main__":
    unittest.main()
