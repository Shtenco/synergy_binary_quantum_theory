import unittest

from scripts.bqg_v11_w5_common_carrier import build_carrier, joint_gram_input


class FakeV09:
    def __init__(self, mismatch=False, small_amplitude=None):
        self.mismatch = mismatch
        self.small_amplitude = small_amplitude

    def run_column(self, node, input_K, Jmax2=10, start=0, count=24):
        spins = (4,) * 10
        a_key = (spins, (input_K, 0, 0, 0, 1))
        b_key = (spins, (input_K, 0, 0, 0, 2))
        b_amp = self.small_amplitude if self.small_amplitude is not None else 2.0 + node
        if self.mismatch:
            b_amp += 0.5
        return {
            "status": "COMPUTED_RAW_W5_KKV_COLUMN",
            "epsilon_terms_completed": 24,
            "sine_wrapper_gate_passed": True,
            "sine_wrapper_diagnostics": {},
            "A_Euclidean_sine_Gauss": [[repr(a_key), 1.0 + node, 0.0]],
            "B_Lorentzian_raw_KKV_Gauss": [[repr(b_key), b_amp, 0.0]],
            "B_Lorentzian_raw_KKV_covariant": [[repr((spins, (-1, 0, 0, 0, 0), 0, 0, 0, 0)), b_amp, 0.0]],
            "A_support": 1,
            "B_support": 1,
            "A_norm": 1.0 + node,
            "B_raw_norm": b_amp,
            "B_non_scalar_covariant_norm": 0.0,
            "B_scalar_projection_rejections": [],
            "term_diagnostics": [{"max_basis_leakage": 0.0}],
        }


class FakeV10:
    def __init__(self, status="PASS_ADJOINT_CONSISTENCY", forward_amplitude=20.0):
        self.status = status
        self.forward_amplitude = forward_amplitude

    @staticmethod
    def w5_input_key(K):
        return ((4,) * 10, (K, 0, 0, 0, 0))

    def run_closure(self, source_nodes, Jmax2, max_states, max_depth, amp_tol):
        if self.status != "PASS_ADJOINT_CONSISTENCY":
            return {"status": self.status, "partial_diagnostics": {}}
        spins = (4,) * 10
        outputs = [(spins, (K, 0, 0, 0, 2)) for K in (0, 2, 4, 6, 8)]
        seeds = [self.w5_input_key(K) for K in (0, 2, 4, 6, 8)]
        keys = sorted(seeds + outputs, key=repr)
        ix = {k: i for i, k in enumerate(keys)}
        cols = [[] for _ in keys]
        sym = [[] for _ in keys]
        anti = [[] for _ in keys]
        for K in (0, 2, 4, 6, 8):
            key = (spins, (K, 0, 0, 0, 2))
            col = ix[self.w5_input_key(K)]
            cols[col] = [[ix[key], self.forward_amplitude, 0.0]]
            sym[col] = [[ix[key], self.forward_amplitude, 0.0]]
        return {
            "status": self.status,
            "basis_keys_repr": [repr(k) for k in keys],
            "basis_id_sha256": "fake-v10-basis",
            "closure_size": len(keys),
            "closure_depth_max": 1,
            "forward_B_columns_sparse": cols,
            "hermitian_ordering_candidate": {
                "B_comm_symmetric_columns": sym,
                "B_comm_antihermitian_columns": anti,
            },
            "diagnostics": {
                "direct_adjoint_relative_residual": 0.0,
                "max_component_leakage": 0.0,
            },
        }


class TestW5CommonCarrier(unittest.TestCase):
    def test_assembles_common_carrier_and_crosschecks_node_sum(self):
        result = build_carrier(
            source_nodes=(0, 1, 2, 3, 4),
            v09_module=FakeV09(),
            v10_module=FakeV10(),
        )
        self.assertEqual(result["status"], "PASS_STRUCTURAL_COMMON_CARRIER")
        self.assertTrue(result["crosschecks"]["all_passed"])
        self.assertTrue(all(x["passed"] for x in result["crosschecks"]["V10_forward_equals_symmetric_plus_anti"]))
        self.assertEqual(len(result["A"][0]), 5)
        self.assertEqual(len(result["B"][0]), 5)
        # A is coherently summed over five source nodes: 1+2+3+4+5 = 15.
        a_rows = result["common_output_basis"]["keys_repr"]
        a_idx = a_rows.index(repr(((4,) * 10, (0, 0, 0, 0, 1))))
        self.assertEqual(result["A"][a_idx][0], [15.0, 0.0])
        self.assertEqual(result["B"][a_rows.index(repr(((4,) * 10, (0, 0, 0, 0, 2))))][0], [20.0, 0.0])

    def test_threshold_is_applied_after_coherent_node_sum(self):
        amp_tol = 1e-12
        node_amp = 0.6 * amp_tol
        result = build_carrier(
            amp_tol=amp_tol,
            v09_module=FakeV09(small_amplitude=node_amp),
            v10_module=FakeV10(forward_amplitude=5 * node_amp),
        )
        self.assertEqual(result["status"], "PASS_STRUCTURAL_COMMON_CARRIER")
        self.assertTrue(result["crosschecks"]["all_passed"])
        self.assertAlmostEqual(
            result["crosschecks"]["V09_per_node_then_node_sum_vs_V10_all_node_forward"][0]
            ["v09_node_sum_vs_v10_forward_relative_residual"],
            0.0,
        )

    def test_amp_tol_cannot_fall_below_v09_serialization_accumulation_floor(self):
        with self.assertRaisesRegex(ValueError, "serialization accumulation floor"):
            build_carrier(
                amp_tol=5e-13,
                v09_module=FakeV09(),
                v10_module=FakeV10(),
            )

    def test_v10_incomplete_closure_blocks_carrier(self):
        result = build_carrier(
            v09_module=FakeV09(),
            v10_module=FakeV10("INCOMPLETE_CLOSURE_LIMIT_REACHED"),
        )
        self.assertEqual(result["status"], "BLOCKED_V10_CLOSURE")

    def test_v09_v10_raw_amplitude_mismatch_blocks_carrier(self):
        result = build_carrier(
            v09_module=FakeV09(mismatch=True),
            v10_module=FakeV10(),
        )
        self.assertEqual(result["status"], "FAIL_V09_V10_RAW_AMPLITUDE_CROSSCHECK")
        self.assertFalse(result["crosschecks"][0]["passed"])

    def test_joint_input_requires_explicit_provenance(self):
        carrier = build_carrier(
            v09_module=FakeV09(),
            v10_module=FakeV10(),
        )
        with self.assertRaisesRegex(ValueError, "normalization_source"):
            joint_gram_input(carrier, 1.0, 1.0, "")

    def test_joint_input_preserves_basis_identity_and_coefficients(self):
        carrier = build_carrier(
            v09_module=FakeV09(),
            v10_module=FakeV10(),
        )
        payload = joint_gram_input(carrier, 1.0, -0.25, "explicit synthetic test convention")
        self.assertEqual(payload["output_basis_id"], payload["A_output_basis_id"])
        self.assertEqual(payload["output_basis_id"], payload["B_output_basis_id"])
        self.assertEqual(payload["normalization"]["c_L"], -0.25)


if __name__ == "__main__":
    unittest.main()
