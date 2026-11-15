import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import peter_weyl_w5_direct_all_M_su2_gate as gate
import peter_weyl_covariant_K_sine_composition_gate as sine_gate


class TestDirectAllMSU2FailClosed(unittest.TestCase):
    def test_empty_commutator_inputs_are_rejected(self):
        with self.assertRaisesRegex(AssertionError, "empty state"):
            gate.relative_error({}, {})

    def test_zero_norm_commutator_inputs_are_rejected(self):
        zero = {"state": 0j}
        with self.assertRaisesRegex(AssertionError, "zero or invalid norm"):
            gate.relative_error(zero, zero)

    def test_nonzero_output_is_required(self):
        with self.assertRaisesRegex(AssertionError, "empty sparse state"):
            gate.require_nonzero_output({}, label="synthetic H_E")
        with self.assertRaisesRegex(AssertionError, "nonzero norm"):
            gate.require_nonzero_output({"state": 0j}, label="synthetic H_E")

    def test_tiny_but_nonzero_operator_output_is_not_misclassified_as_empty(self):
        state = {"state": 1.0e-200 + 0j}
        self.assertEqual(gate.require_nonzero_output(state, label="tiny H_E"), 1.0e-200)
        rel, absolute, scale = gate.relative_error(state, state)
        self.assertEqual(rel, 0.0)
        self.assertEqual(absolute, 0.0)
        self.assertEqual(scale, 1.0e-200)

    def test_nonfinite_output_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "non-finite amplitude"):
            gate.require_nonzero_output({"state": complex(float("inf"), 0)}, label="H_E")

    def test_wrong_sector_fraction_uses_full_output_weight(self):
        state = {
            ((4,) * 10, (("C", 2, 0, 0, 2),)): 1.0 + 0j,
            ((4,) * 10, (("C", 4, 0, 0, 2),)): 1.0 + 0j,
        }
        fraction = gate.wrong_sector_fraction(
            state, node=0, expected_j2=2, expected_m2=0
        )
        self.assertEqual(fraction, 0.5)

    def test_make_input_places_the_charged_sector_on_each_requested_node(self):
        for node in range(5):
            _, labels = gate.make_input(0, node)
            self.assertEqual(labels[node], ("C", 2, 0, 0, 2))
            self.assertEqual(sum(label[0] == "C" for label in labels), 1)

    def test_make_input_rejects_invalid_magnetic_component(self):
        for m2 in (-4, 4, 1, True):
            with self.subTest(m2=m2):
                with self.assertRaises(ValueError):
                    gate.make_input(m2, 0)

    def test_make_input_rejects_invalid_node(self):
        for node in (-1, 5, True, 1.5):
            with self.subTest(node=node):
                with self.assertRaises(ValueError):
                    gate.make_input(0, node)

    def test_run_rejects_duplicate_or_empty_node_selection(self):
        with self.assertRaisesRegex(ValueError, "non-empty subset"):
            gate.run(nodes=())
        with self.assertRaisesRegex(ValueError, "duplicates"):
            gate.run(nodes=(0, 0))

    def test_all_five_nodes_receive_four_adjacent_m_commutator_checks(self):
        def identity_map(key, source_v, Jmax2, charged_nodes):
            self.assertEqual(charged_nodes, (source_v,))
            return {key: 1.0 + 0j}, 0.0, 0.0

        def identity_superposition(state, source_v, Jmax2, charged_nodes):
            self.assertEqual(charged_nodes, (source_v,))
            return dict(state), 0.0, 0.0

        with patch.object(gate.SINE, "complete_HE_sine", side_effect=identity_map), \\
             patch.object(gate.SINE, "complete_HE_sine_state", side_effect=identity_superposition):
            result = gate.run(Jmax2=4, nodes=(0, 1, 2, 3, 4))

        self.assertTrue(result["passed"], result)
        self.assertEqual(result["source_nodes"], [0, 1, 2, 3, 4])
        self.assertEqual(result["input_charged_sector"]["channels_K12_K34"], [[0, 2], [2, 0]])
        self.assertEqual(result["magnetic_components_tested_M2"], [-2, 0, 2])
        self.assertLess(result["max_adjointness_compression_relative_residual"], 1e-12)
        self.assertLess(result["max_matrix_superposition_identity_relative_residual"], 1e-12)
        self.assertEqual(len(result["actual_complex_superposition_rows"]), 5)
        self.assertLess(result["max_actual_complex_superposition_relative_residual"], 1e-12)
        self.assertEqual(len(result["cases"]), 40)
        self.assertEqual(len(result["operator_columns_sparse"]), 30)
        self.assertTrue(all(len(column["entries"]) >= 1 for column in result["operator_columns_sparse"]))
        self.assertEqual(
            {(column["source_node"], column["K12"], column["K34"], column["input_M2"])
             for column in result["operator_columns_sparse"]},
            {(node, k12, k34, m2) for node in range(5)
             for k12, k34 in ((0, 2), (2, 0)) for m2 in (-2, 0, 2)},
        )
        self.assertEqual(
            {case["source_node"] for case in result["cases"]},
            {0, 1, 2, 3, 4},
        )
        self.assertEqual(
            {(case["direction"], case["input_M2"]) for case in result["cases"]},
            {
                ("J_plus", -2), ("J_plus", 0),
                ("J_minus", 0), ("J_minus", 2),
            },
        )

    def test_complex_superposition_helper_propagates_amplitudes_into_branches(self):
        key = gate.make_input(0, 0, 0, 2)

        def branch_from_key(input_key, amp=1.0 + 0j):
            return (input_key, (), complex(amp))

        def apply_sequence(branch, seq, source_v, Jmax2):
            return [(seq, branch[2])], 0.0

        def project_branch(branch, charged_nodes):
            seq, amp = branch
            return {("output", seq): amp}, 0.0

        with patch.object(sine_gate.KC.CK, "branch_from_key", side_effect=branch_from_key), \
             patch.object(sine_gate.KC.PW, "oriented_specs", return_value=[(1, (0, 1, 2, 3))]), \
             patch.object(sine_gate.KC.PW, "T_sequences", return_value=[(1.0 + 0j, ("forward",))]), \
             patch.object(sine_gate.KC.PW, "adjoint_sequence", return_value=("adjoint",)), \
             patch.object(sine_gate.KC.CK, "apply_sequence_to_branch", side_effect=apply_sequence), \
             patch.object(sine_gate.KC.CK, "project_branch_complete_charges", side_effect=project_branch):
            alpha = 0.7 + 0.2j
            beta = -0.3 + 0.9j
            mixed, _, _ = sine_gate.complete_HE_sine_state(
                {key: alpha}, 0, 4, charged_nodes=(0,)
            )
            scaled, _, _ = sine_gate.complete_HE_sine_state(
                {key: beta}, 0, 4, charged_nodes=(0,)
            )

        self.assertEqual(set(mixed), {("output", ("forward",)), ("output", ("adjoint",))})
        for output_key in mixed:
            self.assertAlmostEqual(scaled[output_key], mixed[output_key] * beta / alpha)

    def test_full_gate_rejects_empty_production_operator(self):
        with patch.object(gate.SINE, "complete_HE_sine", return_value=({}, 0.0, 0.0)):
            with self.assertRaisesRegex(AssertionError, "H_E.*M2=-2.*empty"):
                gate.run(Jmax2=10, pruning_mode="production")


if __name__ == "__main__":
    unittest.main(verbosity=2)
