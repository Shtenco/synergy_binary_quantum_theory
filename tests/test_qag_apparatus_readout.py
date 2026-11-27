import unittest
import numpy as np

from qag_apparatus_readout import (
    apparatus_readout_probabilities,
    apparatus_readout_unitary,
    apply_apparatus_readout_dilation,
    conditional_apparatus_readout,
)


class QAGApparatusReadoutTests(unittest.TestCase):
    def test_controlled_addition_is_unitary(self):
        unitary = apparatus_readout_unitary(7)
        np.testing.assert_allclose(unitary.conj().T @ unitary, np.eye(49), atol=1e-12)

    def test_absolute_position_readout_exposes_reference_frame(self):
        n = 6
        unitary = apparatus_readout_unitary(n)
        shift = np.roll(np.eye(n, dtype=complex), 1, axis=0)
        joint_shift = np.kron(shift, shift)
        # The meter's initialized |0> defines an origin; copying absolute
        # apparatus position is not covariant under a joint shift.
        self.assertGreater(np.linalg.norm(unitary @ joint_shift - joint_shift @ unitary), 1e-3)

    def test_dilation_matches_controlled_addition_unitary(self):
        n = 4
        rng = np.random.default_rng(7)
        state = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        state /= np.linalg.norm(state)
        actual = apply_apparatus_readout_dilation(state)
        unitary = apparatus_readout_unitary(n)
        expected = np.zeros((n, n, n), dtype=complex)
        for source_position in range(n):
            initial = np.zeros((n, n), dtype=complex)
            initial[:, 0] = state[source_position, :]
            expected[source_position] = (unitary @ initial.reshape(-1)).reshape(n, n)
        np.testing.assert_allclose(actual, expected, atol=1e-12)

    def test_meter_records_apparatus_position(self):
        state = np.zeros((5, 5), dtype=complex)
        state[1, 3] = 1.0
        dilated = apply_apparatus_readout_dilation(state)
        self.assertAlmostEqual(np.linalg.norm(dilated[:, :, 3]), 1.0, places=12)
        self.assertAlmostEqual(np.linalg.norm(dilated[:, :, 1]), 0.0, places=12)
        self.assertAlmostEqual(dilated[1, 3, 3], 1.0, places=12)

    def test_all_meter_outcomes_follow_born_rule(self):
        state = np.array([[1, 1j, 0], [1, -1j, 1], [0, 1, 2]], dtype=complex)
        state /= np.linalg.norm(state)
        dilated = apply_apparatus_readout_dilation(state)
        probs = apparatus_readout_probabilities(dilated)
        expected = tuple(float(np.sum(np.abs(state[:, k]) ** 2)) for k in range(3))
        np.testing.assert_allclose(probs, expected, atol=1e-12)
        self.assertAlmostEqual(sum(probs), 1.0, places=12)

    def test_conditional_state_matches_position_projection(self):
        state = np.array([[1, 1j], [2, -1]], dtype=complex)
        state /= np.linalg.norm(state)
        dilated = apply_apparatus_readout_dilation(state)
        probability, conditional = conditional_apparatus_readout(dilated, 1)
        expected = np.zeros_like(state)
        expected[:, 1] = state[:, 1]
        expected /= np.linalg.norm(expected)
        self.assertAlmostEqual(probability, float(np.sum(np.abs(state[:, 1]) ** 2)), places=12)
        np.testing.assert_allclose(conditional, expected, atol=1e-12)

    def test_zero_probability_outcome_is_explicit(self):
        state = np.zeros((4, 4), dtype=complex)
        state[0, 0] = 1.0
        dilated = apply_apparatus_readout_dilation(state)
        with self.assertRaises(ZeroDivisionError):
            conditional_apparatus_readout(dilated, 2)

    def test_rejects_invalid_states_and_outcomes(self):
        with self.assertRaises(ValueError):
            apply_apparatus_readout_dilation(np.zeros((2, 3)))
        state = np.zeros((3, 3), dtype=complex)
        state[0, 0] = 1.0
        dilated = apply_apparatus_readout_dilation(state)
        with self.assertRaises(ValueError):
            conditional_apparatus_readout(dilated, 3)


if __name__ == "__main__":
    unittest.main()
