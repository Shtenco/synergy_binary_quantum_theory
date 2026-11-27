import unittest
import numpy as np

from qag_stinespring import (
    complete_stinespring_probabilities,
    gaussian_momentum_state,
    run_stinespring_recoil_model,
)


class QAGStinespringTests(unittest.TestCase):
    def setUp(self):
        eps = 0.1
        self.source = np.array([np.sqrt(1 + eps), np.sqrt(1 - eps)]) / np.sqrt(2)
        self.post = np.array([1.0, -1.0]) / np.sqrt(2)
        self.recoil = gaussian_momentum_state(32, center=0.0, width=4.0)
        self.probe = gaussian_momentum_state(32, center=0.0, width=4.0)

    def test_pre_postselection_unitary_model_normalizes(self):
        result = run_stinespring_recoil_model(
            self.source, self.post, self.recoil, self.probe,
            impulse_units=1, pointer_angle=0.0
        )
        self.assertAlmostEqual(np.trace(result.probe_state).real, 1.0, places=12)
        self.assertGreaterEqual(np.linalg.eigvalsh(result.probe_state).min(), -1e-12)

    def test_recoil_and_probe_shifts_cancel_in_full_conditional_state(self):
        result = run_stinespring_recoil_model(
            self.source, self.post, self.recoil, self.probe,
            impulse_units=1, pointer_angle=0.0
        )
        self.assertAlmostEqual(
            result.total_momentum_conditional,
            result.total_momentum_initial, places=10
        )
        self.assertAlmostEqual(
            result.source_recoil_shift + result.probe_shift, 0.0, places=10
        )

    def test_pointer_overlap_controls_interference(self):
        zero = run_stinespring_recoil_model(
            self.source, self.post, self.recoil, self.probe,
            impulse_units=1, pointer_angle=0.0
        )
        angle = 0.8
        recorded = run_stinespring_recoil_model(
            self.source, self.post, self.recoil, self.probe,
            impulse_units=1, pointer_angle=angle
        )
        self.assertAlmostEqual(
            recorded.pointer_branch_overlap.real, np.cos(angle / 2), places=12
        )
        self.assertLess(
            np.linalg.norm(recorded.probe_state - zero.probe_state), 0.1
        )

    def test_all_source_measurement_outcomes_sum_to_one(self):
        p0, p1 = complete_stinespring_probabilities(
            self.source, np.eye(2), self.recoil, self.probe,
            impulse_units=1, pointer_angle=0.6
        )
        self.assertAlmostEqual(p0 + p1, 1.0, places=12)

    def test_zero_impulse_does_not_transfer_momentum(self):
        result = run_stinespring_recoil_model(
            self.source, self.post, self.recoil, self.probe,
            impulse_units=0, pointer_angle=0.3
        )
        self.assertAlmostEqual(result.source_recoil_shift, 0.0, places=12)
        self.assertAlmostEqual(result.probe_shift, 0.0, places=12)

    def test_rejects_invalid_impulse(self):
        with self.assertRaises(ValueError):
            run_stinespring_recoil_model(
                self.source, self.post, self.recoil, self.probe,
                impulse_units=16
            )


if __name__ == "__main__":
    unittest.main()
