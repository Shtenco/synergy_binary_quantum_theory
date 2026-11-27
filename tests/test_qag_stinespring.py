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
        self.recoil = gaussian_momentum_state(32, center=0.0, width=1.5)
        self.probe = gaussian_momentum_state(32, center=0.0, width=4.0)

    def test_pre_postselection_unitary_model_normalizes(self):
        result = run_stinespring_recoil_model(
            self.source, self.post, self.recoil, self.probe,
            impulse_units=1, pointer_angle=0.0
        )
        self.assertAlmostEqual(np.trace(result.probe_state).real, 1.0, places=12)
        self.assertGreaterEqual(np.linalg.eigvalsh(result.probe_state).min(), -1e-12)
        self.assertLess(result.recoil_wrap_probability, 1e-12)
        self.assertLess(result.probe_wrap_probability, 1e-3)

    def test_complete_ensemble_preserves_modular_momentum_phase(self):
        outcomes = [
            run_stinespring_recoil_model(
                self.source, np.eye(2, dtype=complex)[:, k],
                self.recoil, self.probe, impulse_units=1, pointer_angle=0.0
            )
            for k in range(2)
        ]
        p_total = sum(item.probability for item in outcomes)
        self.assertAlmostEqual(p_total, 1.0, places=12)
        ensemble_phase = sum(
            item.probability * item.total_momentum_phase_conditional
            for item in outcomes
        )
        reference_phase = outcomes[0].total_momentum_phase_initial
        self.assertAlmostEqual(ensemble_phase.real, reference_phase.real, places=12)
        self.assertAlmostEqual(ensemble_phase.imag, reference_phase.imag, places=12)

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
        self.assertGreater(
            np.linalg.norm(recorded.probe_state - zero.probe_state), 1e-10
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
        self.assertEqual(result.recoil_wrap_probability, 0.0)
        self.assertEqual(result.probe_wrap_probability, 0.0)

    def test_wrap_diagnostic_detects_boundary_crossing(self):
        n = 32
        near_upper_edge = gaussian_momentum_state(n, center=15.0, width=0.8)
        near_lower_edge = gaussian_momentum_state(n, center=-15.0, width=0.8)
        result = run_stinespring_recoil_model(
            np.array([0.0, 1.0]), self.post,
            near_lower_edge, near_upper_edge,
            impulse_units=1, pointer_angle=0.0
        )
        self.assertGreater(result.recoil_wrap_probability, 0.1)
        self.assertGreater(result.probe_wrap_probability, 0.1)

    def test_rejects_invalid_impulse(self):
        with self.assertRaises(ValueError):
            run_stinespring_recoil_model(
                self.source, self.post, self.recoil, self.probe,
                impulse_units=16
            )

    def test_rejects_nonfinite_parameters(self):
        with self.assertRaises(ValueError):
            gaussian_momentum_state(32, center=np.nan)
        with self.assertRaises(ValueError):
            run_stinespring_recoil_model(
                self.source, self.post, self.recoil, self.probe,
                pointer_angle=np.inf
            )


if __name__ == "__main__":
    unittest.main()
