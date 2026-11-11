#!/usr/bin/env python3
"""L2 characterization of real V10 KKV state maps.

These tests deliberately exercise the production maps on valid Peter-Weyl
covariant basis states.  Expected failures are *characterizations* of the
currently known state-amplitude pruning; they are not acceptance passes.
"""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_covariant_composition_gate as COMP
import peter_weyl_covariant_K_composition_gate as KCOMP
import peter_weyl_zeroaware_volume_migration_experiment as ZVM


def scale_state(state, scalar):
    return {key: scalar * value for key, value in state.items()}


def add_states(left, right):
    out = dict(left)
    for key, value in right.items():
        out[key] = out.get(key, 0j) + value
    return out


def diff_norm(left, right):
    keys = set(left) | set(right)
    return math.sqrt(sum(abs(left.get(k, 0j) - right.get(k, 0j)) ** 2 for k in keys))


def max_amplitude(state):
    return max((abs(value) for value in state.values()), default=0.0)


def complete_covariant_basis_state(gauss_key, source_v=0):
    covariant = CV.gauss_to_covariant({gauss_key: 1.0 + 0j}, source_v)
    return COMP.covariant_to_complete(covariant, source_v)


class TestL2ProductionMapLinearity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Match the zero-aware volume convention used by the frozen V10
        # covariant-K composition experiment. Clear all dependent caches.
        ZVM.patch_and_clear()
        CK.HE_complete_cached.cache_clear()
        CK.canonical_volume_block_general.cache_clear()

        cls.source_v = 0
        cls.initial_gauss_key = PW.basis_full_jhalf()[0]
        cls.input_state = complete_covariant_basis_state(
            cls.initial_gauss_key, cls.source_v
        )
        if len(cls.input_state) != 1:
            raise AssertionError("Expected a single valid complete covariant basis key")
        cls.input_key = next(iter(cls.input_state))

        # The fixture checks are ordinary hard tests.  They prevent an empty
        # or unsupported sector from being misreported as an expected failure.
        cls.volume_unit = CK.apply_volume_complete_state(
            cls.input_state, cls.source_v
        )
        cls.he_unit, cls.he_vleak, cls.he_bleak = CK.apply_HE_complete_state(
            cls.input_state, cls.source_v, 5
        )
        cls.custom_he_unit, cls.custom_vleak, cls.custom_bleak = (
            KCOMP.apply_HE_complete_state_custom(
                cls.input_state, cls.source_v, 5, (cls.source_v,)
            )
        )

    def test_fixture_reaches_nonzero_volume_output(self):
        self.assertTrue(
            self.volume_unit,
            "Chosen valid J=0 covariant fixture has zero volume response; "
            "select a preregistered nonzero basis fixture before interpreting L2.",
        )

    def test_fixture_reaches_nonzero_HE_output(self):
        self.assertTrue(
            self.he_unit,
            "Chosen valid covariant fixture has zero H_E response; "
            "this is a fixture failure, not evidence of linearity.",
        )
        self.assertTrue(math.isfinite(self.he_vleak))
        self.assertTrue(math.isfinite(self.he_bleak))

    def test_fixture_reaches_nonzero_custom_HE_output(self):
        self.assertTrue(
            self.custom_he_unit,
            "Custom H_E production fixture is empty; no L2 conclusion is valid.",
        )
        self.assertTrue(math.isfinite(self.custom_vleak))
        self.assertTrue(math.isfinite(self.custom_bleak))

    @unittest.expectedFailure
    def test_volume_map_additivity_across_output_threshold(self):
        # CK.apply_volume_complete_state removes output entries <= 1e-11.
        # Choose alpha from the actual unit-column norm so F(alpha*x) lies
        # below the threshold while F(2*alpha*x) crosses it.
        peak = max_amplitude(self.volume_unit)
        self.assertGreater(peak, 0.0)
        alpha = 0.75e-11 / peak
        x = scale_state(self.input_state, alpha)
        twice_x = scale_state(self.input_state, 2.0 * alpha)
        fx = CK.apply_volume_complete_state(x, self.source_v)
        f2x = CK.apply_volume_complete_state(twice_x, self.source_v)
        residual = diff_norm(f2x, add_states(fx, fx))
        self.assertLess(
            residual, 1e-13 * max(1.0, max_amplitude(f2x)),
            f"volume additivity residual across 1e-11 output pruning: {residual}",
        )

    @unittest.expectedFailure
    def test_HE_state_map_homogeneity_across_output_threshold(self):
        # CK.apply_HE_complete_state applies a final 1e-10 state-amplitude
        # cutoff after arbitrary input coefficients have been applied.
        peak = max_amplitude(self.he_unit)
        self.assertGreater(peak, 0.0)
        alpha = (0.75e-10 / peak) * (0.7 + 0.4j)
        scaled = CK.apply_HE_complete_state(
            scale_state(self.input_state, alpha), self.source_v, 5
        )[0]
        expected = scale_state(self.he_unit, alpha)
        residual = diff_norm(scaled, expected)
        self.assertLess(
            residual, 1e-13 * max(1.0, max_amplitude(expected)),
            f"H_E state-map homogeneity residual across 1e-10 pruning: {residual}",
        )

    @unittest.expectedFailure
    def test_custom_HE_map_additivity_across_COMP_threshold(self):
        # The active K-composition wrapper applies arbitrary input amplitudes
        # through COMP.add(TOL=1e-11).  The two half-size actions are pruned
        # separately while their coherent sum crosses the cutoff.
        peak = max_amplitude(self.custom_he_unit)
        self.assertGreater(peak, 0.0)
        alpha = 0.75e-11 / peak
        x = scale_state(self.input_state, alpha)
        twice_x = scale_state(self.input_state, 2.0 * alpha)
        fx = KCOMP.apply_HE_complete_state_custom(
            x, self.source_v, 5, (self.source_v,)
        )[0]
        f2x = KCOMP.apply_HE_complete_state_custom(
            twice_x, self.source_v, 5, (self.source_v,)
        )[0]
        residual = diff_norm(f2x, add_states(fx, fx))
        self.assertLess(
            residual, 1e-13 * max(1.0, max_amplitude(f2x)),
            f"custom H_E additivity residual across COMP.TOL: {residual}",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
