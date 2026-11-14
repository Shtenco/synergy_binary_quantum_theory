#!/usr/bin/env python3
"""L2 tests against real V10 production primitives.

The fixture is a genuine finite Peter-Weyl basis state converted into the
production covariant/complete representation.  Volume, cached H_E, and custom
composition H_E threshold-boundary characterizations are hard acceptance
checks: all three maps must preserve additivity at the tested boundary.
Safe-domain complex homogeneity remains a hard acceptance check.
"""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_covariant_K_composition_gate as KCOMP
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_covariant_composition_gate as COMP
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_zeroaware_volume_migration_experiment as ZVM


def state_scale(state, scale):
    return {key: scale * value for key, value in state.items()}


def output_dict(result):
    # Production functions return (state, diagnostics...) except the volume leg.
    return result if isinstance(result, dict) else result[0]


def norm(state):
    return math.sqrt(sum(abs(value) ** 2 for value in state.values()))


def difference_norm(left, right):
    keys = set(left) | set(right)
    return math.sqrt(sum(abs(left.get(k, 0j) - right.get(k, 0j)) ** 2 for k in keys))


def scaled_map(map_fn, state, scale):
    return output_dict(map_fn(state_scale(state, scale)))


class TestV10ProductionPrimitiveLinearity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Match the zero-aware volume convention used by the V10 composition
        # research path.  Clear affected caches before constructing the fixture.
        ZVM.patch_and_clear()
        initial = PW.basis_full_jhalf()[0]
        covariant = CV.gauss_to_covariant({initial: 1.0 + 0j}, 0)
        cls.complete_state = COMP.covariant_to_complete(covariant, 0)
        if len(cls.complete_state) != 1:
            raise RuntimeError(
                "Expected one frozen source basis key; got "
                f"{len(cls.complete_state)}"
            )

    def test_fixture_volume_map_is_nonzero(self):
        result = CK.apply_volume_complete_state(self.complete_state, 0)
        self.assertTrue(
            result,
            "Frozen Peter-Weyl fixture has no volume output; L2 cannot be inferred.",
        )

    def test_complete_branch_projection_is_homogeneous_below_amplitude_cutoff(self):
        # The tensor/basis projection coefficients are fixed by the source
        # basis key, while the branch amplitude is an arbitrary complex input.
        # Scaling that input below the historical 1e-11 threshold must scale
        # the complete projected state rather than erase it.
        key = next(iter(self.complete_state))
        alpha = 2.5e-14 - 1.75e-14j
        unit, unit_leak = CK.project_branch_complete_charges(
            CK.branch_from_key(key, 1.0 + 0j)
        )
        scaled, scaled_leak = CK.project_branch_complete_charges(
            CK.branch_from_key(key, alpha)
        )
        self.assertTrue(unit, "Reference complete branch projection is empty.")
        self.assertTrue(
            scaled,
            "A nonzero input-scaled branch projection was removed by a fixed amplitude cutoff.",
        )
        expected = {state: alpha * value for state, value in unit.items()}
        residual = difference_norm(scaled, expected)
        relative = residual / max(norm(scaled), norm(expected), 1e-300)
        self.assertLessEqual(
            relative, 5e-12,
            f"complete branch projection complex homogeneity residual={relative:.6e}",
        )
        self.assertLessEqual(abs(unit_leak - scaled_leak), 5e-14)

    def test_fixture_HE_state_map_is_nonzero(self):
        result = output_dict(CK.apply_HE_complete_state(self.complete_state, 0, 5))
        self.assertTrue(
            result,
            "Frozen Peter-Weyl fixture has no H_E output; L2 cannot be inferred.",
        )

    def test_fixture_custom_HE_map_is_nonzero(self):
        result = output_dict(
            KCOMP.apply_HE_complete_state_custom(self.complete_state, 0, 5, (0,))
        )
        self.assertTrue(
            result,
            "Frozen Peter-Weyl fixture has no custom H_E output; L2 cannot be inferred.",
        )

    def test_volume_map_additivity_on_independent_basis_states(self):
        # Use two distinct genuine Peter-Weyl basis keys, not two rescalings
        # of the same vector.  The larger common scale keeps the comparison
        # away from the declared 1e-11 output-amplitude cutoff.
        basis = PW.basis_full_jhalf()
        second_covariant = CV.gauss_to_covariant({basis[16]: 1.0 + 0j}, 0)
        second_complete = COMP.covariant_to_complete(second_covariant, 0)
        self.assertEqual(len(second_complete), 1)
        self.assertNotEqual(next(iter(second_complete)), next(iter(self.complete_state)))

        x = state_scale(self.complete_state, 1e6)
        y = state_scale(second_complete, 1e6)
        fx = CK.apply_volume_complete_state(x, 0)
        fy = CK.apply_volume_complete_state(y, 0)
        self.assertTrue(fx, "First independent basis vector has no volume response.")
        self.assertTrue(fy, "Second independent basis vector has no volume response.")
        combined = dict(x)
        for key, value in y.items():
            combined[key] = combined.get(key, 0j) + value
        fxy = CK.apply_volume_complete_state(combined, 0)

        rhs = dict(fx)
        for key, value in fy.items():
            rhs[key] = rhs.get(key, 0j) + value
        residual = difference_norm(fxy, rhs)
        relative = residual / max(norm(fxy), norm(rhs), 1e-300)
        self.assertLessEqual(
            relative, 5e-12,
            f"CK.apply_volume_complete_state independent-superposition "
            f"additivity residual={relative:.6e}",
        )

    def test_he_state_map_additivity_on_independent_basis_states(self):
        # Test the cached-basis-column extension on two independent
        # Peter-Weyl inputs, not merely homogeneity on one input ray.
        basis = PW.basis_full_jhalf()
        second_covariant = CV.gauss_to_covariant({basis[16]: 1.0 + 0j}, 0)
        second_complete = COMP.covariant_to_complete(second_covariant, 0)
        self.assertEqual(len(second_complete), 1)
        self.assertNotEqual(
            next(iter(second_complete)), next(iter(self.complete_state))
        )

        # CK.apply_HE_complete_state applies cached unit-key columns and
        # combines their input coefficients before its final 1e-10 filter.
        # A large common scale keeps the comparison away from that final
        # output threshold; this does not test threshold-boundary behavior.
        x = state_scale(self.complete_state, 1e6)
        y = state_scale(second_complete, 1e6)
        fx = output_dict(CK.apply_HE_complete_state(x, 0, 5))
        fy = output_dict(CK.apply_HE_complete_state(y, 0, 5))
        self.assertTrue(fx, "First independent basis vector has no H_E response.")
        self.assertTrue(fy, "Second independent basis vector has no H_E response.")

        combined = dict(x)
        for key, value in y.items():
            combined[key] = combined.get(key, 0j) + value
        fxy = output_dict(CK.apply_HE_complete_state(combined, 0, 5))

        rhs = dict(fx)
        for key, value in fy.items():
            rhs[key] = rhs.get(key, 0j) + value
        residual = difference_norm(fxy, rhs)
        relative = residual / max(norm(fxy), norm(rhs), 1e-300)
        self.assertLessEqual(
            relative, 5e-12,
            f"CK.apply_HE_complete_state independent-superposition "
            f"additivity residual={relative:.6e}",
        )

    def test_full_covariant_C_volume_component_additivity_on_independent_states(self):
        # Exercise the complete covariant C_e(V) composition, not only its
        # local volume primitive.  Inputs are two independent Gauss basis
        # states mapped into the production covariant representation.
        basis = PW.basis_full_jhalf()
        cov_x = CV.gauss_to_covariant({basis[0]: 1.0 + 0j}, 0)
        cov_y = CV.gauss_to_covariant({basis[16]: 1.0 + 0j}, 0)
        self.assertEqual(len(cov_x), 1)
        self.assertEqual(len(cov_y), 1)
        self.assertNotEqual(next(iter(cov_x)), next(iter(cov_y)))

        x = state_scale(cov_x, 1e6)
        y = state_scale(cov_y, 1e6)
        fx = output_dict(COMP.C_volume_component(x, 0, 1, 0, 0, 5))
        fy = output_dict(COMP.C_volume_component(y, 0, 1, 0, 0, 5))
        self.assertTrue(fx, "First independent covariant input has no C(V) response.")
        self.assertTrue(fy, "Second independent covariant input has no C(V) response.")

        combined = dict(x)
        for key, value in y.items():
            combined[key] = combined.get(key, 0j) + value
        fxy = output_dict(COMP.C_volume_component(combined, 0, 1, 0, 0, 5))

        rhs = dict(fx)
        for key, value in fy.items():
            rhs[key] = rhs.get(key, 0j) + value
        residual = difference_norm(fxy, rhs)
        relative = residual / max(norm(fxy), norm(rhs), 1e-300)
        self.assertLessEqual(
            relative, 5e-12,
            f"COMP.C_volume_component independent-superposition "
            f"additivity residual={relative:.6e}",
        )

    def _assert_safe_domain_homogeneity(self, map_fn, label):
        alpha = 1.25 + 1.5j
        base_scale = 1e6
        base = output_dict(map_fn(state_scale(self.complete_state, base_scale)))
        if not base:
            self.skipTest(f"{label}: production primitive returned empty reference output")
        transformed = output_dict(
            map_fn(state_scale(self.complete_state, alpha * base_scale))
        )
        expected = {key: alpha * value for key, value in base.items()}
        residual = difference_norm(transformed, expected)
        relative = residual / max(norm(transformed), norm(expected), 1e-300)
        self.assertLessEqual(relative, 5e-12, f"{label}: complex homogeneity residual={relative:.6e}")

    def _assert_threshold_additivity_characterization(self, map_fn, threshold, label):
        unit = output_dict(map_fn(self.complete_state))
        if not unit:
            self.skipTest(f"{label}: production primitive returned empty reference output")
        peak = max(abs(value) for value in unit.values())
        scale = 0.6 * threshold / peak
        x = state_scale(self.complete_state, scale)
        y = state_scale(self.complete_state, scale)
        lhs = output_dict(map_fn(state_scale(self.complete_state, 2.0 * scale)))
        fx = output_dict(map_fn(x))
        fy = output_dict(map_fn(y))
        rhs = dict(fx)
        for key, value in fy.items():
            rhs[key] = rhs.get(key, 0j) + value
        residual = difference_norm(lhs, rhs)
        relative = residual / max(norm(lhs), norm(rhs), 1e-300)
        self.assertLessEqual(
            relative, 1e-12,
            f"{label}: expected threshold nonlinearity characterization; "
            f"relative additivity residual={relative:.6e}",
        )

    def test_volume_safe_domain_complex_homogeneity(self):
        self._assert_safe_domain_homogeneity(
            lambda state: CK.apply_volume_complete_state(state, 0),
            "CK.apply_volume_complete_state",
        )

    def test_he_state_safe_domain_complex_homogeneity(self):
        self._assert_safe_domain_homogeneity(
            lambda state: CK.apply_HE_complete_state(state, 0, 5),
            "CK.apply_HE_complete_state",
        )

    def test_composition_he_safe_domain_complex_homogeneity(self):
        self._assert_safe_domain_homogeneity(
            lambda state: KCOMP.apply_HE_complete_state_custom(state, 0, 5, (0,)),
            "COMP.apply_HE_complete_state_custom",
        )

    def test_volume_threshold_additivity_characterization(self):
        self._assert_threshold_additivity_characterization(
            lambda state: CK.apply_volume_complete_state(state, 0),
            1e-11,
            "CK.apply_volume_complete_state",
        )

    def test_he_state_threshold_additivity_characterization(self):
        self._assert_threshold_additivity_characterization(
            lambda state: CK.apply_HE_complete_state(state, 0, 5),
            1e-10,
            "CK.apply_HE_complete_state",
        )

    def test_composition_he_threshold_additivity_characterization(self):
        self._assert_threshold_additivity_characterization(
            lambda state: KCOMP.apply_HE_complete_state_custom(state, 0, 5, (0,)),
            COMP.TOL,
            "COMP.apply_HE_complete_state_custom",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
