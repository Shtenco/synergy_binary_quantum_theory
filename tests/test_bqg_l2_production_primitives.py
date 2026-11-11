#!/usr/bin/env python3
"""L2 tests against real V10 production primitives.

The fixture is a genuine finite Peter-Weyl basis state converted into the
production covariant/complete representation.  The threshold-boundary tests
are expected failures until state-dependent output pruning is removed.  The
safe-domain complex-homogeneity tests are hard acceptance checks.
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

    def _assert_safe_domain_homogeneity(self, map_fn, label):
        alpha = 0.375 + 0.625j
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

    @unittest.expectedFailure
    def test_volume_threshold_additivity_characterization(self):
        self._assert_threshold_additivity_characterization(
            lambda state: CK.apply_volume_complete_state(state, 0),
            1e-11,
            "CK.apply_volume_complete_state",
        )

    @unittest.expectedFailure
    def test_he_state_threshold_additivity_characterization(self):
        self._assert_threshold_additivity_characterization(
            lambda state: CK.apply_HE_complete_state(state, 0, 5),
            1e-10,
            "CK.apply_HE_complete_state",
        )

    @unittest.expectedFailure
    def test_composition_he_threshold_additivity_characterization(self):
        self._assert_threshold_additivity_characterization(
            lambda state: KCOMP.apply_HE_complete_state_custom(state, 0, 5, (0,)),
            COMP.TOL,
            "COMP.apply_HE_complete_state_custom",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
