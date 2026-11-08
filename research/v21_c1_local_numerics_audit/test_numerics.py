import math
import unittest
from collections import defaultdict


def stable_complex_sum(values):
    return complex(
        math.fsum(float(complex(v).real) for v in values),
        math.fsum(float(complex(v).imag) for v in values),
    )


def stable_l2(values):
    norm = 0.0
    for value in values:
        norm = math.hypot(norm, abs(value))
    return norm


def sum_states_exact(states):
    buckets = defaultdict(list)
    for state in states:
        for key, amplitude in state.items():
            buckets[key].append(amplitude)
    result = {}
    for key, values in buckets.items():
        total = stable_complex_sum(values)
        if total != 0j:
            result[key] = total
    return result


def sum_states_legacy_pruned(states, threshold):
    total = {}
    for state in states:
        for key, amplitude in state.items():
            value = total.get(key, 0j) + amplitude
            if abs(value) > threshold:
                total[key] = value
            elif key in total:
                del total[key]
    return total


def compare_states(a, b):
    delta = stable_l2(a.get(k, 0j) - b.get(k, 0j) for k in set(a) | set(b))
    norm_a, norm_b = stable_l2(a.values()), stable_l2(b.values())
    denominator = max(norm_a, norm_b)
    if not all(math.isfinite(v) for v in (delta, norm_a, norm_b)):
        relative = math.inf
    else:
        relative = delta / denominator if denominator > 0 else (0.0 if delta == 0 else math.inf)
    return relative, delta, norm_a, norm_b


class LocalNumericsAudit(unittest.TestCase):
    def test_real_cancellation_order_independence(self):
        values = [1e16, 1.0, -1e16]
        self.assertEqual(stable_complex_sum(values), 1 + 0j)
        self.assertEqual(stable_complex_sum(list(reversed(values))), 1 + 0j)

    def test_complex_channels_are_summed_separately(self):
        values = [1e16 + 1e16j, 1 - 2j, -1e16 + 1e16j]
        self.assertEqual(stable_complex_sum(values), 1 + 2e16j)
        self.assertEqual(stable_complex_sum(list(reversed(values))), 1 + 2e16j)

    def test_underflow_safe_norm(self):
        self.assertGreater(stable_l2([1e-300 + 0j]), 0.0)

    def test_overflow_safe_norm(self):
        norm = stable_l2([1e308 + 0j, 1e308 + 0j])
        self.assertTrue(math.isfinite(norm))
        self.assertGreater(norm, 1e308)

    def test_relative_error_tiny_scale_is_not_hidden(self):
        relative, delta, *_ = compare_states({"x": 1e-200 + 0j}, {"x": 2e-200 + 0j})
        self.assertAlmostEqual(relative, 0.5)
        self.assertGreater(delta, 0.0)

    def test_exact_zero_and_nonzero_are_distinct(self):
        self.assertEqual(sum_states_exact([{"x": 1 + 2j}, {"x": -1 - 2j}]), {})
        self.assertNotEqual(sum_states_exact([{"x": 1e-300 + 0j}]), {})

    def test_legacy_incremental_pruning_is_order_sensitive(self):
        sequence = [-1.5, -1.5, -1.0]
        forward = sum_states_legacy_pruned([{"x": x} for x in sequence], 1.0)
        reverse = sum_states_legacy_pruned([{"x": x} for x in reversed(sequence)], 1.0)
        self.assertNotEqual(forward, reverse)

    def test_exact_lane_preserves_subthreshold_residual(self):
        tolerance = 1e-11
        result = sum_states_exact([{"x": 0.4 * tolerance}, {"x": 0.4 * tolerance}])
        self.assertAlmostEqual(result["x"].real, 0.8 * tolerance)
        self.assertNotEqual(result, {})

    def test_nonfinite_comparison_fails_closed(self):
        relative, *_ = compare_states({"x": complex(math.inf, 0)}, {})
        self.assertEqual(relative, math.inf)

    def test_empty_states_compare_as_zero(self):
        self.assertEqual(compare_states({}, {}), (0.0, 0.0, 0.0, 0.0))

    def test_complex_l2_uses_modulus(self):
        self.assertAlmostEqual(stable_l2([3 + 4j]), 5.0)

    def test_l2_is_order_independent_for_finite_values(self):
        values = [3e200, 4e200, 1e-200, 7e100]
        self.assertEqual(stable_l2(values), stable_l2(list(reversed(values))))

    def test_exact_state_aggregation_is_order_independent(self):
        states = [{"a": 1e16 + 1j}, {"a": 1 + 2j}, {"a": -1e16 + 3j}, {"b": 4j}]
        forward = sum_states_exact(states)
        reverse = sum_states_exact(list(reversed(states)))
        self.assertEqual(forward, reverse)
        self.assertEqual(forward["a"], 1 + 6j)

    def test_legacy_pruning_uses_strict_greater_than(self):
        result = sum_states_legacy_pruned([{"x": 1.0}], 1.0)
        self.assertEqual(result, {})

    def test_relative_error_of_two_zero_states_is_zero(self):
        self.assertEqual(compare_states({}, {"x": 0j})[0], 0.0)

    def test_nan_comparison_fails_closed(self):
        relative, *_ = compare_states({"x": complex(math.nan, 0)}, {})
        self.assertEqual(relative, math.inf)

    def test_state_union_includes_missing_keys(self):
        relative, delta, norm_a, norm_b = compare_states({"left": 3 + 0j}, {"right": 4 + 0j})
        self.assertAlmostEqual(delta, 5.0)
        self.assertAlmostEqual(norm_a, 3.0)
        self.assertAlmostEqual(norm_b, 4.0)
        self.assertAlmostEqual(relative, 1.25)


if __name__ == "__main__":
    unittest.main(verbosity=2)
