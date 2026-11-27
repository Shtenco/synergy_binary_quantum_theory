import math
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from bqg_neumaier_sparse_accumulator import NeumaierSparseAccumulator

EPS = 2.220446049250313e-16


def _fsum_map(items):
    parts = {}
    for key, value in items:
        parts.setdefault(key, []).append(value)
    return {
        key: complex(math.fsum(z.real for z in values), math.fsum(z.imag for z in values))
        for key, values in parts.items()
    }


def _scale_map(items):
    parts = {}
    for key, value in items:
        parts[key] = parts.get(key, 0.0) + abs(value)
    return parts


def _assert_maps_roundoff_bounded(actual, expected, items, factor=32.0):
    scale = _scale_map(items)
    for key in set(actual) | set(expected):
        error = abs(actual.get(key, 0j) - expected.get(key, 0j))
        bound = factor * EPS * scale.get(key, 0.0) + 1e-300
        assert error <= bound, (key, error, bound, actual.get(key), expected.get(key))


def test_seeded_single_pass_and_partition_merge_match_fsum():
    rng = random.Random(20261127)
    for _ in range(3000):
        n = rng.randint(3, 32)
        items = []
        for _ in range(n):
            key = rng.randrange(8)
            mag = 10 ** rng.uniform(-13, 3)
            phase = rng.uniform(-math.pi, math.pi)
            items.append((key, mag * complex(math.cos(phase), math.sin(phase))))
        expected = _fsum_map(items)
        one = NeumaierSparseAccumulator()
        parts = [NeumaierSparseAccumulator() for _ in range(4)]
        for i, (key, value) in enumerate(items):
            one.add(key, value)
            parts[i % len(parts)].add(key, value)
        merged = NeumaierSparseAccumulator()
        for part in reversed(parts):
            merged.merge(part)
        _assert_maps_roundoff_bounded(one.finalize(), expected, items)
        _assert_maps_roundoff_bounded(merged.finalize(), expected, items)


def test_permutation_invariance_on_seeded_cancellation_sets():
    rng = random.Random(17)
    for _ in range(1000):
        items = [(rng.randrange(5), complex(rng.uniform(-1, 1), rng.uniform(-1, 1))) for _ in range(40)]
        a = NeumaierSparseAccumulator()
        b = NeumaierSparseAccumulator()
        for key, value in items:
            a.add(key, value)
        shuffled = list(items)
        rng.shuffle(shuffled)
        for key, value in shuffled:
            b.add(key, value)
        expected = _fsum_map(items)
        _assert_maps_roundoff_bounded(a.finalize(), expected, items)
        _assert_maps_roundoff_bounded(b.finalize(), expected, items)


def test_no_intermediate_magnitude_pruning_under_adversarial_cancellation():
    tau = 1e-11
    expected = 1200 * tau
    for key in range(128):
        acc = NeumaierSparseAccumulator()
        for _ in range(1000):
            acc.add(key, 1.2 * tau)
            acc.add(key, -0.7 * tau)
            acc.add(key, 0.7 * tau)
        result = acc.finalize()
        assert key in result
        assert abs(result[key].real - expected) <= 8 * EPS * 384000 * tau


def test_complex_scaling_and_add_state():
    x = {"a": 1.25 - 2.5j, "b": -0.75 + 0.125j}
    y = {"a": -0.5 + 0.25j, "c": 2.0j}
    scale = -0.3 + 0.7j
    items = [(k, scale * v) for k, v in list(x.items()) + list(y.items())]
    acc = NeumaierSparseAccumulator()
    acc.add_state(x, scale)
    acc.add_state(y, scale)
    _assert_maps_roundoff_bounded(acc.finalize(), _fsum_map(items), items)


def test_extreme_exponent_range_and_self_merge():
    rng = random.Random(2026112701)
    for _ in range(5000):
        items = []
        acc = NeumaierSparseAccumulator()
        for _ in range(rng.randint(4, 24)):
            exponent = rng.uniform(-300, 250)
            magnitude = 10.0 ** exponent
            phase = rng.uniform(-math.pi, math.pi)
            value = magnitude * complex(math.cos(phase), math.sin(phase))
            key = rng.randrange(4)
            items.append((key, value))
            acc.add(key, value)
        expected = _fsum_map(items)
        _assert_maps_roundoff_bounded(acc.finalize(), expected, items, factor=64.0)

    single = NeumaierSparseAccumulator()
    single.add("x", 1.0 + 2.0j)
    single.merge(single)
    assert single.finalize()["x"] == 2.0 + 4.0j


def test_exact_zero_policy_is_finalization_only():
    acc = NeumaierSparseAccumulator()
    acc.add("z", 1.0 + 2.0j)
    acc.add("z", -1.0 - 2.0j)
    assert "z" not in acc.finalize()
    assert acc.finalize(keep_exact_zeros=True)["z"] == 0j
