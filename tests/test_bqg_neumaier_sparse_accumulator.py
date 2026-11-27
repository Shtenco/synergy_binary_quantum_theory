import math
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from bqg_neumaier_sparse_accumulator import NeumaierSparseAccumulator


def _fsum_map(items):
    parts = {}
    for key, value in items:
        parts.setdefault(key, []).append(value)
    return {
        key: complex(math.fsum(z.real for z in values), math.fsum(z.imag for z in values))
        for key, values in parts.items()
    }


def _assert_maps_close(actual, expected, atol=1e-25):
    for key in set(actual) | set(expected):
        assert abs(actual.get(key, 0j) - expected.get(key, 0j)) <= atol


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
        left = NeumaierSparseAccumulator()
        right = NeumaierSparseAccumulator()
        cut = rng.randrange(1, n)
        for i, (key, value) in enumerate(items):
            one.add(key, value)
            (left if i < cut else right).add(key, value)
        left.merge(right)
        _assert_maps_close(one.finalize(), expected, atol=1e-10)
        _assert_maps_close(left.finalize(), expected, atol=1e-10)


def test_permutation_invariance_on_seeded_cancellation_sets():
    rng = random.Random(17)
    for _ in range(1000):
        items = [(rng.randrange(5), complex(rng.uniform(-1, 1), rng.uniform(-1, 1))) for _ in range(40)]
        a = NeumaierSparseAccumulator()
        b = NeumaierSparseAccumulator()
        for key, value in items:
            a.add(key, value)
        rng.shuffle(items)
        for key, value in items:
            b.add(key, value)
        _assert_maps_close(a.finalize(), b.finalize(), atol=1e-13)


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
        assert abs(result[key].real - expected) < 1e-20


def test_complex_scaling_and_add_state():
    x = {"a": 1.25 - 2.5j, "b": -0.75 + 0.125j}
    y = {"a": -0.5 + 0.25j, "c": 2.0j}
    scale = -0.3 + 0.7j
    acc = NeumaierSparseAccumulator()
    acc.add_state(x, scale)
    acc.add_state(y, scale)
    expected = _fsum_map([(k, scale * v) for k, v in list(x.items()) + list(y.items())])
    _assert_maps_close(acc.finalize(), expected, atol=1e-14)
