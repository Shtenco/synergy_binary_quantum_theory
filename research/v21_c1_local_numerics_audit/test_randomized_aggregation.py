import math
import random
import unittest


def stable_complex_sum(values):
    real = math.fsum(complex(value).real for value in values)
    imag = math.fsum(complex(value).imag for value in values)
    return complex(real, imag)


def sum_states_exact(states):
    buckets = {}
    for state in states:
        for key, amplitude in state.items():
            buckets.setdefault(key, []).append(amplitude)
    total = {key: stable_complex_sum(values) for key, values in buckets.items()}
    return {key: value for key, value in total.items() if value != 0j}


class RandomizedAggregationAudit(unittest.TestCase):
    def test_seeded_cancellation_sweep_is_order_invariant(self):
        rng = random.Random(20261108)
        states = []
        for _ in range(80):
            state = {}
            for key_index in range(5):
                scale = 1e12 if rng.random() < 0.5 else 1.0
                sign = -1.0 if rng.random() < 0.5 else 1.0
                real = sign * scale + rng.choice((-1.0, 0.0, 1.0))
                imag = rng.uniform(-2.0, 2.0)
                state[((key_index,), (0,))] = complex(real, imag)
            states.append(state)
        expected = sum_states_exact(states)
        for _ in range(1000):
            permuted = list(states)
            rng.shuffle(permuted)
            permuted = [dict(reversed(list(state.items()))) for state in permuted]
            self.assertEqual(sum_states_exact(permuted), expected)


if __name__ == "__main__":
    unittest.main(verbosity=2)
