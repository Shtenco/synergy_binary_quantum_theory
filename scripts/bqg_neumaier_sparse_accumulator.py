#!/usr/bin/env python3
"""Streaming, threshold-free complex sparse accumulator candidate for BQG.

Research candidate only: not integrated into the production operator. Uses
Neumaier compensation independently for real and imaginary components.
No magnitude pruning occurs during accumulation.
"""
from __future__ import annotations


class NeumaierSparseAccumulator:
    __slots__ = ("_data",)

    def __init__(self):
        # key -> [sum_re, correction_re, sum_im, correction_im]
        self._data = {}

    @staticmethod
    def _add_scalar(s: float, c: float, x: float) -> tuple[float, float]:
        t = s + x
        if abs(s) >= abs(x):
            c += (s - t) + x
        else:
            c += (x - t) + s
        return t, c

    def add(self, key, amplitude: complex, scale: complex = 1.0) -> None:
        z = scale * amplitude
        v = self._data.get(key)
        if v is None:
            v = [0.0, 0.0, 0.0, 0.0]
            self._data[key] = v
        v[0], v[1] = self._add_scalar(v[0], v[1], float(z.real))
        v[2], v[3] = self._add_scalar(v[2], v[3], float(z.imag))

    def add_state(self, state, scale: complex = 1.0) -> None:
        for key, amplitude in state.items():
            self.add(key, amplitude, scale)

    def merge(self, other: "NeumaierSparseAccumulator") -> None:
        if not isinstance(other, NeumaierSparseAccumulator):
            raise TypeError("other must be a NeumaierSparseAccumulator")
        # Snapshot is required for self-merge because each key's correction
        # fields are mutated by add().
        items = [(key, tuple(v)) for key, v in other._data.items()]
        for key, v in items:
            self.add(key, complex(v[0], v[2]))
            self.add(key, complex(v[1], v[3]))

    def finalize(self, *, keep_exact_zeros: bool = False) -> dict:
        out = {}
        for key, v in self._data.items():
            z = complex(v[0] + v[1], v[2] + v[3])
            if keep_exact_zeros or z != 0j:
                out[key] = z
        return out
