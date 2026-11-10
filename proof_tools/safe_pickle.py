"""Restricted pickle loading for numeric proof artifacts.

Pickle is executable by default. These helpers allow only inert Python
containers/scalars, and optionally NumPy numeric arrays/scalars. Object arrays
and every unlisted global are rejected before the object is returned.
"""
from __future__ import annotations

import builtins
import importlib
import io
import pickle
from typing import BinaryIO

import numpy as np


_SAFE_BUILTINS = {
    "bool", "bytearray", "bytes", "complex", "dict", "float", "frozenset",
    "int", "list", "set", "slice", "str", "tuple",
}
_SAFE_NUMPY_GLOBALS = {
    ("numpy", "ndarray"),
    ("numpy", "dtype"),
    ("numpy.core.multiarray", "_reconstruct"),
    ("numpy.core.multiarray", "scalar"),
    ("numpy.core.numeric", "_frombuffer"),
    ("numpy._core.multiarray", "_reconstruct"),
    ("numpy._core.multiarray", "scalar"),
    ("numpy._core.numeric", "_frombuffer"),
}
_NUMERIC_KINDS = frozenset("biufc")
_MAX_CONTAINER_DEPTH = 128
_MAX_CONTAINER_NODES = 10_000_000


class RestrictedUnpickler(pickle.Unpickler):
    def __init__(self, file: BinaryIO, *, allow_numpy: bool = False):
        super().__init__(file)
        self.allow_numpy = allow_numpy

    def find_class(self, module: str, name: str):
        if module == "builtins" and name in _SAFE_BUILTINS:
            return getattr(builtins, name)
        if self.allow_numpy and (module, name) in _SAFE_NUMPY_GLOBALS:
            imported = importlib.import_module(module)
            return getattr(imported, name)
        raise pickle.UnpicklingError(
            f"forbidden pickle global: {module}.{name}"
        )


def _validate_inert_tree(root, *, allow_numpy: bool) -> None:
    stack = [(root, 0)]
    seen = set()
    visited = 0
    while stack:
        value, depth = stack.pop()
        visited += 1
        if visited > _MAX_CONTAINER_NODES:
            raise pickle.UnpicklingError("pickle object graph exceeds node limit")
        if depth > _MAX_CONTAINER_DEPTH:
            raise pickle.UnpicklingError("pickle object graph exceeds depth limit")

        value_type = type(value)
        if value is None or value_type in (bool, int, float, complex, str, bytes, bytearray):
            continue

        if allow_numpy and value_type is np.ndarray:
            if value.dtype.hasobject or value.dtype.kind not in _NUMERIC_KINDS:
                raise pickle.UnpicklingError(
                    f"non-numeric or object NumPy dtype is forbidden: {value.dtype}"
                )
            continue

        if allow_numpy and isinstance(value, np.generic):
            if value.dtype.hasobject or value.dtype.kind not in _NUMERIC_KINDS:
                raise pickle.UnpicklingError(
                    f"non-numeric NumPy scalar is forbidden: {value.dtype}"
                )
            continue

        if value_type in (list, tuple, set, frozenset):
            identity = id(value)
            if identity in seen:
                continue
            seen.add(identity)
            stack.extend((item, depth + 1) for item in value)
            continue

        if value_type is dict:
            identity = id(value)
            if identity in seen:
                continue
            seen.add(identity)
            for key, item in value.items():
                stack.append((key, depth + 1))
                stack.append((item, depth + 1))
            continue

        raise pickle.UnpicklingError(
            f"forbidden object type in pickle payload: {value_type.__module__}.{value_type.__name__}"
        )


def safe_load(file: BinaryIO, *, allow_numpy: bool = False):
    value = RestrictedUnpickler(file, allow_numpy=allow_numpy).load()
    _validate_inert_tree(value, allow_numpy=allow_numpy)
    return value


def safe_load_bytes(payload: bytes, *, allow_numpy: bool = False):
    return safe_load(io.BytesIO(payload), allow_numpy=allow_numpy)


def safe_load_builtin_pickle(file: BinaryIO):
    return safe_load(file, allow_numpy=False)


def safe_load_numpy_pickle(file: BinaryIO):
    return safe_load(file, allow_numpy=True)
