#!/usr/bin/env python3
"""Pure V21 Lorentzian workflow dependency/census checks; no SU(2) physics."""
from itertools import product
import re


def operations(targets, mode, aux):
    a, b, c = targets
    i, j, k = aux
    if mode == "forward":
        return [("CV", c, k, i), ("CK", b, j, k), ("CK", a, i, j)]
    if mode == "adjoint":
        return [("CK", a, j, i), ("CK", b, k, j), ("CV", c, i, k)]
    raise ValueError("INVALID_MODE")


def all_32_specs(abc):
    a, b, c = abc
    for mode, order in product(("forward", "adjoint"), ("p", "q")):
        targets = (a, b, c) if order == "p" else (b, a, c)
        for aux in product((0, 1), repeat=3):
            yield {"mode": mode, "order": order, "aux": list(aux),
                   "targets": targets, "operations": operations(targets, mode, aux)}


def first_leg_census(abc):
    specs = list(all_32_specs(abc))
    first_ops = {tuple(row["operations"][0]) for row in specs}
    pairs = {(row["mode"], row["order"], tuple(row["aux"])) for row in specs}
    if len(pairs) != 32 or len(first_ops) != 12:
        raise ValueError("INVALID_32_PATH_CENSUS")
    return len(first_ops), len(pairs)


def _sha(value):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def predecessor_ok(prev, cur, stage, first_operation, second_operation=None):
    if not isinstance(prev, dict) or not isinstance(cur, dict):
        raise ValueError("BAD_STAGE_METADATA")
    if stage not in (2, 3):
        raise ValueError("INVALID_STAGE")
    if (not isinstance(stage, int) or isinstance(stage, bool)
            or not isinstance(prev.get("stage"), int)
            or isinstance(prev.get("stage"), bool)
            or prev.get("stage") != stage - 1 or prev.get("passed") is not True
            or prev.get("status") != "REAL_ORIGINAL_STAGE_ONLY"):
        raise ValueError("STAGE_NOT_PASSED")
    previous_guard = prev.get("guard")
    if (not isinstance(previous_guard, dict)
            or previous_guard.get("stage_physical_basis_pass") is not True):
        raise ValueError("PREDECESSOR_GUARD_NOT_PASSED")
    for key in ("source_commit", "source_node", "omit", "cycle", "mode",
                "Jmax2", "volume_backend", "volume_backend_sha256",
                "seed_state_sha256"):
        if prev.get(key) != cur.get(key):
            raise ValueError("MISMATCH_" + key)
    if not (_sha(prev.get("seed_state_sha256"))
            and _sha(cur.get("seed_state_sha256"))):
        raise ValueError("BAD_SEED_STATE_DIGEST")
    previous_output = prev.get("output_state_sha256")
    previous_npz = prev.get("npz_sha256")
    current_input = cur.get("input_state_sha256")
    if not (_sha(previous_output) and _sha(previous_npz) and _sha(current_input)):
        raise ValueError("MISSING_CHECKPOINT_CHAIN_DIGEST")
    if previous_output != previous_npz:
        raise ValueError("PREDECESSOR_OUTPUT_DIGEST_MISMATCH")
    if current_input != previous_output:
        raise ValueError("BROKEN_DIRECT_CHECKPOINT_CHAIN")
    if stage == 2:
        if tuple(prev.get("operator_signature", ())) != tuple(first_operation):
            raise ValueError("MISMATCH_FIRST_OPERATOR_SIGNATURE")
    else:
        for key in ("order", "aux"):
            if prev.get(key) != cur.get(key):
                raise ValueError("MISMATCH_" + key)
        if prev.get("post_second_scalar_filter") is not True:
            raise ValueError("STAGE2_FILTER_NOT_APPLIED")
        if tuple(prev.get("operator_signature", ())) != tuple(second_operation):
            raise ValueError("MISMATCH_SECOND_OPERATOR_SIGNATURE")
    return True
