#!/usr/bin/env python3
"""Fail-closed quantum HDA habitat residual certificate.

This gate is deliberately downstream of target construction.  It certifies only
an ACTUAL quantum residual packet on one declared graph-changing habitat:

    Delta_kl(lambda) = [H_k,H_l] - D^target_kl
                     = r0 + lambda r1 + lambda^2 r2,

where
    r0 = HH_EE - D0,
    r1 = HH_EL + HH_LE - D1,
    r2 = HH_LL - D2.

The scalar lambda = 1 + beta^2 is never fitted by this gate.  If the three
component residuals vanish separately, closure is beta-independent.  Otherwise
lambda must be supplied as a preregistered microscopic/convention value and is
only evaluated, never optimized.

Sparse vectors are JSON mappings basis_key -> [real, imag] (real scalars are
also accepted).  Aggregate norms concatenate declared (pair,column) residuals;
there is no assumption that different constraint pairs are physical Hilbert
states to be superposed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple

SCHEMA_PACKET = "BQG_QUANTUM_HDA_HABITAT_PACKET_V1"
SCHEMA_CERT = "BQG_QUANTUM_HDA_RESIDUAL_CERTIFICATE_V1"

Sparse = Dict[str, complex]


def _complex(value) -> complex:
    if isinstance(value, (int, float)):
        return complex(float(value), 0.0)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) == 2:
        return complex(float(value[0]), float(value[1]))
    if isinstance(value, Mapping):
        if "re" in value or "im" in value:
            return complex(float(value.get("re", 0.0)), float(value.get("im", 0.0)))
    raise ValueError(f"invalid complex coefficient: {value!r}")


def _sparse(obj) -> Sparse:
    if obj is None:
        return {}
    if not isinstance(obj, Mapping):
        raise ValueError("sparse state must be a JSON object")
    out: Sparse = {}
    for key, value in obj.items():
        z = _complex(value)
        if z != 0.0j:
            out[str(key)] = z
    return out


def _lincomb(*terms: Tuple[complex, Mapping[str, complex]]) -> Sparse:
    out: Sparse = {}
    for factor, state in terms:
        for key, value in state.items():
            z = out.get(key, 0.0j) + factor * value
            if z == 0.0j:
                out.pop(key, None)
            else:
                out[key] = z
    return out


def _inner(a: Mapping[str, complex], b: Mapping[str, complex]) -> complex:
    # Iterate the smaller map. Missing basis coefficients are exactly zero.
    if len(a) > len(b):
        return _inner(b, a).conjugate()
    return sum((value.conjugate() * b.get(key, 0.0j) for key, value in a.items()), 0.0j)


def _norm2(a: Mapping[str, complex]) -> float:
    return float(sum((value.real * value.real + value.imag * value.imag for value in a.values()), 0.0))


def _norm(a: Mapping[str, complex]) -> float:
    return math.sqrt(max(0.0, _norm2(a)))


def _canonical_sha256(obj) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _expected_cross_product(packet) -> Tuple[set, set, set]:
    coverage = packet.get("coverage", {})
    pairs = [str(x) for x in coverage.get("expected_pairs", [])]
    columns = [str(x) for x in coverage.get("expected_columns", [])]
    expected = {(p, c) for p in pairs for c in columns}
    return set(pairs), set(columns), expected


def _fit_refinement(series):
    """Fit y=A*x^p+B using a preregistered deterministic p-grid.

    This is a diagnostic/refinement certificate, never a beta fit.  x must tend
    toward zero under refinement.  The linear A,B fit is solved exactly for
    each p; the p with minimum SSE is reported.
    """
    pts = [(float(row["regulator"]), float(row["residual_norm"])) for row in series]
    if len(pts) < 3:
        raise ValueError("refinement_series needs at least three points")
    if any(x <= 0.0 or y < 0.0 for x, y in pts):
        raise ValueError("refinement regulator must be >0 and residual_norm >=0")
    pts.sort(reverse=True)  # coarse -> fine for readable diagnostics
    best = None
    # Frozen grid: p in [0.10, 6.00] by 0.01.
    for ip in range(10, 601):
        p = ip / 100.0
        u = [x ** p for x, _ in pts]
        y = [yy for _, yy in pts]
        n = float(len(pts))
        su = sum(u); sy = sum(y)
        suu = sum(v * v for v in u); suy = sum(v * yy for v, yy in zip(u, y))
        den = n * suu - su * su
        if abs(den) < 1e-30:
            continue
        A = (n * suy - su * sy) / den
        B = (sy - A * su) / n
        sse = sum((yy - (A * v + B)) ** 2 for v, yy in zip(u, y))
        candidate = (sse, p, A, B)
        if best is None or candidate < best:
            best = candidate
    if best is None:
        raise ValueError("refinement fit is singular")
    sse, p, A, B = best
    return {
        "model": "residual=A*regulator^p+B",
        "p_grid": [0.10, 6.00, 0.01],
        "points": [{"regulator": x, "residual_norm": y} for x, y in pts],
        "A": A,
        "p": p,
        "B": B,
        "sse": sse,
    }


def certify(packet: Mapping) -> Mapping:
    if packet.get("schema") != SCHEMA_PACKET:
        raise ValueError(f"expected schema {SCHEMA_PACKET}")

    required_hashes = ("habitat_hash", "domain_hash", "constraint_packet_hash", "convention_hash")
    hashes = {name: str(packet.get(name, "")) for name in required_hashes}
    hashes_present = all(bool(v) for v in hashes.values())

    tol = packet.get("tolerances", {})
    residual_tol = float(tol.get("residual_norm", 1e-12))
    leakage_tol = float(tol.get("leakage", 1e-12))
    recoupling_tol = float(tol.get("recoupling", 1e-12))
    refinement_B_tol = float(tol.get("refinement_intercept", residual_tol))
    refinement_p_min = float(tol.get("refinement_min_exponent", 0.0))

    errors = packet.get("numerical_errors", {})
    cutoff_leakage = float(errors.get("cutoff_leakage", math.inf))
    habitat_leakage = float(errors.get("habitat_leakage", math.inf))
    recoupling_error = float(errors.get("recoupling_error", math.inf))
    numerical_controlled = (
        cutoff_leakage <= leakage_tol
        and habitat_leakage <= leakage_tol
        and recoupling_error <= recoupling_tol
    )

    pairs, columns, expected = _expected_cross_product(packet)
    records = packet.get("records", [])
    if not isinstance(records, list):
        raise ValueError("records must be a list")

    seen = set()
    duplicate_records = []
    unexpected_records = []
    record_results = []
    sum00 = 0.0
    sum11 = 0.0
    sum22 = 0.0
    sum01 = 0.0j
    sum02 = 0.0j
    sum12 = 0.0j
    support = {"r0": 0, "r1": 0, "r2": 0}
    record_hash_consistent = True

    for row in records:
        pair = str(row.get("pair", ""))
        column = str(row.get("column", ""))
        key = (pair, column)
        if key in seen:
            duplicate_records.append([pair, column])
        seen.add(key)
        if expected and key not in expected:
            unexpected_records.append([pair, column])

        for name, value in hashes.items():
            if name in row and str(row[name]) != value:
                record_hash_consistent = False

        ee = _sparse(row.get("HH_EE"))
        el = _sparse(row.get("HH_EL"))
        le = _sparse(row.get("HH_LE"))
        ll = _sparse(row.get("HH_LL"))
        d0 = _sparse(row.get("D0"))
        d1 = _sparse(row.get("D1"))
        d2 = _sparse(row.get("D2"))

        r0 = _lincomb((1.0, ee), (-1.0, d0))
        r1 = _lincomb((1.0, el), (1.0, le), (-1.0, d1))
        r2 = _lincomb((1.0, ll), (-1.0, d2))

        n00 = _norm2(r0); n11 = _norm2(r1); n22 = _norm2(r2)
        i01 = _inner(r0, r1); i02 = _inner(r0, r2); i12 = _inner(r1, r2)
        sum00 += n00; sum11 += n11; sum22 += n22
        sum01 += i01; sum02 += i02; sum12 += i12
        support["r0"] += len(r0); support["r1"] += len(r1); support["r2"] += len(r2)
        record_results.append({
            "pair": pair,
            "column": column,
            "norm_r0": math.sqrt(max(0.0, n00)),
            "norm_r1": math.sqrt(max(0.0, n11)),
            "norm_r2": math.sqrt(max(0.0, n22)),
            "support_r0": len(r0),
            "support_r1": len(r1),
            "support_r2": len(r2),
        })

    missing_records = sorted([list(x) for x in (expected - seen)]) if expected else []
    domain_complete = bool(expected) and not missing_records and not duplicate_records and not unexpected_records

    c0 = sum00
    c1 = 2.0 * sum01.real
    c2 = sum11 + 2.0 * sum02.real
    c3 = 2.0 * sum12.real
    c4 = sum22
    coeffs = [c0, c1, c2, c3, c4]

    component_norms = {
        "r0_EE_minus_D0": math.sqrt(max(0.0, sum00)),
        "r1_EL_plus_LE_minus_D1": math.sqrt(max(0.0, sum11)),
        "r2_LL_minus_D2": math.sqrt(max(0.0, sum22)),
    }
    beta_independent_closure = all(value <= residual_tol for value in component_norms.values())

    lambda_value = packet.get("lambda_preregistered", None)
    lambda_preregistered = lambda_value is not None
    evaluated_norm = None
    lambda_residual_pass = False
    if lambda_preregistered:
        lam = float(lambda_value)
        if lam < 0.0 or not math.isfinite(lam):
            raise ValueError("lambda_preregistered must be finite and >=0")
        n2 = c0 + c1 * lam + c2 * lam**2 + c3 * lam**3 + c4 * lam**4
        # Round-off may produce a tiny negative value for an exact cancellation.
        if n2 < 0.0 and abs(n2) <= 100.0 * sys.float_info.epsilon * max(1.0, sum(abs(x) for x in coeffs)):
            n2 = 0.0
        if n2 < 0.0:
            raise ValueError("negative residual norm^2: packet/inner-product inconsistency")
        evaluated_norm = math.sqrt(n2)
        lambda_residual_pass = evaluated_norm <= residual_tol

    refinement_required = bool(packet.get("require_refinement_certificate", False))
    refinement_fit = None
    refinement_pass = not refinement_required
    if packet.get("refinement_series"):
        refinement_fit = _fit_refinement(packet["refinement_series"])
        refinement_pass = (
            abs(float(refinement_fit["B"])) <= refinement_B_tol
            and float(refinement_fit["p"]) >= refinement_p_min
        )
    elif refinement_required:
        refinement_pass = False

    closure_at_declared_lambda = lambda_preregistered and lambda_residual_pass
    residual_physics_pass = beta_independent_closure or closure_at_declared_lambda

    preconditions = {
        "hashes_present": hashes_present,
        "record_hash_consistent": record_hash_consistent,
        "domain_complete": domain_complete,
        "numerical_controlled": numerical_controlled,
        "residual_physics_pass": residual_physics_pass,
        "refinement_pass": refinement_pass,
        "lambda_was_not_fitted": True,
    }
    certified = all(preconditions.values())

    if beta_independent_closure:
        science_status = "BETA_INDEPENDENT_QUANTUM_HDA_CLOSURE" if certified else "STRUCTURAL_CLOSURE_BUT_CERTIFICATE_INCOMPLETE"
    elif closure_at_declared_lambda:
        science_status = "PREREGISTERED_LAMBDA_QUANTUM_HDA_CLOSURE" if certified else "DECLARED_LAMBDA_RESIDUAL_SMALL_BUT_CERTIFICATE_INCOMPLETE"
    else:
        science_status = "OPEN_QUANTUM_HDA_RESIDUAL"

    return {
        "schema": SCHEMA_CERT,
        "source_packet_schema": SCHEMA_PACKET,
        "source_packet_sha256": _canonical_sha256(packet),
        "hashes": hashes,
        "jmax": packet.get("jmax"),
        "coverage": {
            "expected_pair_count": len(pairs),
            "expected_column_count": len(columns),
            "expected_record_count": len(expected),
            "actual_record_count": len(records),
            "missing_records": missing_records,
            "duplicate_records": duplicate_records,
            "unexpected_records": unexpected_records,
        },
        "numerical_errors": {
            "cutoff_leakage": cutoff_leakage,
            "habitat_leakage": habitat_leakage,
            "recoupling_error": recoupling_error,
        },
        "tolerances": {
            "residual_norm": residual_tol,
            "leakage": leakage_tol,
            "recoupling": recoupling_tol,
            "refinement_intercept": refinement_B_tol,
            "refinement_min_exponent": refinement_p_min,
        },
        "component_residual_norms": component_norms,
        "component_support_sum": support,
        "norm_squared_polynomial": {
            "variable": "lambda=1+beta^2",
            "coefficients_c0_to_c4": coeffs,
            "meaning": "||r0+lambda*r1+lambda^2*r2||_aggregate^2",
            "lambda_fit_performed": False,
            "lambda_fit_used_for_certification": False,
        },
        "lambda_evaluation": {
            "lambda_preregistered": lambda_value,
            "evaluated_residual_norm": evaluated_norm,
            "pass": lambda_residual_pass if lambda_preregistered else False,
        },
        "beta_independent_closure": beta_independent_closure,
        "refinement": {
            "required": refinement_required,
            "fit": refinement_fit,
            "pass": refinement_pass,
        },
        "record_results": record_results,
        "preconditions": preconditions,
        "quantum_habitat_residual_certified": certified,
        "science_status": science_status,
        "forbidden_inference": "A target-only Dtarget certificate or a fitted lambda cannot certify quantum habitat closure.",
    }


def _synthetic_packet(mode: str):
    base = {
        "schema": SCHEMA_PACKET,
        "habitat_hash": "habitat:test:v1",
        "domain_hash": "domain:test:v1",
        "constraint_packet_hash": "constraints:test:v1",
        "convention_hash": "convention:test:v1",
        "jmax": 3.5,
        "coverage": {"expected_pairs": ["0,1"], "expected_columns": ["psi0"]},
        "numerical_errors": {"cutoff_leakage": 0.0, "habitat_leakage": 0.0, "recoupling_error": 0.0},
        "tolerances": {"residual_norm": 1e-12, "leakage": 1e-12, "recoupling": 1e-12},
        "records": [{
            "pair": "0,1", "column": "psi0",
            "HH_EE": {"a": [1.0, 0.0]}, "D0": {"a": [1.0, 0.0]},
            "HH_EL": {"b": [0.5, 0.0]}, "HH_LE": {"b": [-0.5, 0.0]}, "D1": {},
            "HH_LL": {"c": [2.0, 0.0]}, "D2": {"c": [2.0, 0.0]},
        }],
    }
    if mode == "exact":
        return base
    if mode == "wrong-target":
        base["records"][0]["D0"] = {"a": [0.0, 0.0]}
        return base
    if mode == "missing":
        base["coverage"]["expected_columns"] = ["psi0", "psi1"]
        return base
    if mode == "hash-mismatch":
        base["records"][0]["domain_hash"] = "domain:wrong"
        return base
    if mode == "declared-lambda":
        # r0=1, r1=-1, r2=0 closes only at preregistered lambda=1.
        base["records"][0] = {
            "pair": "0,1", "column": "psi0",
            "HH_EE": {"a": 1.0}, "D0": {},
            "HH_EL": {"a": -1.0}, "HH_LE": {}, "D1": {},
            "HH_LL": {}, "D2": {},
        }
        base["lambda_preregistered"] = 1.0
        return base
    raise ValueError(mode)


def self_test() -> None:
    exact = certify(_synthetic_packet("exact"))
    assert exact["quantum_habitat_residual_certified"] is True
    assert exact["beta_independent_closure"] is True

    wrong = certify(_synthetic_packet("wrong-target"))
    assert wrong["quantum_habitat_residual_certified"] is False
    assert wrong["science_status"] == "OPEN_QUANTUM_HDA_RESIDUAL"

    missing = certify(_synthetic_packet("missing"))
    assert missing["quantum_habitat_residual_certified"] is False
    assert missing["preconditions"]["domain_complete"] is False

    mismatch = certify(_synthetic_packet("hash-mismatch"))
    assert mismatch["quantum_habitat_residual_certified"] is False
    assert mismatch["preconditions"]["record_hash_consistent"] is False

    declared = certify(_synthetic_packet("declared-lambda"))
    assert declared["quantum_habitat_residual_certified"] is True
    assert declared["beta_independent_closure"] is False
    assert declared["lambda_evaluation"]["pass"] is True
    assert declared["norm_squared_polynomial"]["lambda_fit_performed"] is False

    print("PASS quantum_hda_habitat_residual_gate self-test")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--packet", type=Path)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return 0
    if args.packet is None or args.output is None:
        ap.error("--packet and --output are required unless --self-test is used")

    packet = json.loads(args.packet.read_text(encoding="utf-8"))
    certificate = certify(packet)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "science_status": certificate["science_status"],
        "quantum_habitat_residual_certified": certificate["quantum_habitat_residual_certified"],
        "component_residual_norms": certificate["component_residual_norms"],
        "lambda_evaluation": certificate["lambda_evaluation"],
    }, indent=2, sort_keys=True))
    # Fail closed as a command-line gate.
    return 0 if certificate["quantum_habitat_residual_certified"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
