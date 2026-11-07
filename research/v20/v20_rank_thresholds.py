#!/usr/bin/env python3
"""BQG V20: reproducible finite W5 kernel-nullity NECESSARY bounds.

This code replays the genuine archived complex K5 W5 Gram numerically.
It does NOT compute a new SU(2) amplitude, physical Lorentzian H, or HDA.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

import numpy as np

EXPECTED_SOURCE_SHA256 = "950423420850c887f8facb747be66a44cc4ed007655c6b04aac8ce21b347dbc9"
EXPECTED_GIT_BLOB_SHA = "0231c54d2e9a7fc9401995c999d3989f3f1f7648"


def source_read(path):
    raw = Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != EXPECTED_SOURCE_SHA256:
        raise ValueError("SOURCE_SHA256_MISMATCH")
    src = json.loads(raw, parse_float=Decimal)
    if src.get("source_git_blob_sha") != EXPECTED_GIT_BLOB_SHA:
        raise ValueError("SOURCE_GIT_BLOB_MISMATCH")
    if not src.get("known_HDA_not_proved"):
        raise ValueError("FALSE_PHYSICAL_CLAIM")
    return src


def exact_printed_gershgorin(src):
    """Exact rational lower eigenvalue bound on printed decimal entries."""
    data = src["complex_gram_re_im"]
    if len(data) != 5 or any(len(row) != 5 for row in data):
        raise ValueError("WRONG_W5_SHAPE")
    z = [[(Fraction(re), Fraction(im)) for re, im in row] for row in data]
    bounds = []
    for i in range(5):
        if z[i][i][1] != 0:
            raise ValueError("DIAGONAL_NOT_REAL")
        for j in range(5):
            if z[i][j][0] != z[j][i][0] or z[i][j][1] != -z[j][i][1]:
                raise ValueError("NOT_HERMITIAN")
        b = z[i][i][0] - sum(
            abs(z[i][j][0]) + abs(z[i][j][1])
            for j in range(5) if i != j
        )
        bounds.append(b)
    if min(bounds) <= 0:
        raise ValueError("POSITIVITY_NOT_CERTIFIED")
    return min(bounds)


def replay(path):
    src = source_read(path)
    z = src["complex_gram_re_im"]
    G = np.array([[complex(float(re), float(im)) for re, im in row]
                  for row in z], dtype=np.complex128)
    if not np.all(np.isfinite(G)):
        raise ValueError("NONFINITE_GRAM")
    if np.linalg.norm(G - G.conj().T) != 0:
        raise ValueError("NOT_FLOAT_HERMITIAN")
    eig = np.linalg.eigvalsh(G)
    np.testing.assert_allclose(eig, [float(x) for x in src["published_eigvals"]],
                               atol=1e-11, rtol=1e-12)
    exact_bound = exact_printed_gershgorin(src)
    with localcontext() as context:
        context.prec = 65
        bound_decimal = str(Decimal(exact_bound.numerator) /
                            Decimal(exact_bound.denominator))
    a = 2.0 / 3.0
    thresholds = [
        {
            "k": k,
            "necessary_Q_operator_norm": a * math.sqrt(eig[k-1]),
            "necessary_Q_Frobenius_norm": a * math.sqrt(float(np.sum(eig[:k]))),
            "not_sufficient": True
        } for k in range(1,6)
    ]
    cmin = float(src["commutator_smin"])
    if src["commutator_rank"] != 10 or cmin <= 0:
        raise ValueError("BAD_PL16_COMMUTATOR_EVIDENCE")
    residuals = [
        {"maximum_rank_K": r,
         "necessary_HDA_residual_Frobenius_lower": math.sqrt(10-r) * (4/9) * cmin,
         "necessary_HDA_residual_operator_lower": (4/9) * cmin}
        for r in range(10)
    ]
    return {
        "version": "V20_REPLAY",
        "source_git_blob_sha": EXPECTED_GIT_BLOB_SHA,
        "source_snapshot_sha256": EXPECTED_SOURCE_SHA256,
        "K5_W5": {
            "eigenvalues": eig.tolist(),
            "Gershgorin_exact_printed_fraction": f"{exact_bound.numerator}/{exact_bound.denominator}",
            "Gershgorin_decimal": bound_decimal,
            "kernel_nullity_compensator_thresholds": thresholds
        },
        "PL16_ten_EE_probes": {
            "minimum_singular_value": cmin,
            "rank_limited_compensator_residual_floors": residuals
        },
        "never_stack_K5_with_PL16": True,
        "new_physical_Lorentzian_amplitudes": False,
        "full_HDA_proved": False
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, default=Path(__file__).with_name("SOURCE_SNAPSHOT.json"))
    p.add_argument("--out", type=Path, default=Path("RESULTS_V20.json"))
    opts = p.parse_args()
    result = replay(opts.source)
    opts.out.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
