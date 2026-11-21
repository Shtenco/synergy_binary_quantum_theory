#!/usr/bin/env python3
"""Reference/negative control for gravitational background normalization.

Shows that one global geometry-independent normalization preserves the distinct
lambda^2 and lambda^4 pieces of -log Z(lambda), whereas separately normalizing
Z(lambda) to one at every background erases the entire effective action.

This is bookkeeping/algebra only, not a BQG vacuum-energy calculation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def design(lam: np.ndarray) -> np.ndarray:
    return np.column_stack([lam**2, lam**4, np.ones_like(lam)])


def fit_coefficients(lam: np.ndarray, W: np.ndarray):
    c, *_ = np.linalg.lstsq(design(lam), W, rcond=None)
    return c


def run() -> dict:
    lam = np.array([0.7, 0.85, 1.0, 1.2, 1.45, 1.75], float)
    c2 = 0.37
    c4 = -0.021
    c0 = 1.13  # geometry-independent normalization constant in W=-log Z

    W = c2 * lam**2 + c4 * lam**4 + c0
    Z = np.exp(-W)

    # Global normalization: multiply every Z by one common constant. This only
    # shifts W by a constant and must leave c2,c4 unchanged.
    common_factor = 7.3
    Z_global = common_factor * Z
    W_global = -np.log(Z_global)
    fit_global = fit_coefficients(lam, W_global)

    # Ratio to one frozen reference geometry removes the geometry-independent
    # constant while retaining both scaling terms.
    iref = 2
    W_ratio = -np.log(Z / Z[iref])
    Xratio = np.column_stack([lam**2 - lam[iref]**2, lam**4 - lam[iref]**4])
    fit_ratio, *_ = np.linalg.lstsq(Xratio, W_ratio, rcond=None)

    # Illegal per-background normalization: choosing N(lambda)=1/Z(lambda)
    # forces every normalized amplitude to unity and destroys all background
    # information.
    Z_per_background = Z * (1.0 / Z)
    W_per_background = -np.log(Z_per_background)
    fit_bad = fit_coefficients(lam, W_per_background)

    global_err = max(abs(fit_global[0] - c2), abs(fit_global[1] - c4))
    ratio_err = max(abs(fit_ratio[0] - c2), abs(fit_ratio[1] - c4))
    bad_signal = max(abs(fit_bad[0]), abs(fit_bad[1]))

    # Extensive control: a per-cell vacuum-like weight z0^N gives a free-energy
    # density independent of N. This checks that the extensive observable is
    # -log Z_N / N rather than a separately normalized Z_N=1.
    z0 = np.exp(-0.083)
    Ns = np.array([8, 16, 32, 64, 128], int)
    ZN = z0 ** Ns
    fN = -np.log(ZN) / Ns
    extensive_err = float(np.max(np.abs(fN - 0.083)))

    checks = {
        "global_normalization_preserves_c2_c4": bool(global_err < 2e-12),
        "reference_ratio_preserves_c2_c4": bool(ratio_err < 2e-12),
        "per_background_normalization_erases_signal": bool(np.max(np.abs(W_per_background)) < 1e-14 and bad_signal < 1e-14),
        "extensive_free_energy_density_stable": bool(extensive_err < 2e-14),
    }

    return {
        "status": "background absolute-normalization reference and negative control",
        "passed": bool(all(checks.values())),
        "science_status": "BACKGROUND_NORMALIZATION_FIREWALL",
        "synthetic_true_coefficients": {"lambda2": c2, "lambda4": c4, "constant": c0},
        "global_normalization_fit": {
            "lambda2": float(fit_global[0]),
            "lambda4": float(fit_global[1]),
            "constant": float(fit_global[2]),
            "coefficient_error": float(global_err),
        },
        "reference_geometry_ratio_fit": {
            "reference_lambda": float(lam[iref]),
            "lambda2": float(fit_ratio[0]),
            "lambda4": float(fit_ratio[1]),
            "coefficient_error": float(ratio_err),
        },
        "illegal_per_background_normalization": {
            "max_abs_W": float(np.max(np.abs(W_per_background))),
            "fitted_lambda2": float(fit_bad[0]),
            "fitted_lambda4": float(fit_bad[1]),
        },
        "extensive_control": {
            "N": [int(x) for x in Ns],
            "free_energy_density": [float(x) for x in fN],
            "max_error": extensive_err,
        },
        "checks": checks,
        "interpretation": (
            "A common geometry-independent normalization changes only an additive constant in W. A ratio to one frozen reference background preserves geometry-dependent lambda^2/lambda^4 information. Separately normalizing every geometry to Z=1 erases exactly the information that could contain a gravitational vacuum term."
        ),
        "claim_boundary": "Reference algebra only; c2 and c4 are synthetic and are not BQG predictions.",
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = run()
    text = json.dumps(out, indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
