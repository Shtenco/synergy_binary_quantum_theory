#!/usr/bin/env python3
"""Exact homogeneity discriminator: Einstein/Regge curvature versus 4-volume.

This is a structural background diagnostic, not a physical FLRW calculation.
Under a uniform four-dimensional length rescaling l -> lambda*l,

  hinge areas A_h -> lambda^2 A_h,
  Regge deficit angles stay dimensionless,
  4-simplex volumes -> lambda^4 V_4.

Equivalently g -> lambda^2 g gives sqrt(g) R -> lambda^2 sqrt(g) R
and sqrt(g) -> lambda^4 sqrt(g).  The gate freezes these distinct homogeneity
exponents so a future physical Gamma[g] cannot mistake the existing bare
Einstein/Regge term for an independently generated cosmological-volume term.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sympy as sp


def log_slope(expr, lam):
    return sp.simplify(lam * sp.diff(expr, lam) / expr)


def run():
    lam = sp.symbols("lambda", positive=True)
    A0, delta0, V0 = sp.symbols("A0 delta0 V0", nonzero=True)

    regge = lam**2 * A0 * delta0
    volume = lam**4 * V0
    regge_slope = sp.simplify(log_slope(regge, lam))
    volume_slope = sp.simplify(log_slope(volume, lam))

    # Continuum Weyl-scaling bookkeeping in D=4:
    # sqrt(g) -> lambda^4 sqrt(g), R -> lambda^-2 R.
    sqrtg_scale = lam**4
    R_scale = lam**-2
    EH_scale = sp.simplify(sqrtg_scale * R_scale)

    checks = {
        "hinge_Regge_term_has_uniform_length_exponent_2": regge_slope == 2,
        "four_volume_term_has_uniform_length_exponent_4": volume_slope == 4,
        "continuum_sqrtgR_has_uniform_metric_length_exponent_2": EH_scale == lam**2,
        "curvature_and_volume_exponents_are_distinct": regge_slope != volume_slope,
    }

    return {
        "status": "exact Regge/EH versus cosmological-volume homogeneity discriminator",
        "passed": bool(all(checks.values())),
        "uniform_length_rescaling": "l -> lambda l (equivalently g -> lambda^2 g)",
        "Regge_EH_scaling": "sum_h A_h delta_h -> lambda^2 sum_h A_h delta_h",
        "four_volume_scaling": "sum_sigma V4_sigma -> lambda^4 sum_sigma V4_sigma",
        "Regge_EH_logarithmic_slope": int(regge_slope),
        "four_volume_logarithmic_slope": int(volume_slope),
        "checks": checks,
        "interpretation": (
            "The already-tested bare Regge/EH bridge supplies the curvature homogeneity class, "
            "not an independently generated cosmological-volume term. A future physical-history "
            "Gamma may generate such a relevant term, but it must be read from that Gamma rather "
            "than inferred from the existing curvature reconstruction."
        ),
        "claim_boundary": (
            "Uniform four-dimensional Weyl/length scaling is not FLRW time evolution and this "
            "gate does not compute rho_hist(a), Lambda, H(a), vacuum energy, or dark energy. It "
            "only supplies an exact action-term discriminator for the future Gamma_FLRW analysis."
        ),
    }


def main():
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
