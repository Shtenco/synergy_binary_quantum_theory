#!/usr/bin/env python3
"""Exact symbolic discriminator for scalar IR response classes.

This is an interpretation/reference gate for a future derived physical scalar
kernel.  It proves the distinct small-k behavior of:

1. analytic noncritical K0 + K2 k^2;
2. Poisson/constraint Zk k^2;
3. a frequency-dependent propagating scalar kernel.

No BQG coefficient is supplied by this gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sympy as sp


def run():
    k, w = sp.symbols("k w", real=True)
    K0, K2, Zk, Zw, m2 = sp.symbols("K0 K2 Zk Zw m2", nonzero=True)

    analytic = K0 + K2 * k**2
    analytic_inv = sp.series(1 / analytic, k, 0, 5).removeO()
    analytic_limit = sp.simplify(sp.limit(1 / analytic, k, 0))
    analytic_k2_inv_limit = sp.simplify(sp.limit(k**2 / analytic, k, 0))

    poisson = Zk * k**2
    poisson_inv = sp.simplify(1 / poisson)
    poisson_residue = sp.simplify(sp.limit(k**2 / poisson, k, 0))

    dynamic = Zw * w**2 - Zk * k**2 - m2
    roots_w2 = sp.solve(sp.Eq(dynamic, 0), w**2)
    expected_root = sp.simplify((m2 + Zk * k**2) / Zw)
    dK_dw2 = sp.simplify(sp.diff(dynamic, w) / (2 * w))

    checks = {
        "analytic_inverse_finite_at_k0": analytic_limit == 1 / K0,
        "analytic_channel_has_no_kminus2_residue": analytic_k2_inv_limit == 0,
        "poisson_inverse_exactly_kminus2": poisson_inv == 1 / (Zk * k**2),
        "poisson_static_residue_is_inverse_Zk": poisson_residue == 1 / Zk,
        "dynamic_kernel_has_frequency_pole": len(roots_w2) == 1 and sp.simplify(roots_w2[0] - expected_root) == 0,
        "dynamic_w2_slope_is_Zw": sp.simplify(dK_dw2 - Zw) == 0,
    }

    return {
        "status": "exact scalar IR response-class discriminator",
        "passed": bool(all(checks.values())),
        "analytic_noncritical": {
            "kernel": "K0 + K2 k^2",
            "inverse_series_through_k4": str(analytic_inv),
            "inverse_k0_limit": "1/K0",
            "k2_times_inverse_k0_limit": "0",
            "interpretation": "finite/contact-like IR response; no k^-2 long-range law",
        },
        "poisson_constraint": {
            "kernel": "Zk k^2",
            "inverse": "1/(Zk k^2)",
            "k2_residue": "1/Zk",
            "interpretation": "long-range static response; not by itself a new propagating degree of freedom",
        },
        "propagating_scalar": {
            "kernel": "Zw w^2 - Zk k^2 - m2",
            "pole_equation": "w^2=(m2+Zk k^2)/Zw",
            "w2_derivative": "Zw",
            "interpretation": "candidate extra physical scalar only after projector/gauge reduction and residue/stability tests",
        },
        "checks": checks,
        "claim_boundary": (
            "This symbolic gate classifies possible infrared forms. It does not derive "
            "K0, K2, Zk, Zw, m2, mu, Sigma, dark matter or dark energy from BQG."
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
