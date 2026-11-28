#!/usr/bin/env python3
"""Bounded end-to-end C(K)/KKV structural pilot on one W5 seed.

This enumerates conservative joint label support only. It does not evaluate
matrix coefficients, cancellations, or physical W5 amplitudes. A guard hit is
reported as an expansion diagnostic and never truncates the support.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import bqg_w5_su2_structural_basis as SB


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--K", type=int, choices=SB.W5_K, default=0)
    ap.add_argument("--source", type=int, default=0)
    ap.add_argument("--target", type=int, default=1)
    ap.add_argument("--triple", nargs=3, type=int, default=(1, 2, 3))
    ap.add_argument("--Jmax2", type=int, default=10)
    ap.add_argument("--max-states", type=int, default=5000)
    ap.add_argument("--mode", choices=("C_K", "forward", "adjoint", "node"),
                    default="C_K")
    ap.add_argument("--i", type=int, choices=(0, 1), default=0)
    ap.add_argument("--j", type=int, choices=(0, 1), default=0)
    args = ap.parse_args()

    key = SB.w5_seed_as_covariant(args.K, args.source)
    support = None
    overflow = None
    try:
        if args.mode == "C_K":
            support = SB.C_K_structural_support(
                key, args.source, args.target, args.Jmax2,
                args.i, args.j, max_states=args.max_states
            )
        elif args.mode == "forward":
            support = SB.triple_forward_structural_support(
                key, args.source, *args.triple, args.Jmax2,
                max_states=args.max_states
            )
        elif args.mode == "adjoint":
            support = SB.triple_adjoint_structural_support(
                key, args.source, *args.triple, args.Jmax2,
                max_states=args.max_states
            )
        else:
            support = SB.node_kkv_structural_support(
                key, args.source, args.Jmax2,
                max_states=args.max_states
            )
    except RuntimeError as exc:
        overflow = str(exc)

    result = {
        "status": "bounded joint structural support pilot",
        "mode": args.mode,
        "seed": f"W5 all-link j=2 seed K={args.K}",
        "source": args.source,
        "target": args.target if args.mode == "C_K" else None,
        "triple": list(args.triple) if args.mode in ("forward", "adjoint") else None,
        "matrix_component": [args.i, args.j] if args.mode == "C_K" else None,
        "Jmax2": args.Jmax2,
        "max_states": args.max_states,
        "support_size": len(support) if overflow is None else None,
        "state_limit_reached": overflow is not None,
        "guard_message": overflow,
        "spin_sector_count": len({item[0] for item in support}) if support is not None else None,
        "max_link_spin2": max((max(item[0]) for item in support), default=None)
            if support is not None else None,
        "amplitudes_evaluated": False,
        "threshold_pruning": False,
        "physical_W5_claim": False,
    }
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
