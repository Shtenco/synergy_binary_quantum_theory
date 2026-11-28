#!/usr/bin/env python3
"""Bounded structural-support pilot; no amplitudes or tolerance pruning."""
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
    ap.add_argument("--Jmax2", type=int, default=2,
                    help="maximum doubled link spin")
    ap.add_argument("--max-states", type=int, default=5000)
    ap.add_argument("--source", type=int, default=0)
    ap.add_argument("--triple", nargs=3, type=int, default=(1, 2, 3))
    ap.add_argument("--w5-K", type=int, choices=SB.W5_K, default=None,
                    help="use an actual all-link j=2 W5 seed instead of the minimal diagnostic seed")
    ap.add_argument("--include-adjoint", action="store_true", default=True)
    args = ap.parse_args()

    # Minimal nontrivial seed: all links spin zero, all five nodes Gauss.
    # This is a diagnostic pilot, not the physical W5 j=2 seed.
    if args.w5_K is None:
        key = ((0,) * 10, tuple(("G", 0) for _ in range(5)))
        seed_name = "K5 all-link j=0 Gauss seed (diagnostic only)"
    else:
        key = SB.w5_seed_as_covariant(args.w5_K, args.source)
        seed_name = f"W5 all-link j=2 seed K={args.w5_K}"
    overflow = None
    try:
        support = SB.t_sequence_structural_support(
            key, args.source, *args.triple, args.Jmax2,
            charged_nodes=(args.source, args.triple[0]),
            include_adjoint=args.include_adjoint,
            max_states=args.max_states,
            volume_v=args.source,
        )
    except RuntimeError as exc:
        support = set()
        overflow = str(exc)
    result = {
        "status": "bounded structural support pilot",
        "seed": seed_name,
        "Jmax2": args.Jmax2,
        "source": args.source,
        "triple": list(args.triple),
        "include_adjoint": args.include_adjoint,
        "max_states": args.max_states,
        "support_size": len(support) if overflow is None else None,
        "state_limit_reached": overflow is not None,
        "guard_message": overflow,
        "spin_sector_count": len({key[0] for key in support}),
        "max_link_spin2": max((max(key[0]) for key in support), default=0),
        "amplitudes_evaluated": False,
        "threshold_pruning": False,
        "physical_W5_claim": False,
    }
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
