#!/usr/bin/env python3
"""Bounded factorized operator-level W5 KKV support pilot.

Reports conservative marginal support envelopes only. It does not evaluate
coefficients, preserve cross-node correlations, or certify physical KKV output.
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import bqg_w5_su2_structural_basis as SB
import bqg_w5_factorized_operator_support as FO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--K", type=int, default=0)
    parser.add_argument("--Jmax2", type=int, default=10)
    parser.add_argument("--max-factor-labels", type=int, default=100000)
    parser.add_argument("--source", type=int, default=0)
    parser.add_argument("--triple", type=int, nargs=3, default=(1, 2, 3))
    parser.add_argument("--mode", choices=("forward", "adjoint", "both", "node"),
                        default="both")
    args = parser.parse_args()

    key = SB.w5_seed_as_covariant(args.K, args.source)
    lifted = SB.lift_key_to_all_charged(key)
    initial = {
        lifted[0]: tuple(frozenset((lab,)) for lab in lifted[1])
    }
    modes = ("forward", "adjoint") if args.mode == "both" else (args.mode,)
    report = {
        "seed": "W5",
        "K": args.K,
        "source": args.source,
        "triple": list(args.triple),
        "Jmax2": args.Jmax2,
        "max_factor_labels": args.max_factor_labels,
        "amplitudes_evaluated": False,
        "joint_keys_enumerated": False,
        "cross_node_correlations_preserved": False,
        "interpretation": "conservative Cartesian-product upper bound from local marginals",
        "results": [],
    }

    for mode in modes:
        item = {"mode": mode}
        try:
            if mode == "node":
                support = FO.factorized_node_kkv(
                    initial, args.source, args.Jmax2,
                    max_factor_labels=args.max_factor_labels,
                )
            else:
                support = FO.factorized_triple(
                    initial, args.source, *args.triple, args.Jmax2,
                    adjoint=(mode == "adjoint"),
                    max_factor_labels=args.max_factor_labels,
                )
            item["status"] = "completed"
            item.update(FO.factorized_support_profile(support))
            item["max_doubled_link_spin"] = max(
                (max(spins) for spins in support), default=0
            )
        except RuntimeError as exc:
            item["status"] = "guard_hit"
            item["guard_message"] = str(exc)
        report["results"].append(item)

    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
