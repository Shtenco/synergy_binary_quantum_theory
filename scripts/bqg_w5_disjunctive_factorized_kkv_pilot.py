#!/usr/bin/env python3
"""Bounded branch-aware disjunctive W5 KKV support pilot."""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import bqg_w5_disjunctive_factorized_support as DF
import bqg_w5_su2_structural_basis as SB


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--K", type=int, default=0)
    parser.add_argument("--Jmax2", type=int, default=10)
    parser.add_argument("--source", type=int, default=0)
    parser.add_argument("--triple", type=int, nargs=3, default=(1, 2, 3))
    parser.add_argument("--limit", type=int, default=10000)
    parser.add_argument("--mode", choices=("forward", "adjoint", "both"),
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
        "entry_limit": args.limit,
        "representation": "disjunctive alternatives; each branch stores five node marginals",
        "branch_identity_preserved": True,
        "within_branch_cross_node_correlations_preserved": False,
        "joint_keys_enumerated": False,
        "amplitudes_evaluated": False,
        "results": [],
    }
    for mode in modes:
        item = {"mode": mode}
        try:
            branches = DF.factorized_triple_branches(
                [initial], args.source, *args.triple, args.Jmax2,
                adjoint=(mode == "adjoint"), limit=args.limit
            )
            item["status"] = "completed"
            item.update(DF.profile(branches))
        except RuntimeError as exc:
            item["status"] = "guard_hit"
            item["guard_message"] = str(exc)
        report["results"].append(item)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
