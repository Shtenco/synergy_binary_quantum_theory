#!/usr/bin/env python3
"""Bounded branch-aware pairwise W5 KKV support pilot."""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import bqg_w5_disjunctive_pairwise_support as DP
import bqg_w5_su2_structural_basis as SB
import bqg_w5_pairwise_factor_graph as PG


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--K", type=int, default=0)
    parser.add_argument("--Jmax2", type=int, default=10)
    parser.add_argument("--source", type=int, default=0)
    parser.add_argument("--triple", type=int, nargs=3, default=(1, 2, 3))
    parser.add_argument("--limit", type=int, default=25000)
    parser.add_argument("--mode", choices=("forward", "adjoint", "both"),
                        default="both")
    args = parser.parse_args()

    key = SB.w5_seed_as_covariant(args.K, args.source)
    initial = PG._initial(key)
    modes = ("forward", "adjoint") if args.mode == "both" else (args.mode,)
    report = {
        "seed": "W5",
        "K": args.K,
        "source": args.source,
        "triple": list(args.triple),
        "Jmax2": args.Jmax2,
        "entry_limit": args.limit,
        "representation": "disjunctive oriented-spec alternatives; each branch is a pairwise CSP",
        "branch_identity_preserved": True,
        "within_branch_pair_correlations_preserved": True,
        "higher_order_correlations_preserved": False,
        "joint_keys_enumerated": False,
        "amplitudes_evaluated": False,
        "results": [],
    }
    for mode in modes:
        item = {"mode": mode}
        try:
            branches = DP.factor_graph_triple_branches(
                [initial], args.source, *args.triple, args.Jmax2,
                adjoint=(mode == "adjoint"), limit=args.limit
            )
            item["status"] = "completed"
            item.update(DP.profile(branches))
        except RuntimeError as exc:
            item["status"] = "guard_hit"
            item["guard_message"] = str(exc)
        report["results"].append(item)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
