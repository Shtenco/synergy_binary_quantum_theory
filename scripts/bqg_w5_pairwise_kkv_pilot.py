#!/usr/bin/env python3
"""Bounded W5 pairwise factor-graph KKV support pilot."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bqg_w5_su2_structural_basis as SB
import bqg_w5_pairwise_factor_graph as FG


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--K", type=int, default=0)
    p.add_argument("--Jmax2", type=int, default=10)
    p.add_argument("--source", type=int, default=0)
    p.add_argument("--triple", type=int, nargs=3, default=(1, 2, 3))
    p.add_argument("--limit", type=int, default=50000)
    p.add_argument("--mode", choices=("forward", "adjoint", "both"), default="both")
    a = p.parse_args()
    key = SB.w5_seed_as_covariant(a.K, a.source)
    initial = FG._initial(key)
    modes = ("forward", "adjoint") if a.mode == "both" else (a.mode,)
    report = {
        "seed": "W5", "K": a.K, "source": a.source,
        "triple": list(a.triple), "Jmax2": a.Jmax2,
        "entry_limit": a.limit,
        "representation": "five unary domains plus ten pair relations per spin sector",
        "cross_node_pair_correlations_preserved": True,
        "higher_order_correlations_preserved": False,
        "amplitudes_evaluated": False,
        "joint_keys_enumerated": False,
        "results": [],
    }
    for mode in modes:
        item = {"mode": mode}
        try:
            support = FG.factor_graph_triple(
                initial, a.source, *a.triple, a.Jmax2,
                adjoint=(mode == "adjoint"), limit=a.limit
            )
            item["status"] = "completed"
            item.update(FG.profile(support))
        except RuntimeError as exc:
            item["status"] = "guard_hit"
            item["guard_message"] = str(exc)
        report["results"].append(item)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
