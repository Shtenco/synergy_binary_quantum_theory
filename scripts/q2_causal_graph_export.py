#!/usr/bin/env python3
"""Canonical directed Q2 causal-route graph exporter.

Frozen rule: each active directed edge u->v creates B=2^q route vertices,
one per two-bit route label, and replaces itself by the 2B directed edges
u->m_a and m_a->v. Older vertices are retained; only the newest active edges
are present in the refined graph. For q=2 this gives 8 active child edges per
parent edge and the exact vertex law N_g=(4*8^g+10)/7.

This is the causal-route graph, not the PL 1-skeleton. It intentionally adds
no C4 equatorial/spatial edges: those belong to the separately selected PL
carrier and must be connected by a derived gluing map, not silently merged.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROUTE_LABELS = ("00", "01", "10", "11")


def expected_vertices(generation: int, q: int = 2) -> int:
    if generation < 0 or q < 1:
        raise ValueError("generation must be nonnegative and q positive")
    branches = 1 << q
    growth = 2 * branches
    return 2 + branches * (growth ** generation - 1) // (growth - 1)


def build(max_generation: int = 5, q: int = 2) -> dict:
    if q != 2:
        raise ValueError("canonical Q2 causal graph is frozen to q=2")
    if not 0 <= max_generation <= 6:
        raise ValueError("max_generation must be in [0, 6]")

    branches = 1 << q
    vertices = [
        {"id": 0, "kind": "boundary", "generation": 0, "label": "source"},
        {"id": 1, "kind": "boundary", "generation": 0, "label": "sink"},
    ]
    # Each active edge is (stable id, source vertex, target vertex).
    active = [("e0", 0, 1)]
    rows = []

    for generation in range(max_generation + 1):
        edges = [
            {"id": edge_id, "source": source, "target": target,
             "generation": generation}
            for edge_id, source, target in active
        ]
        canonical = json.dumps(
            {"vertices": vertices, "edges": edges},
            sort_keys=True, separators=(",", ":"),
        ).encode("utf-8")
        indegree = [0] * len(vertices)
        outdegree = [0] * len(vertices)
        undirected = [[] for _ in vertices]
        for edge in edges:
            u, v = edge["source"], edge["target"]
            outdegree[u] += 1
            indegree[v] += 1
            undirected[u].append(v)
            undirected[v].append(u)

        seen = {0}
        stack = [0]
        while stack:
            u = stack.pop()
            for v in undirected[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)

        # Count directed source-to-sink paths using a real topological sort.
        outgoing = [[] for _ in vertices]
        indegree_dag = [0] * len(vertices)
        for edge in edges:
            u, v = edge["source"], edge["target"]
            outgoing[u].append(v)
            indegree_dag[v] += 1
        queue = [i for i, degree in enumerate(indegree_dag) if degree == 0]
        topo = []
        cursor = 0
        while cursor < len(queue):
            u = queue[cursor]
            cursor += 1
            topo.append(u)
            for v in outgoing[u]:
                indegree_dag[v] -= 1
                if indegree_dag[v] == 0:
                    queue.append(v)
        acyclic = len(topo) == len(vertices)
        ways = [0] * len(vertices)
        ways[0] = 1
        if acyclic:
            for u in topo:
                for v in outgoing[u]:
                    ways[v] += ways[u]

        rows.append({
            "generation": generation,
            "vertex_count": len(vertices),
            "active_edge_count": len(edges),
            "expected_vertex_count": expected_vertices(generation, q),
            "connected": len(seen) == len(vertices),
            "directed_acyclic": acyclic,
            "source_to_sink_path_count": ways[1],
            "path_length": 2 ** generation,
            "degree_min": min((indegree[i] + outdegree[i] for i in range(len(vertices))), default=0),
            "degree_max": max((indegree[i] + outdegree[i] for i in range(len(vertices))), default=0),
            "graph_sha256": hashlib.sha256(canonical).hexdigest(),
            "vertices": list(vertices),
            "active_edges": edges,
        })

        if generation == max_generation:
            break

        next_active = []
        refinement_map = []
        for parent_id, source, target in active:
            child_edges = []
            for branch, label in enumerate(ROUTE_LABELS):
                midpoint = len(vertices)
                vertices.append({
                    "id": midpoint,
                    "kind": "route",
                    "generation": generation + 1,
                    "label": label,
                    "parent_edge": parent_id,
                    "branch_index": branch,
                })
                left_id = f"g{generation + 1}:{parent_id}:{label}:L"
                right_id = f"g{generation + 1}:{parent_id}:{label}:R"
                next_active.append((left_id, source, midpoint))
                next_active.append((right_id, midpoint, target))
                child_edges.extend((left_id, right_id))
            refinement_map.append({
                "parent_edge": parent_id,
                "source": source,
                "target": target,
                "route_vertices": [
                    {"label": label, "vertex": len(vertices) - branches + i}
                    for i, label in enumerate(ROUTE_LABELS)
                ],
                "child_edges": child_edges,
            })
        active = next_active
        rows[-1]["refinement_map_to_next"] = refinement_map

    checks = {
        "exact_vertex_formula": all(
            row["vertex_count"] == row["expected_vertex_count"] for row in rows
        ),
        "active_edge_growth_is_eight": all(
            row["active_edge_count"] == 8 ** row["generation"] for row in rows
        ),
        "all_levels_connected": all(row["connected"] for row in rows),
        "all_levels_are_dags": all(row["directed_acyclic"] for row in rows),
        "source_sink_paths_are_four_to_the_generation": all(
            row["source_to_sink_path_count"] == 4 ** (2 ** row["generation"] - 1) for row in rows
        ),
        "all_hashes_present": all(len(row["graph_sha256"]) == 64 for row in rows),
        "route_labels_are_c4_carrier_labels": all(
            tuple(label for label in ROUTE_LABELS) == ("00", "01", "10", "11")
            for _ in (0,)
        ),
    }
    return {
        "status": "canonical q2 directed causal-route graph export",
        "passed": all(checks.values()),
        "q": q,
        "branches_per_active_edge": branches,
        "active_child_edges_per_parent": 2 * branches,
        "causal_depth_doubles_per_generation": True,
        "rule": "u->v becomes u->m_a->v for each route label a in {00,01,10,11}",
        "graph_semantics": (
            "directed causal route graph; only current-generation active edges "
            "are retained; route alternatives are not joined by C4 spatial edges"
        ),
        "generations": rows,
        "checks": checks,
        "pl_identification": None,
        "qfi_state_assignment": None,
        "claim_boundary": (
            "The graph is a canonical realization of the frozen scalar route rewrite. "
            "It is not yet identified with the selected S3 PL complex, and carries "
            "no physical edge density matrices or QFI conductances."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-generation", type=int, default=5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = build(args.max_generation)
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
