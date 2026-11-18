#!/usr/bin/env python3
"""Candidate global correspondence: PL tetrahedra -> active Q2 causal edges.

The map is many-to-one. A PL tetrahedron is assigned to one active causal
edge; every PL-dual adjacency must map to two causal edges sharing an endpoint.
Route labels are selected by the colour of the parent face opposite the first
vertex in a barycentric flag. The L/R child is selected by the side assigned
to the parent face opposite the final vertex.

For a shared parent face:
* if its two parent tetrahedra map to distinct causal edges, those edges must
  share exactly one endpoint; each side is forced by which endpoint is shared;
* if they map to the same causal edge, the oriented PL boundary fixes
  opposite sides on the two incident tetrahedra.

This is a finite combinatorial correspondence gate, not a metric, QFI, or
continuum-limit derivation. It is audited only through PL generation 2.
"""
from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from q2_pl_graph_export import (
    barycentric_subdivision,
    cross_polytope_boundary_4,
    dual_graph_payload,
)
from q2_causal_graph_export import ROUTE_LABELS


def _permutation_data(chain, faces, parent):
    """Recover the ordered parent-vertex permutation from a barycentric flag."""
    permutation = []
    previous = set()
    for face_id in chain:
        current = set(faces[face_id])
        added = current - previous
        if len(added) != 1:
            raise ValueError("barycentric flag is not a maximal simplex chain")
        permutation.append(next(iter(added)))
        previous = current
    if len(permutation) != 4 or set(permutation) != set(parent):
        raise ValueError("flag does not enumerate its parent tetrahedron")
    return permutation


def _vertex_face_colours(dual):
    colours = [dict() for _ in dual["vertices"]]
    for edge in dual["edges"]:
        face = set(edge["shared_face"])
        for endpoint in (edge["source"], edge["target"]):
            tetrahedron = set(dual["vertices"][endpoint]["tetrahedron"])
            omitted = tetrahedron - face
            if len(omitted) != 1:
                raise ValueError("dual edge does not identify one opposite vertex")
            vertex = next(iter(omitted))
            old = colours[endpoint].get(vertex)
            if old is not None and old != edge["route_colour"]:
                raise ValueError("inconsistent route colour on a tetrahedral face")
            colours[endpoint][vertex] = edge["route_colour"]
    if any(sorted(mapping.values()) != [0, 1, 2, 3] for mapping in colours):
        raise ValueError("each tetrahedron must have one face of every route colour")
    return colours


def _boundary_side(tetrahedron, face, orientation):
    omitted = next(iter(set(tetrahedron) - set(face)))
    position = tuple(sorted(tetrahedron)).index(omitted)
    sign = orientation * (-1 if position % 2 else 1)
    return "L" if sign > 0 else "R"


def _parent_face_sides(tetrahedra, parent_map, causal_edges, orientations, dual):
    """Derive L/R for every oriented parent face; return diagnostics too."""
    sides = {}
    same_edge_faces = 0
    distinct_edge_faces = 0
    for edge in dual["edges"]:
        a, b = edge["source"], edge["target"]
        ta = tuple(tetrahedra[a])
        tb = tuple(tetrahedra[b])
        face = tuple(sorted(edge["shared_face"]))
        ea, eb = parent_map[ta], parent_map[tb]
        if ea == eb:
            same_edge_faces += 1
            side_a = _boundary_side(ta, face, orientations[ta])
            side_b = _boundary_side(tb, face, orientations[tb])
            if side_a == side_b:
                raise ValueError("oriented PL boundary sides do not cancel")
        else:
            distinct_edge_faces += 1
            ua, va = causal_edges[ea]
            ub, vb = causal_edges[eb]
            common = set((ua, va)) & set((ub, vb))
            if len(common) != 1:
                raise ValueError("adjacent PL parents map to non-incident causal edges")
            shared = next(iter(common))
            side_a = "L" if ua == shared else "R"
            side_b = "L" if ub == shared else "R"
        sides[(ta, face)] = side_a
        sides[(tb, face)] = side_b
    return sides, {
        "same_parent_edge_faces": same_edge_faces,
        "distinct_incident_parent_edge_faces": distinct_edge_faces,
        "face_side_assignments": len(sides),
    }


def _causal_children(causal_edges, generation, next_vertex):
    children = {}
    for parent_id, (source, target) in causal_edges.items():
        for label in ROUTE_LABELS:
            midpoint = next_vertex
            next_vertex += 1
            children[f"g{generation}:{parent_id}:{label}:L"] = (source, midpoint)
            children[f"g{generation}:{parent_id}:{label}:R"] = (midpoint, target)
    return children, next_vertex


def _summarise(generation, tetrahedra, parent_map, causal_edges, dual):
    fibres = Counter(parent_map[tuple(tet)] for tet in tetrahedra)
    adjacency_ok = True
    for edge in dual["edges"]:
        a = tuple(dual["vertices"][edge["source"]]["tetrahedron"])
        b = tuple(dual["vertices"][edge["target"]]["tetrahedron"])
        ea, eb = parent_map[a], parent_map[b]
        if not (set(causal_edges[ea]) & set(causal_edges[eb])):
            adjacency_ok = False
            break
    histogram = Counter(fibres.values())
    return {
        "generation": generation,
        "pl_tetrahedra": len(tetrahedra),
        "causal_active_edges": len(causal_edges),
        "surjective": set(fibres) == set(causal_edges),
        "fibre_min": min(fibres.values()),
        "fibre_max": max(fibres.values()),
        "fibre_mean": len(tetrahedra) / len(causal_edges),
        "fibre_size_histogram": {str(k): histogram[k] for k in sorted(histogram)},
        "pl_dual_adjacencies": len(dual["edges"]),
        "adjacency_preserved": adjacency_ok,
    }


def run(max_generation=2):
    if not 0 <= max_generation <= 2:
        raise ValueError("audited PL exporter supports generations 0..2")

    tetrahedra = cross_polytope_boundary_4()
    parent_map = {tuple(tet): "e0" for tet in tetrahedra}
    orientations = {
        tuple(tet): (-1 if sum(vertex % 2 for vertex in tet) % 2 else 1)
        for tet in tetrahedra
    }
    causal_edges = {"e0": (0, 1)}
    next_vertex = 2
    previous_faces = None
    chain_by_tet = None
    levels = []
    refinement_checks = []
    face_diagnostics = []

    for generation in range(max_generation + 1):
        dual = dual_graph_payload(
            generation, tetrahedra, chain_by_tet, previous_faces
        )
        _vertex_face_colours(dual)
        levels.append(_summarise(
            generation, tetrahedra, parent_map, causal_edges, dual
        ))
        if generation == max_generation:
            break

        vertex_colours = _vertex_face_colours(dual)
        face_sides, face_info = _parent_face_sides(
            tetrahedra, parent_map, causal_edges, orientations, dual
        )
        face_diagnostics.append({
            "generation": generation,
            **face_info,
        })

        children, faces, child_chains = barycentric_subdivision(tetrahedra)
        parent_index = {tuple(tet): i for i, tet in enumerate(tetrahedra)}
        child_causal_edges, next_vertex = _causal_children(
            causal_edges, generation + 1, next_vertex
        )
        child_map = {}
        child_orientations = {}
        naturality_ok = True

        for child in children:
            chain = child_chains[tuple(child)]
            parent = tuple(faces[chain[-1]])
            pidx = parent_index[parent]
            permutation = _permutation_data(chain, faces, parent)
            route_colour = vertex_colours[pidx][permutation[0]]
            omitted_last = permutation[-1]
            parent_face = tuple(sorted(v for v in parent if v != omitted_last))
            side = face_sides[(parent, parent_face)]
            label = ROUTE_LABELS[route_colour]
            causal_id = f"g{generation + 1}:{parent_map[parent]}:{label}:{side}"
            child_map[tuple(child)] = causal_id
            if causal_id not in child_causal_edges:
                naturality_ok = False

            positions = [tuple(sorted(parent)).index(v) for v in permutation]
            inversions = sum(
                positions[i] > positions[j]
                for i in range(4) for j in range(i + 1, 4)
            )
            child_orientations[tuple(child)] = (
                orientations[parent] * (-1 if inversions % 2 else 1)
            )

        child_dual = dual_graph_payload(
            generation + 1, children, child_chains, faces
        )
        child_level = _summarise(
            generation + 1, children, child_map, child_causal_edges, child_dual
        )
        refinement_checks.append({
            "generation": generation + 1,
            "every_child_maps_to_a_child_of_its_parent": naturality_ok,
        })
        # The next loop recomputes the same level summary from the carried map.
        tetrahedra = children
        parent_map = child_map
        orientations = child_orientations
        causal_edges = child_causal_edges
        previous_faces = faces
        chain_by_tet = child_chains

    checks = {
        "all_levels_surjective": all(row["surjective"] for row in levels),
        "all_levels_preserve_dual_adjacency": all(
            row["adjacency_preserved"] for row in levels
        ),
        "all_fibres_nonempty": all(row["fibre_min"] > 0 for row in levels),
        "refinement_naturality": all(
            row["every_child_maps_to_a_child_of_its_parent"]
            for row in refinement_checks
        ),
        "all_parent_face_side_assignments_valid": True,
        "four_route_face_colours_present": True,
    }
    return {
        "status": "candidate q2 causal-to-PL edge correspondence",
        "passed": all(checks.values()),
        "map_semantics": (
            "PL tetrahedra map many-to-one onto active causal edges; every PL "
            "dual adjacency maps to endpoint incidence of causal edges"
        ),
        "side_rule": (
            "distinct parent causal edges: choose the child side containing "
            "their shared endpoint; equal parent edge: use opposite oriented "
            "PL boundary sides"
        ),
        "levels": levels,
        "parent_face_diagnostics": face_diagnostics,
        "refinement_checks": refinement_checks,
        "checks": checks,
        "claim_boundary": (
            "Finite combinatorial graph correspondence only. Nonuniform fibres "
            "at generation 2 are reported, not hidden or reweighted. This does "
            "not derive a physical measure, QFI conductance, spectral dimension, "
            "semiclassical dynamics, or spacetime."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-generation", type=int, default=2)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args.max_generation)
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
