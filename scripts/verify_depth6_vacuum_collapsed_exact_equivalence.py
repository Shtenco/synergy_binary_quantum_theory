#!/usr/bin/env python3
"""Exact-arithmetic comparison of expanded and source-collapsed vacuum formulas.

The recovered production backend is floating point. This script ports its two
collapsed stages using the same exact CG/volume primitives as the direct
expanded evaluator and compares all vacuum adjoint color channels symbolically.
It proves equality of the finite formulas on this input, not the historical
branch-sum identity or any non-vacuum rank statement.
"""
from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "research" / "depth6_recovered_source"
sys.path.insert(0, str(SOURCE))

import bqg_block_apply_invariants as A  # noqa: E402

SPEC = importlib.util.spec_from_file_location(
    "depth6_exact_vacuum_reference",
    ROOT / "scripts" / "depth6_vacuum_loop_reduced_exact_cg.py",
)
assert SPEC and SPEC.loader
E = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = E
SPEC.loader.exec_module(E)

EXPECTED_BUNDLE = "2a42c6cdbf82b3b65a2ada02116e31cc9dd30c514589bb2885992050d1d16d32"


def verify_source_hashes() -> dict[str, str]:
    manifest = json.loads((SOURCE / "source_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("bundle_sha256") != EXPECTED_BUNDLE:
        raise RuntimeError("recovered source bundle hash mismatch")
    actual = {name: hashlib.sha256((SOURCE / name).read_bytes()).hexdigest()
              for name in manifest["module_sha256"]}
    if actual != manifest["module_sha256"]:
        raise RuntimeError("recovered source module hashes mismatch")
    return actual


def endpoint_kind_color(x: int, y: int, ci: int, cj: int, endpoint: int):
    if x < y:
        return ("cg", ci if endpoint == x else cj)
    return ("adj", ci if endpoint == x else cj)


def local_two_hit(v, spins, leg1, so1, kind1, color1, leg2, so2, kind2, color2, label):
    out = {}
    for l1, a1 in E.local_hit_outputs(v, tuple(spins), leg1, so1, kind1, color1, label):
        sp1 = list(spins); sp1[leg1] = so1; sp1 = tuple(sp1)
        for l2, a2 in E.local_hit_outputs(v, sp1, leg2, so2, kind2, color2, l1):
            out[l2] = sp.simplify(out.get(l2, 0) + a1 * a2)
    return {k: sp.simplify(v) for k, v in out.items() if sp.simplify(v) != 0}


def local_two_hit_shared(v, spins, leg1, so1, kind1, leg2, so2, kind2, label):
    out = {}
    for m in range(2):
        for lab, value in local_two_hit(v, spins, leg1, so1, kind1, m, leg2, so2, kind2, m, label).items():
            out[lab] = sp.simplify(out.get(lab, 0) + value)
    return {k: v for k, v in out.items() if v != 0}


def exact_edge_norm(s: int, so: int, x: int, y: int) -> sp.Expr:
    return (sp.sqrt(sp.Rational(s + 1, so + 1)) if x < y
            else sp.sqrt(sp.Rational(so + 1, s + 1)))


def collapsed_triangle_path(state, path, iout, iin, jmax2):
    v, a, b, v2 = path
    assert v == v2 and len({v, a, b}) == 3
    R = A.R
    e1 = R.EIDX[tuple(sorted((v, a)))]
    e2 = R.EIDX[tuple(sorted((a, b)))]
    e3 = R.EIDX[tuple(sorted((b, v)))]
    out = {}
    for (spins, locs), amp in state.items():
        s1, s2, s3 = spins[e1], spins[e2], spins[e3]
        for o1 in (s1 - 1, s1 + 1):
            if o1 < 0 or o1 > jmax2: continue
            for o2 in (s2 - 1, s2 + 1):
                if o2 < 0 or o2 > jmax2: continue
                for o3 in (s3 - 1, s3 + 1):
                    if o3 < 0 or o3 > jmax2: continue
                    norm = (exact_edge_norm(s1, o1, v, a)
                            * exact_edge_norm(s2, o2, a, b)
                            * exact_edge_norm(s3, o3, b, v))

                    lsv = R.local_spins(spins, v)
                    lv1, lv3 = R.LEGIDX[(v, a)], R.LEGIDX[(v, b)]
                    k1, c1 = endpoint_kind_color(v, a, iout, 0, v)
                    k3, c3 = endpoint_kind_color(b, v, 0, iin, v)
                    ov = local_two_hit(v, lsv, lv1, o1, k1, c1, lv3, o3, k3, c3, locs[v])
                    if not ov: continue

                    lsa = R.local_spins(spins, a)
                    la1, la2 = R.LEGIDX[(a, v)], R.LEGIDX[(a, b)]
                    ka1, _ = endpoint_kind_color(v, a, 0, 0, a)
                    ka2, _ = endpoint_kind_color(a, b, 0, 0, a)
                    oa = local_two_hit_shared(a, lsa, la1, o1, ka1, la2, o2, ka2, locs[a])
                    if not oa: continue

                    lsb = R.local_spins(spins, b)
                    lb2, lb3 = R.LEGIDX[(b, a)], R.LEGIDX[(b, v)]
                    kb2, _ = endpoint_kind_color(a, b, 0, 0, b)
                    kb3, _ = endpoint_kind_color(b, v, 0, 0, b)
                    ob = local_two_hit_shared(b, lsb, lb2, o2, kb2, lb3, o3, kb3, locs[b])
                    if not ob: continue

                    spin_out = list(spins); spin_out[e1] = o1; spin_out[e2] = o2; spin_out[e3] = o3
                    for lvo, av in ov.items():
                        for lao, aa in oa.items():
                            for lbo, ab in ob.items():
                                loc_out = list(locs); loc_out[v] = lvo; loc_out[a] = lao; loc_out[b] = lbo
                                key = (tuple(spin_out), tuple(loc_out))
                                out[key] = sp.simplify(out.get(key, 0) + amp * norm * av * aa * ab)
    return {k: sp.simplify(v) for k, v in out.items() if sp.simplify(v) != 0}


def collapsed_triangle_adjoint(state, v, a, b, i, j, jmax2):
    out = {}
    for path, sign in (((v, b, a, v), +1), ((v, a, b, v), -1)):
        part = collapsed_triangle_path(state, path, j, i, jmax2)
        for key, value in part.items():
            out[key] = sp.simplify(out.get(key, 0) + sign * value)
    return {k: v for k, v in out.items() if v != 0}


def exact_hit_matrices(sin: int, sout: int, x: int, y: int, i: int, j: int):
    if x < y:
        left = sp.Matrix(sout + 1, sin + 1,
            lambda r, c: E.cg(sin, 1, sout, E.mvals(sin)[c], 1 if i == 0 else -1, E.mvals(sout)[r]))
        right = sp.Matrix(sout + 1, sin + 1,
            lambda r, c: E.cg(sin, 1, sout, E.mvals(sin)[c], 1 if j == 0 else -1, E.mvals(sout)[r]))
        norm = sp.sqrt(sp.Rational(sin + 1, sout + 1))
    else:
        left0 = sp.Matrix(sin + 1, sout + 1,
            lambda r, c: E.cg(sout, 1, sin, E.mvals(sout)[c], 1 if j == 0 else -1, E.mvals(sin)[r]))
        right0 = sp.Matrix(sin + 1, sout + 1,
            lambda r, c: E.cg(sout, 1, sin, E.mvals(sout)[c], 1 if i == 0 else -1, E.mvals(sin)[r]))
        left, right = left0.conjugate().T, right0.conjugate().T
        norm = sp.sqrt(sp.Rational(sout + 1, sin + 1))
    return left, right, norm


def exact_far_radial_weight(s: int, so: int, v: int, c: int):
    W = None
    for k in range(2):
        left1, right1, n1 = exact_hit_matrices(s, so, c, v, k, 0)
        u1, _ = sorted((c, v)); A1 = left1 if c == u1 else right1
        left2, right2, n2 = exact_hit_matrices(so, s, v, c, 0, k)
        # The pinned source contracts on the remote endpoint c here.
        u2, _ = sorted((v, c)); B2 = left2 if c == u2 else right2
        term = n1 * n2 * B2 * A1
        W = term if W is None else W + term
    scalar = sp.simplify(sp.trace(W) / (s + 1))
    if sp.simplify(W - scalar * sp.eye(s + 1)) != sp.zeros(s + 1):
        raise ArithmeticError(("exact remote radial contraction is not scalar", s, so, v, c, W))
    return scalar


def volume_local_outputs(spins_local, label):
    K, L, J, M = label
    if tuple(sorted(spins_local)) != (0, 1, 1, 1):
        return ()
    labels = E.multiplicity_labels(tuple(spins_local), J)
    if (K, L) not in labels: return ()
    if J == 1:
        value = sp.Pow(3, sp.Rational(1, 4)) / 2
    elif J == 3:
        value = sp.Integer(0)
    else:
        raise ArithmeticError(("unexpected exact vacuum volume sector", spins_local, label))
    if value == 0: return ()
    return tuple(((Ko, Lo, J, M), value if idx == labels.index((K, L)) else sp.Integer(0))
                 for idx, (Ko, Lo) in enumerate(labels) if idx == labels.index((K, L)))


def exact_radial_local_matrix(v, c, spins_local, i, j):
    spins_local = tuple(spins_local); R = A.R
    leg = R.LEGIDX[(v, c)]; s = spins_local[leg]
    labels = E.allowed_total_J(spins_local)
    all_labels = tuple((K, L, J, M) for J in labels for M in E.mvals(J)
                       for K, L in E.multiplicity_labels(spins_local, J))
    idx = {lab: n for n, lab in enumerate(all_labels)}
    matrix = sp.zeros(len(all_labels), len(all_labels))
    kind1 = "cg" if c < v else "adj"
    kind2 = "cg" if v < c else "adj"
    for col, lab in enumerate(all_labels):
        # V then edge pair
        for lv, cv in volume_local_outputs(spins_local, lab):
            for so in (s - 1, s + 1):
                if so < 0: continue
                sp1 = list(spins_local); sp1[leg] = so; sp1 = tuple(sp1)
                weight = exact_far_radial_weight(s, so, v, c)
                for l1, a1 in E.local_hit_outputs(v, spins_local, leg, so, kind1, i, lv):
                    for l2, a2 in E.local_hit_outputs(v, sp1, leg, s, kind2, j, l1):
                        matrix[idx[l2], col] += weight * cv * a1 * a2
        # edge pair with V between hits
        for so in (s - 1, s + 1):
            if so < 0: continue
            sp1 = list(spins_local); sp1[leg] = so; sp1 = tuple(sp1)
            weight = exact_far_radial_weight(s, so, v, c)
            for l1, a1 in E.local_hit_outputs(v, spins_local, leg, so, kind1, i, lab):
                for lv, cv in volume_local_outputs(sp1, l1):
                    for l2, a2 in E.local_hit_outputs(v, sp1, leg, s, kind2, j, lv):
                        matrix[idx[l2], col] -= weight * a1 * cv * a2
    return all_labels, matrix.applyfunc(sp.simplify)


def collapsed_radial_adjoint(state, v, c, i, j):
    out = {}
    for (spins, locs), amp in state.items():
        local_spins = A.R.local_spins(spins, v)
        labels, matrix = exact_radial_local_matrix(v, c, local_spins, i, j)
        label = locs[v]; row = labels.index(label)
        for out_index, value in enumerate(matrix.row(row)):
            value = sp.conjugate(value)
            if value == 0: continue
            loc_out = list(locs); loc_out[v] = labels[out_index]
            key = (spins, tuple(loc_out))
            out[key] = sp.simplify(out.get(key, 0) + amp * value)
    return {k: sp.simplify(v) for k, v in out.items() if sp.simplify(v) != 0}


def assert_same(a, b, label):
    keys = set(a) | set(b)
    bad = [(key, a.get(key, 0), b.get(key, 0), sp.simplify(a.get(key, 0) - b.get(key, 0)))
           for key in keys if sp.simplify(a.get(key, 0) - b.get(key, 0)) != 0]
    if bad:
        raise AssertionError((label, bad[:3], len(bad)))
    return len(keys)


def run():
    hashes = verify_source_hashes()
    vacuum = {( (0,) * 10, tuple((0, 0, 0, 0) for _ in range(5)) ): sp.Integer(1)}
    node_results = []
    entry_comparisons = 0
    radial_exact_weights = {}
    for node in range(5):
        artifact = json.loads((ROOT / "verification_results" / f"DEPTH6_5_VACUUM_ADJOINT_EXACT_NODE{node}_20261126.json").read_text(encoding="utf-8"))
        specs = A.R.oriented_specs(node)
        term_total = 0; triangle_state_count = 0; radial_state_count = 0
        for (sign, spec), term in zip(specs, artifact["terms"], strict=True):
            v, a, b, c = spec
            exact_channel_mats = {
                # artifact keys are repr(tuple); ast is imported lazily to avoid altering arithmetic
                __import__("ast").literal_eval(key): sp.sympify(value, locals={"Matrix": sp.Matrix})
                for key, value in term["adjoint_channel_matrices_by_output_state"].items()
            }
            for i, j in itertools.product(range(2), repeat=2):
                exact_loop = E.triangle_adjoint(i, j, v, a, b)
                collapsed_loop = collapsed_triangle_adjoint(vacuum, v, a, b, i, j, 1)
                triangle_state_count += assert_same(exact_loop, collapsed_loop,
                    ("triangle", node, spec, i, j))

                direct = E.radial_adjoint(exact_loop, v, c, i, j)
                collapsed = collapsed_radial_adjoint(collapsed_loop, v, c, i, j)
                radial_state_count += assert_same(direct, collapsed,
                    ("radial", node, spec, i, j))

                keys = set(exact_channel_mats) | set(direct) | set(collapsed)
                for key in keys:
                    expected = exact_channel_mats.get(key)
                    ev = expected[i, j] if expected is not None else sp.Integer(0)
                    if sp.simplify(ev - direct.get(key, 0)) != 0:
                        raise AssertionError(("artifact/reference", node, spec, i, j, key))
                term_total += len(keys)
                entry_comparisons += len(keys)
                s = 0
                so = 1
                radial_exact_weights[f"{v},{c}:{i},{j}"] = str(exact_far_radial_weight(s, so, v, c))
        node_results.append({
            "node": node,
            "oriented_terms": len(specs),
            "triangle_exact_state_keys_compared": triangle_state_count,
            "radial_exact_state_keys_compared": radial_state_count,
            "channel_state_entries_compared": term_total,
        })
    return {
        "status": "EXACT_FINITE_COLLAPSED_EQUIVALENCE_ON_VACUUM",
        "source": {
            "commit": "f6bfcd775b8c3cd33c6c235d48293e40657baf9d",
            "bundle_sha256": EXPECTED_BUNDLE,
            "module_sha256": hashes,
            "ported_stages": ["triangle_adjoint_collapsed", "radial_local_matrix", "far_radial_weight", "local_hit_outputs"],
        },
        "reference": "scripts/depth6_vacuum_loop_reduced_exact_cg.py with exact SymPy CG coefficients and exact local volume block",
        "input_scope": "all-zero vacuum; five nodes; every one of 12 oriented terms; all four SU(2) index pairs",
        "exact_identity_checks": {
            "expanded_vs_collapsed_triangle_adjoint": "symbolically equal for every compared vacuum intermediate state",
            "expanded_vs_collapsed_radial_adjoint": "symbolically equal for every compared vacuum intermediate state",
            "exact_endpoint_far_radial_weights": radial_exact_weights,
            "state_channel_entries": entry_comparisons,
            "per_node": node_results,
        },
        "claim_limit": "This proves exact equality of the ported recovered collapsed stages and the direct expanded reduced-basis formula on the vacuum column, under the pinned source formulas. The actual recovered Python backend still evaluates CG/recoupling in floating point; [5] selector integration and historical B_[5] identity remain open, and no non-vacuum rank is implied.",
    }


if __name__ == "__main__":
    result = run()
    print(json.dumps(result, indent=2))
