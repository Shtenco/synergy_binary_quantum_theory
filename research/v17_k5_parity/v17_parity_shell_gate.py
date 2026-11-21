#!/usr/bin/env python3
"""V17 parity, graph-sector and restricted master audit from real saved K5 columns.

Does not fill missing columns with zeros. Canonical Git blob SHA is required
unless --allow-noncanonical is specified for research-only diagnostics.
"""
from __future__ import annotations
import argparse, collections, hashlib, itertools, json, sys
from pathlib import Path
import numpy as np

OFFICIAL_K5_BLOB = "e7e1e0705c4b833dcee167a2036ef5f8929a4388"
EDGES = tuple(itertools.combinations(range(5), 2))

def gitblob(data):
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()

def kread(x):
    return tuple(x[0]), tuple(x[1])

def zread(x):
    return complex(float.fromhex(x[0]), float.fromhex(x[1]))

def signature(k):
    return tuple(EDGES[i] for i, spin in enumerate(k[0]) if spin != 0)

def get_column(path, k, v, sha256):
    d = json.loads(path.read_text())
    if d["backend_sha256"] != sha256 or d["node"] != v or kread(d["input"]) != k:
        raise ValueError("invalid K5 column provenance: " + str(path))
    result = {}
    for r in d["rows"]:
        target = kread(r["key"])
        if target in result:
            raise ValueError("duplicate output " + repr(target))
        coeff = r.get("coefficient", r.get("amplitude"))
        result[target] = zread(coeff)
    return result

def load_complete_source_set(backend, directory, allow_noncanonical=False, extra_dirs=()):
    raw = Path(backend).read_bytes()
    blob = gitblob(raw)
    if blob != OFFICIAL_K5_BLOB and not allow_noncanonical:
        raise RuntimeError("K5 original Git blob differs: " + blob)
    sha256 = hashlib.sha256(raw).hexdigest()
    sys.path.insert(0, str(Path(backend).parent.resolve()))
    import k5_peter_weyl_safe_hda_column as K5
    origins = K5.basis_full_jhalf()
    if len(origins) != 32:
        raise ValueError("original sector must have 32 states")
    dirs = [Path(directory)] + [Path(x) for x in extra_dirs]
    originals = {}
    for i, k in enumerate(origins):
        for v in range(5):
            candidates = []
            for dr in dirs:
                candidates += [dr / f"source_{i:02d}_node_{v}_wall5.json",
                               dr / f"source_{i:02d}_vertex_{v}_wall5.json"]
            matches = [p for p in candidates if p.is_file()]
            if len(matches) != 1:
                raise RuntimeError(f"incomplete/ambiguous original column {i},{v}")
            originals[i, v] = get_column(matches[0], k, v, sha256)
    return origins, originals, sha256, blob

def check_support(originals):
    first = set().union(*(set(c) for c in originals.values()))
    if len(first) != 720:
        raise AssertionError(f"expected 720 first-shell states, got {len(first)}")
    counts = collections.Counter(sum(s == 0 for s in k[0]) for k in first)
    graphs = collections.Counter(signature(k) for k in first)
    expected = {0: 320, 1: 240, 2: 120, 3: 40}
    if dict(counts) != expected:
        raise AssertionError(("unexpected single-step zero-link histogram", counts))
    by_missing = collections.Counter(10 - len(g) for g in graphs)
    if dict(by_missing) != {0: 1, 1: 10, 2: 30, 3: 10}:
        raise AssertionError(("unexpected sector histogram", by_missing))
    for key in first:
        if sum(key[0]) % 2 != 1:
            raise AssertionError("single H_v must flip doubled-spin parity")
        missing = [set(EDGES[i]) for i, s in enumerate(key[0]) if s == 0]
        if len(missing) == 2 and not missing[0].intersection(missing[1]):
            raise AssertionError("two missing edges are not adjacent")
        if len(missing) == 3 and len(set.union(*missing)) != 3:
            raise AssertionError("missing three edges are not a triangle")
    return first, {"first_shell_states": len(first),
                   "graph_signatures": len(graphs),
                   "state_count_by_absent_links": dict(sorted(counts.items())),
                   "graph_count_by_absent_links": dict(sorted(by_missing.items()))}

def load_selected(path_dirs, sha256):
    by_key = collections.defaultdict(dict)
    for directory in map(Path, path_dirs):
        if not directory.exists():
            continue
        for path in directory.glob("*.json"):
            d = json.loads(path.read_text())
            if not {"input", "node", "rows", "backend_sha256"} <= d.keys():
                continue
            source = kread(d["input"])
            v = int(d["node"])
            col = get_column(path, source, v, sha256)
            if v in by_key[source]:
                raise ValueError("duplicate sampled incoming K5 column")
            by_key[source][v] = col
    return {k: cs for k, cs in by_key.items() if len(cs) == 5}

def restricted_audit(origins, originals, selected):
    old = set(origins)
    selected = {k: cs for k, cs in selected.items() if k not in old}
    inputs = list(origins) + sorted(selected, key=repr)
    n = len(inputs)
    sparse = []
    for i, state in enumerate(inputs):
        by_v = {v: originals[i, v] for v in range(5)} if i < 32 else selected[state]
        sparse.append({(v, k): z for v in range(5) for k, z in by_v[v].items()})
    M = np.zeros((n, n), complex)
    for i in range(n):
        for j in range(i + 1):
            z = sum(a.conjugate() * sparse[j].get(k, 0j) for k, a in sparse[i].items())
            M[i, j] = z
            M[j, i] = z.conjugate()
    eig = np.linalg.eigvalsh(M)
    parity = [sum(k[0]) % 2 for k in inputs]
    opposite = max([abs(M[i, j]) for i in range(n) for j in range(n)
                    if parity[i] != parity[j]] or [0])
    if opposite > 1e-10:
        raise AssertionError("master violated exact doubled-spin parity")
    errors = []
    nontrivial = 0
    for state, cs in selected.items():
        for v in range(5):
            for i, origin in enumerate(origins):
                a = originals[i, v].get(state, 0j)
                b = cs[v].get(origin, 0j)
                errors.append(abs(a - b.conjugate()))
                nontrivial += abs(a) > 1e-12 or abs(b) > 1e-12
    return {
        "input_dimension": n, "selected_first_odd_inputs": parity[32:].count(1),
        "selected_second_even_inputs": parity[32:].count(0),
        "rank": int(np.count_nonzero(eig > max(float(eig[-1]), 1.0) * 1e-10)),
        "min_eigenvalue": float(eig[0]), "max_eigenvalue": float(eig[-1]),
        "opposite_parity_master_max": float(opposite),
        "adjoint_pairs": len(errors), "adjoint_nontrivial_pairs": int(nontrivial),
        "adjoint_max_difference": float(max(errors, default=0.0)),
        "all_first_shell_inputs_computed": False,
        "all_second_shell_inputs_computed": False,
        "full_BQG_physical_projector": False,
    }, M

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", type=Path, default=Path("scripts/k5_peter_weyl_safe_hda_column.py"))
    ap.add_argument("--original-dir", type=Path, required=True)
    ap.add_argument("--more-original-dirs", type=Path, nargs="*", default=[])
    ap.add_argument("--selected-dirs", type=Path, nargs="*", default=[])
    ap.add_argument("--allow-noncanonical", action="store_true")
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    origins, originals, sha256, blob = load_complete_source_set(
        args.backend, args.original_dir, args.allow_noncanonical, args.more_original_dirs)
    first, report = check_support(originals)
    samples = load_selected(args.selected_dirs, sha256)
    nums, M = restricted_audit(origins, originals, samples)
    report.update(nums)
    report["backend_git_blob_sha"] = blob
    report["backend_sha256"] = sha256
    report["canonical_backend_verified"] = blob == OFFICIAL_K5_BLOB
    report["note"] = "Finite selected-input Gram only. No physical BQG projector or 3+1 spacetime proof."
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "V17_RESULT.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    np.savez_compressed(args.output_dir / "V17_GRAM.npz", M=M)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
