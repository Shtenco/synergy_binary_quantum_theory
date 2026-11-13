#!/usr/bin/env python3
"""BQG V16 full initial 32-state K5 certificate.

Runs against checkpointed H_v source columns from the canonical V15 runner.
Refuses false zeros for incomplete columns and refuses noncanonical upstream
unless explicitly asked for a research-only diagnostic.
"""
from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np

OFFICIAL_BLOB_SHA = "e7e1e0705c4b833dcee167a2036ef5f8929a4388"

def git_blob_sha(data):
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()

def parse_column(rows):
    out = {}
    for x in rows:
        k = (tuple(x["key"][0]), tuple(x["key"][1]))
        z = x.get("amplitude", x.get("coefficient"))
        if z is None:
            raise ValueError("missing amplitude")
        amp = complex(float.fromhex(z[0]), float.fromhex(z[1]))
        if not np.isfinite(amp.real) or not np.isfinite(amp.imag):
            raise ValueError("nonfinite coefficient")
        if k in out:
            raise ValueError("duplicate output basis state")
        out[k] = amp
    return out

def certify(backend, cache_dirs, allow_noncanonical=False):
    backend = Path(backend)
    raw = backend.read_bytes()
    git_sha = git_blob_sha(raw)
    sha256 = hashlib.sha256(raw).hexdigest()
    if git_sha != OFFICIAL_BLOB_SHA and not allow_noncanonical:
        raise ValueError(f"not canonical GitHub K5 blob: {git_sha}")
    sys.path.insert(0, str(backend.parent))
    import k5_peter_weyl_safe_hda_column as K5
    sources = K5.basis_full_jhalf()
    if len(sources) != 32:
        raise ValueError("expected 32 initial K5 basis states")
    cols = {}
    outputs = set(sources)
    for i, source in enumerate(sources):
        for v in K5.VERT:
            candidates = []
            for root in map(Path, cache_dirs):
                candidates += [
                    root / f"source_{i:02d}_node_{v}_wall5.json",
                    root / f"source_{i:02d}_vertex_{v}_wall5.json",
                ]
            hits = [p for p in candidates if p.exists()]
            if len(hits) != 1:
                raise RuntimeError(f"INCOMPLETE or ambiguous source={i} H_v={v}; do not insert zero")
            d = json.loads(hits[0].read_text())
            if (d["backend_sha256"] != sha256 or
                d["input"] != [list(source[0]), list(source[1])] or
                d["node"] != v or d["wall"] != 5):
                raise ValueError(f"column provenance error: {hits[0]}")
            vals = parse_column(d["rows"])
            cols[i, v] = vals
            outputs.update(vals)
    if len(cols) != 160:
        raise RuntimeError("only 160/160 columns establish the 32 input Gram")
    outkeys = sorted(outputs, key=repr)
    outidx = {k: i for i, k in enumerate(outkeys)}
    C = np.zeros((5 * len(outkeys), 32), complex)
    sparse = []
    inside = 0.0
    changed = 0.0
    for i in range(32):
        di = {}
        for v in range(5):
            for k, a in cols[i, v].items():
                C[v * len(outkeys) + outidx[k], i] = a
                di[(v, k)] = a
                if k in sources: inside += abs(a)**2
                if 0 in k[0]: changed += abs(a)**2
        sparse.append(di)
    direct = C.conj().T @ C
    independent = np.zeros((32, 32), complex)
    for j in range(32):
        for i in range(j + 1):
            z = sum(a.conjugate() * sparse[j].get(k, 0j)
                    for k, a in sparse[i].items())
            independent[i, j] = z
            independent[j, i] = z.conjugate()
    error = float(np.linalg.norm(independent - direct))
    if error > 5e-10:
        raise AssertionError(f"two Gram assembly methods disagree: {error}")
    eig = np.linalg.eigvalsh(direct)
    sv = np.linalg.svd(C, compute_uv=False)
    r = int(np.count_nonzero(sv > 1e-10 * max(float(sv[0]), 1.0)))
    trace = float(np.trace(direct).real)
    return {
        "scope": "finite Euclidean 5-H_v K5 initial 32 sector; NOT full physical projector",
        "github_blob_sha": git_sha,
        "backend_sha256": sha256,
        "canonical_blob_verified": git_sha == OFFICIAL_BLOB_SHA,
        "columns_complete": len(cols),
        "sources_complete": len(sources),
        "original_input_graph_sector_count": 1,
        "distinct_output_keys_including_32_original": len(outkeys),
        "trace": trace, "inside_original32_norm2": inside,
        "changed_graph_norm2_fraction": changed / trace if trace else 0.0,
        "smallest_eigenvalue": float(eig[0]),
        "largest_eigenvalue": float(eig[-1]),
        "minimum_singular_value": float(sv[-1]),
        "rank": r,
        "numerical_nullity": 32-r,
        "independent_gram_difference": error,
        "eigenvalues": eig.tolist(),
        "physical_projector_computed": False,
        "cross_graph_physical_amplitude_computed": False,
    }, direct

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--backend", type=Path, default=Path("scripts/k5_peter_weyl_safe_hda_column.py"))
    p.add_argument("--cache-dirs", type=Path, nargs="+", required=True)
    p.add_argument("--allow-noncanonical", action="store_true")
    p.add_argument("--out-dir", type=Path, required=True)
    a = p.parse_args()
    result, M = certify(a.backend, a.cache_dirs, a.allow_noncanonical)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    (a.out_dir/"V16_CERTIFICATE.json").write_text(json.dumps(result, indent=2) + "\n")
    np.savez_compressed(a.out_dir/"V16_GRAM_32.npz", M=M)
    print(json.dumps({k:v for k,v in result.items() if k != "eigenvalues"}, indent=2))
if __name__ == "__main__":
    main()
