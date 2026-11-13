#!/usr/bin/env python3
"""V15: resumable ORIGINAL K5 backend 32-input leakage-safe master audit.

Runs actual scripts/k5_peter_weyl_safe_hda_column.py, never inserts zero
columns for computations that have not completed. JSON/float.hex checkpointing.
A restricted Gram matrix is not the complete BQG physical projector.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE.parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))
import k5_peter_weyl_safe_hda_column as K5

def hash_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(p, data):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    t = p.with_suffix(".tmp")
    t.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    t.replace(p)

def keyout(k):
    return [list(k[0]), list(k[1])]

def keyin(a):
    return tuple(a[0]), tuple(a[1])

def coeffout(a):
    return [float(a.real).hex(), float(a.imag).hex()]

def coeffin(a):
    return complex(float.fromhex(a[0]), float.fromhex(a[1]))

def load_checked(path, source, v, wall, backend_hash):
    d = json.loads(path.read_text())
    if (d["backend_sha256"] != backend_hash or
        d["input"] != keyout(source) or
        d["node"] != v or d["wall"] != wall):
        raise RuntimeError(f"Cache provenance mismatch: {path}")
    return {keyin(r["key"]): coeffin(r["amplitude"]) for r in d["rows"]}

def run(output_dir, wall=5, max_new=160, source_limit=32):
    output_dir = Path(output_dir)
    sources = K5.basis_full_jhalf()
    assert len(sources) == 32
    assert 1 <= source_limit <= 32
    backend_hash = hash_file(K5.__file__)
    columns = {}
    fresh = 0
    for i, source in enumerate(sources[:source_limit]):
        for v in K5.VERT:
            cache = output_dir / "columns" / f"source_{i:02d}_node_{v}_wall{wall}.json"
            if cache.exists():
                img = load_checked(cache, source, v, wall, backend_hash)
            elif fresh < max_new:
                started = time.monotonic()
                img = K5.prune_state(
                    K5.apply_H_cached_state({source: 1.0 + 0.0j}, v, wall), 1e-10)
                save(cache, {
                    "backend_sha256": backend_hash, "input": keyout(source),
                    "input_index": i, "node": v, "wall": wall,
                    "elapsed_seconds": time.monotonic() - started,
                    "rows": [
                        {"key": keyout(k), "amplitude": coeffout(a)}
                        for k, a in sorted(img.items(), key=lambda kv: repr(kv[0]))
                    ],
                })
                fresh += 1
                print(f"source={i}/31 node={v} support={len(img)}", flush=True)
            else:
                return analyse(output_dir, columns, sources, backend_hash, wall)
            columns[i, v] = img
    return analyse(output_dir, columns, sources, backend_hash, wall)

def analyse(output_dir, columns, sources, backend_hash, wall):
    complete = [
        i for i in range(len(sources))
        if all((i, v) in columns for v in K5.VERT)
    ]
    all_keys = sorted(
        set(sources).union(
            *(set(columns[i, v]) for i in complete for v in K5.VERT)),
        key=repr,
    )
    rows = {k: j for j, k in enumerate(all_keys)}
    original = {k: j for j, k in enumerate(sources)}
    C = np.zeros((len(K5.VERT) * len(all_keys), len(complete)), complex)
    Cold = np.zeros((len(K5.VERT) * len(sources), len(complete)), complex)
    per_source = []
    for j, i in enumerate(complete):
        total = changed = inside = 0.0
        for v in K5.VERT:
            for target, amp in columns[i, v].items():
                C[v * len(all_keys) + rows[target], j] = amp
                if target in original:
                    Cold[v * len(sources) + original[target], j] = amp
                    inside += abs(amp)**2
                if 0 in target[0]:
                    changed += abs(amp)**2
                total += abs(amp)**2
        per_source.append({
            "source_index": i, "norm2": total,
            "graph_changed_norm2": changed,
            "inside_original_32_norm2": inside,
        })
    sv = np.linalg.svd(C, compute_uv=False)
    sv_old = np.linalg.svd(Cold, compute_uv=False)
    eig = np.linalg.eigvalsh(C.conj().T @ C)
    old_eig = np.linalg.eigvalsh(Cold.conj().T @ Cold)
    restricted = {
        "backend_sha256": backend_hash, "Jmax2": wall,
        "completed_inputs": complete, "completed_columns": len(complete)*len(K5.VERT),
        "total_original_inputs": 32, "total_required_columns": 160,
        "output_union_dimension": len(all_keys),
        "outside_original_32_dimension": len(all_keys)-32,
        "source_norms": per_source,
        "leakage_safe_singular_values": sv.tolist(),
        "old_compressed_singular_values": sv_old.tolist(),
        "leakage_safe_gram_eigenvalues": eig.tolist(),
        "old_compressed_gram_eigenvalues": old_eig.tolist(),
        "no_uncomputed_column_filled_with_zero": True,
        "input_cross_graph_blocks_available": False,
        "physical_projector_computed": False,
        "claim_scope": ("Gram matrix for complete H_v images of selected "
                        "input sources only; no full BQG projector"),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    save(output_dir / "V15_RESTRICTED_RESULTS.json", restricted)
    if complete:
        np.savez_compressed(output_dir / "V15_RESTRICTED_GRAM.npz",
            safe=C.conj().T@C, old=Cold.conj().T@Cold,
            complete_sources=np.asarray(complete, dtype=int))
    return restricted

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--wall", type=int, default=5)
    p.add_argument("--max-new", type=int, default=160)
    p.add_argument("--source-limit", type=int, default=32)
    a = p.parse_args()
    result = run(a.output_dir, a.wall, a.max_new, a.source_limit)
    print(json.dumps(result, ensure_ascii=False, indent=2))
