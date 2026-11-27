#!/usr/bin/env python3
"""V19 canonical K5 direct-column replay with strict source/cache provenance.

Run on the repository's real scripts/k5_peter_weyl_safe_hda_column.py;
mirror data are comparison fixtures, NEVER silently promoted to canonical.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, time
from pathlib import Path

EXPECTED_BLOB = "e7e1e0705c4b833dcee167a2036ef5f8929a4388"
MIRROR_SHA256 = "9babf9e985f066c9f76346f1f58e2683a7ab2bb4d212056d143e90fd1e8fd081"

def git_blob(data: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()

def load_k5(path: Path):
    raw = path.read_bytes()
    actual = git_blob(raw)
    if actual != EXPECTED_BLOB:
        raise RuntimeError(f"Refusing noncanonical K5: {actual} != {EXPECTED_BLOB}")
    spec = importlib.util.spec_from_file_location("bqg_v19_upstream_k5", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, hashlib.sha256(raw).hexdigest()

def unpack(rows):
    out = {}
    for row in rows:
        key = (tuple(row["key"][0]), tuple(row["key"][1]))
        coeff = row.get("coefficient", row.get("amplitude"))
        if coeff is None or key in out:
            raise ValueError("Malformed/duplicate column row")
        out[key] = complex(float.fromhex(coeff[0]), float.fromhex(coeff[1]))
    return out

def matching_cache(base_dirs, source_index, node):
    name = f"source_{source_index:02d}_vertex_{node}_wall5.json"
    matches = [Path(root)/name for root in base_dirs if (Path(root)/name).is_file()]
    if len(matches) != 1:
        raise RuntimeError(f"Need exactly one mirror fixture for {source_index=}, {node=}; got {len(matches)}")
    fixture = json.loads(matches[0].read_text())
    if fixture["backend_sha256"] != MIRROR_SHA256:
        raise RuntimeError("Fixture is not the declared V14 research mirror")
    if fixture["node"] != node or fixture["source_index"] != source_index or fixture["wall"] != 5:
        raise RuntimeError("Fixture identity mismatch")
    return fixture

def replay(k5, source_index, node, dirs, outdir):
    input_key = k5.basis_full_jhalf()[source_index]
    fixture = matching_cache(dirs, source_index, node)
    if fixture["input"] != [list(input_key[0]), list(input_key[1])]:
        raise RuntimeError("Fixture K5 input key differs")
    started = time.monotonic()
    computed = k5.prune_state(
        k5.apply_H_cached_state({input_key: 1.+0j}, node, 5), 1e-10)
    old = unpack(fixture["rows"])
    keys = computed.keys() | old.keys()
    diff = max((abs(computed.get(key, 0j)-old.get(key, 0j)) for key in keys), default=0.0)
    same = computed.keys() == old.keys()
    record = {
        "source_index": source_index, "node": node,
        "canonical_support": len(computed), "mirror_support": len(old),
        "same_output_keys": same, "max_abs_difference": float(diff),
        "seconds": time.monotonic()-started,
    }
    outdir.mkdir(parents=True, exist_ok=True)
    if not same or diff > 1e-12:
        raise AssertionError("Canonical/research discrepancy: " + json.dumps(record))
    serial = [{
        "key": [list(key[0]), list(key[1])],
        "coefficient": [float(amp.real).hex(), float(amp.imag).hex()],
    } for key, amp in sorted(computed.items(),key=lambda pair: repr(pair[0]))]
    payload = {
        "canonical_git_blob": EXPECTED_BLOB,
        "research_mirror_fixture_sha256": MIRROR_SHA256,
        "record": record, "complete_canonical_cold_replay": True,
        "rows": serial,
    }
    result_path = outdir / f"canonical_s{source_index:02d}_v{node}.json"
    result_path.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    return record

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--backend",type=Path,default=Path("scripts/k5_peter_weyl_safe_hda_column.py"))
    p.add_argument("--mirror-column-dirs",type=Path,nargs="+",required=True)
    p.add_argument("--pairs",nargs="+",default=["0:0","0:1","7:2","15:3","23:4","31:0"])
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    k5,sha=load_k5(a.backend)
    all_records=[]
    for token in a.pairs:
        i,v=map(int,token.split(":"))
        if not (0<=i<32 and 0<=v<5):
            raise ValueError("Only initial 32 K5 sources, nodes 0..4")
        record=replay(k5,i,v,a.mirror_column_dirs,a.output_dir)
        all_records.append(record)
        print(json.dumps(record),flush=True)
    report={
        "canonical_git_blob_verified":True,
        "canonical_sha256":sha,
        "canonical_replayed_initial_columns":len(all_records),
        "full_160_canonical_coverage":len(set((r["source_index"],r["node"]) for r in all_records)) == 160,
        "complete_first_shell_incoming":False,
        "physical_projector_computed":False,
        "results":all_records,
    }
    (a.output_dir/"V19_REPORT.json").write_text(json.dumps(report,indent=2)+"\n")
    print("PASS; no uncomputed column filled with zero")
if __name__=="__main__":
    main()
