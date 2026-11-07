#!/usr/bin/env python3
"""Build SHA-256 pins for V21 stage files; store resulting manifest independently."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def attest_manifest(manifest):
    if not isinstance(manifest, dict) or not isinstance(manifest.get("paths"), list):
        raise ValueError("BAD_V21_MANIFEST")
    for index, item in enumerate(manifest["paths"]):
        if not isinstance(item, dict):
            raise ValueError("BAD_PATH_ENTRY")
        pins = {}
        for stage in ("stage1", "stage2", "stage3"):
            if stage not in item:
                raise ValueError("MISSING_" + stage.upper())
            prefix = Path(item[stage])
            meta_path = prefix.with_suffix(".json")
            npz_path = prefix.with_suffix(".npz")
            metadata = json.loads(meta_path.read_text("utf-8"))
            if metadata.get("stage") != int(stage[-1]):
                raise ValueError("STAGE_NUMBER_MISMATCH_" + stage.upper())
            if metadata.get("passed") is not True:
                raise ValueError("STAGE_NOT_PASSED_" + stage.upper())
            npz_digest = sha256(npz_path)
            if metadata.get("npz_sha256") != npz_digest:
                raise ValueError("NPZ_DIGEST_MISMATCH_" + stage.upper())
            if metadata.get("output_state_sha256") != npz_digest:
                raise ValueError("OUTPUT_DIGEST_MISMATCH_" + stage.upper())
            pins[stage] = {
                "metadata_sha256": sha256(meta_path),
                "npz_sha256": npz_digest,
            }
        item["stage_attestations"] = pins
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True,
                        help="JSON manifest with paths and stage1/stage2/stage3 prefixes")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text("utf-8"))
    result = attest_manifest(manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({"status": "ATTESTATIONS_CREATED",
                      "paths": len(result["paths"]),
                      "warning": "Store output manifest independently; a co-located mutable manifest is not an authenticity anchor."},
                     indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
