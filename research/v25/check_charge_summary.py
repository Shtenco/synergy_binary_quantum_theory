#!/usr/bin/env python3
"""BQG V25 K5 charge diagnostic. This is a summary-only replay, not a full projector."""
import json, math
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def run():
    d=json.loads((ROOT/"charge_12_extracted.json").read_text())
    rows=d["rows"]
    assert len(rows)==12 and {r["index"] for r in rows}==set(range(12))
    assert d["full_all_orientation_complex_sum_available"] is False
    impossible=[r["index"] for r in rows if r["raw2"]==0 and r["projected2"]>0]
    unstable=[r["index"] for r in rows if r["numerical_stability"]=="CANCELLATION_DOMINATED_NOT_CERTIFIED"]
    r6=rows[6]
    return {
      "source_blob_sha":d["source_git_blob_sha"],
      "orientation_count":12,
      "sum_separate_raw2":math.fsum(r["raw2"] for r in rows),
      "sum_separate_primitive_incoherent2":math.fsum(r["incoherent2"] for r in rows),
      "inexact_zero_projection_violations":impossible,
      "unstable_cancellation_orientations":unstable,
      "orientation6_incoherent_to_raw_norm_ratio":math.sqrt(r6["incoherent2"]/r6["raw2"]),
      "full_coherent_all_orientation_charged_leak_unknown":True,
      "not_a_physical_HDA_result":True
    }

if __name__=="__main__":print(json.dumps(run(),indent=2,ensure_ascii=False))
