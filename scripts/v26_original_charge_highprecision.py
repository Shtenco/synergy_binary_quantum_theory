#!/usr/bin/env python3
"""Re-evaluate original one-orientation SU(2) HE^sine with independent raw norm.

Original SU(2) branch generator is NEVER replaced. Only raw_gram_norm2
postprocessing is swapped to direct mpmath coordinate aggregation.
This is not all-12 orientation coherence and not a Lorentzian HDA proof.
"""
from __future__ import annotations
import argparse,hashlib,json,re,subprocess,sys,traceback
from pathlib import Path
from v26_direct_branch_norm import squared_norm

def sha256(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def attest(here):
    original=here/"bqg_he_sine_true_coherent_charge_orientation.py"
    if not original.is_file():raise RuntimeError("ORIGINAL_SU2_PRODUCER_UNAVAILABLE")
    commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=here,text=True,timeout=12).strip()
    if re.fullmatch("[a-f0-9]{40}",commit) is None:
        raise RuntimeError("UNVERIFIABLE_COMMIT")
    # Refuse to claim the HEAD identifies the actual source when tracked code differs.
    dirty=subprocess.check_output(["git","status","--porcelain","--untracked-files=no"],
                                  cwd=here,text=True,timeout=12).strip()
    if dirty:raise RuntimeError("TRACKED_SOURCE_DIRTY")
    return {"source_commit":commit,"original_producer_sha256":sha256(original),
            "precision_kernel_sha256":sha256(here/"v26_direct_branch_norm.py")}

def execute(a):
    here=Path(__file__).resolve().parent
    meta=attest(here)
    if a.expected_commit and meta["source_commit"]!=a.expected_commit:
        raise RuntimeError("SOURCE_COMMIT_MISMATCH")
    import bqg_he_sine_true_coherent_charge_orientation as original
    replay=[]
    original_raw=original.raw_gram_norm2
    def replacement(branches,block=64):
        data=squared_norm(branches,original.PW.VERT,digits=a.digits,
                          max_group_dim=a.max_dim,max_work=a.max_work)
        replay.append(data)
        return (float(data["squared_norm"]),data["number_of_spin_sectors"],
                max((sector["branches"] for sector in data["sectors"]),default=0))
    original.raw_gram_norm2=replacement
    try:
        native=original.run(a.source,a.target,a.orientation,a.Jmax2,a.basis_rank)
    finally:
        original.raw_gram_norm2=original_raw
    if len(replay)!=1:raise RuntimeError("ORIGINAL_RAW_NORM_NOT_INVOKED_EXACTLY_ONCE")
    return {"status":"ORIGINAL_SU2_ORIENTATION_RAW_NORM_REPLAYED",
            "source_attestation":meta,
            "original_orientation_result":native,
            "new_high_precision_raw_norm":replay[0],
            "binary64_SU2_element_errors_unbounded":True,
            "all_12_orientation_coherent_sum_available":False,
            "charge_covariance_proved":False,
            "new_lorentzian_V4_paths":0,
            "full_HDA_proved":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--orientation",type=int,choices=range(12),required=True)
    p.add_argument("--source",type=int,default=1)
    p.add_argument("--target",type=int,default=0)
    p.add_argument("--Jmax2",type=int,default=3)
    p.add_argument("--basis-rank",type=int,default=0)
    p.add_argument("--digits",type=int,default=80)
    p.add_argument("--max-dim",type=int,default=4096)
    p.add_argument("--max-work",type=int,default=2000000)
    p.add_argument("--expected-commit")
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    try:
        report=execute(a);rc=0
    except Exception as exc:
        report={"status":"INCOMPLETE_NOT_PHYSICAL",
                "passed":False,"error":str(exc),"error_type":type(exc).__name__,
                "charge_covariance_proved":False,"full_HDA_proved":False}
        rc=2
    a.output.write_text(json.dumps(report,ensure_ascii=False,allow_nan=False,indent=2)+"\n")
    print(json.dumps(report,ensure_ascii=False,allow_nan=False,indent=2))
    return rc

if __name__=="__main__":sys.exit(main())
