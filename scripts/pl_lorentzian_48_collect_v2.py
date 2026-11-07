#!/usr/bin/env python3
"""Provenance-locked collector for tetrahedral-volume PL Lorentzian V2 terms."""
from __future__ import annotations
import argparse,hashlib,json,re,subprocess,traceback
from pathlib import Path
import pl_lorentzian_48_collect as BASE
VERSION='tetrahedral-charged-volume-v2'

def preflight(root):
    rows=[]
    for p in Path(root).rglob('term_*.json'):
        d=json.loads(p.read_text(encoding='utf-8'));rows.append((p,d))
    if len(rows)!=48:raise RuntimeError(f'need exactly 48 V2 term metadata files, got {len(rows)}')
    bad=[str(p) for p,d in rows if d.get('operator_version')!=VERSION]
    if bad:raise RuntimeError(f'non-V2/superseded term metadata present: {bad[:8]}')
    commits={d.get('source_commit') for p,d in rows}
    backend_ids={d.get('volume_backend') for p,d in rows}
    backend_shas={d.get('volume_backend_sha256') for p,d in rows}
    if len(commits)!=1 or not all(isinstance(x,str) and re.fullmatch(r'[0-9a-f]{40}',x) for x in commits):
        raise RuntimeError('V2_TERM_SOURCE_COMMITS_MISSING_OR_MIXED')
    if backend_ids!={'tetrahedral-volume-v1'} or len(backend_shas)!=1:
        raise RuntimeError('V2_TERM_VOLUME_BACKEND_IDENTITY_MISSING_OR_MIXED')
    try:
        actual_head=subprocess.check_output(['git','rev-parse','HEAD'],text=True,timeout=10).strip()
        dirty=subprocess.check_output(['git','status','--porcelain','--untracked-files=no'],text=True,timeout=10)
    except (OSError,subprocess.SubprocessError) as exc:
        raise RuntimeError('V2_COLLECTOR_SOURCE_IDENTITY_NOT_VERIFIABLE') from exc
    if actual_head not in commits:
        raise RuntimeError('V2_COLLECTOR_SOURCE_COMMIT_MISMATCH')
    if dirty.strip():
        raise RuntimeError('V2_COLLECTOR_TRACKED_WORKTREE_NOT_CLEAN')
    backend_file=Path(__file__).with_name('tetrahedral_volume_backend.py')
    actual_backend_sha=hashlib.sha256(backend_file.read_bytes()).hexdigest()
    if backend_shas!={actual_backend_sha}:
        raise RuntimeError('V2_COLLECTOR_BACKEND_SHA_MISMATCH')
    return rows

def run(root):
    preflight(root);L,Ld,S,out=BASE.run(root)
    out['operator_version']=VERSION
    out['charged_volume_definition']='V_tet=sqrt(abs((1/4) sum_r (-1)^r q_hat_r)) with production zero-aware spectrum'
    out['superseded_input_rejected']=True
    out['science_status']='AMPLITUDE_PRECURSOR_S_NODE0_V2'
    return L,Ld,S,out

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--json-output',type=Path,required=True);p.add_argument('--state-output',type=Path,required=True);a=p.parse_args()
    try:L,Ld,S,out=run(a.root);code=0 if out['passed'] else 1
    except Exception as exc:L=Ld=S={};out={'status':'V2 collector exception','passed':False,'operator_version':VERSION,'error_type':type(exc).__name__,'error':str(exc),'traceback':traceback.format_exc()};code=1
    a.json_output.parent.mkdir(parents=True,exist_ok=True);a.json_output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8');BASE.save_bundle(a.state_output,L,Ld,S);print(json.dumps(out,indent=2,sort_keys=True));return code
if __name__=='__main__':raise SystemExit(main())
