#!/usr/bin/env python3
"""Independent production V2 reference for C1, not the staged V4 implementation."""
from __future__ import annotations
import argparse,hashlib,json,re,subprocess,traceback
from pathlib import Path

def verify(commit):
    if not re.fullmatch(r'[a-f0-9]{40}',commit):raise ValueError('INVALID_COMMIT')
    if subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()!=commit:
        raise ValueError('CHECKOUT_SHA_MISMATCH')

def run(mode,index,source,commit,out):
    verify(commit)
    from tetrahedral_volume_backend import install_tetrahedral_volume_backend
    import pl_lorentzian_triple_worker as W
    from pl_dual_complex import DualComplex,seed_16cell_boundary
    backend_sha=hashlib.sha256((Path(__file__).parent/'tetrahedral_volume_backend.py').read_bytes()).hexdigest()
    with install_tetrahedral_volume_backend():state,meta=W.run(index,mode,source)
    if meta.get('passed') is not True or not all(meta.get('checks',{}).values()):
        raise ValueError('GENUINE_V2_NOT_ACCEPTED')
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    d=DualComplex(seed_16cell_boundary())
    W.save_state(out.with_suffix('.npz'),state,len(d.dual_edges()),d.n_tets)
    meta.update(producer_kind='v22_independent_original_v2_worker',source_commit=commit,
      volume_backend_sha256=backend_sha,V4_staged_code_invoked=False,
      npz_sha256=hashlib.sha256(out.with_suffix('.npz').read_bytes()).hexdigest(),HDA_proved=False)
    out.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta

def main():
    p=argparse.ArgumentParser();p.add_argument('--mode',required=True,choices=['forward','adjoint'])
    p.add_argument('--index',type=int,required=True);p.add_argument('--source',type=int,default=0)
    p.add_argument('--commit',required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    try:r=run(a.mode,a.index,a.source,a.commit,a.out);rc=0
    except Exception as e:r=dict(passed=False,status='V2_ERROR_NOT_PHYSICS',error=str(e),traceback=traceback.format_exc());rc=2
    print(json.dumps(r,indent=2));return rc
if __name__=='__main__':raise SystemExit(main())
