#!/usr/bin/env python3
"""One-term V2 PL-S3 Lorentzian worker with tetrahedral charged-volume backend.

This is an execution-only refinement of the preregistered V2 experiment.  The
operator, orientation coefficient, Jmax=7/2 wall, exact-zero semantics and hard
acceptance thresholds are inherited unchanged from the V2/group worker.  Each
process evaluates exactly one of the frozen 24 forward or 24 direct-adjoint
ordered terms so a slow term cannot cancel five completed terms at job timeout.
"""
from __future__ import annotations
import argparse,hashlib,json,subprocess,traceback
from pathlib import Path

import peter_weyl_lorentzian_logical_projection_gate as LP
import pl_lorentzian_group_worker as BASE
import pl_lorentzian_triple_worker as W
from pl_covariant_backend import install_pl_graph
from pl_dual_complex import DualComplex,seed_16cell_boundary
from pl_peter_weyl_euclidean import PLPeterWeylEuclidean
from tetrahedral_volume_backend import install_tetrahedral_volume_backend

OPERATOR_VERSION='tetrahedral-charged-volume-v2'
VOLUME_DEFINITION='V_tet=sqrt(abs((1/4) sum_r (-1)^r q_hat_r)) with production zero-aware spectrum'


def run(mode,index,outdir,source=0):
    try:
        source_commit=subprocess.check_output(
            ['git','rev-parse','HEAD'],text=True,timeout=10).strip()
        dirty=subprocess.check_output(
            ['git','status','--porcelain','--untracked-files=no'],text=True,timeout=10)
    except (OSError,subprocess.SubprocessError) as exc:
        raise ValueError('V2_SOURCE_IDENTITY_NOT_VERIFIABLE') from exc
    if len(source_commit)!=40 or any(ch not in '0123456789abcdef' for ch in source_commit):
        raise ValueError('BAD_V2_SOURCE_COMMIT')
    if dirty.strip():
        raise ValueError('V2_TRACKED_WORKTREE_NOT_CLEAN')
    backend_file=Path(__file__).with_name('tetrahedral_volume_backend.py')
    backend_sha=hashlib.sha256(backend_file.read_bytes()).hexdigest()
    D=DualComplex(seed_16cell_boundary());G=PLPeterWeylEuclidean(D)
    # Validate the frozen orbit index before heavy work.
    W.ordered_spec(D,source,index)
    seed=((1,)*len(G.EDGES),(0,)*D.n_tets)
    oldJ=LP.JMAX2
    with install_tetrahedral_volume_backend():
      with install_pl_graph(G):
        LP.JMAX2=W.JMAX2
        restore,caches=LP.install_sine_cached_stack()
        try:
            state,meta=BASE.evaluate(D,G,seed,source,index,mode,caches)
        finally:
            restore();LP.JMAX2=oldJ
    meta['operator_version']=OPERATOR_VERSION
    meta['volume_definition']=VOLUME_DEFINITION
    meta['source_commit']=source_commit
    meta['volume_backend']='tetrahedral-volume-v1'
    meta['volume_backend_sha256']=backend_sha
    meta['execution_granularity']='one frozen ordered term per process'
    meta['execution_only_note']='No operator, coefficient, cutoff, pruning rule, orientation sign or hard threshold differs from preregistered V2.'
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    W.save_state(outdir/f'term_{mode}_{index}.npz',state,len(G.EDGES),D.n_tets)
    (outdir/f'term_{mode}_{index}.json').write_text(json.dumps(meta,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    return meta


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=('forward','adjoint'),required=True)
    p.add_argument('--index',type=int,choices=range(24),required=True)
    p.add_argument('--source',type=int,default=0)
    p.add_argument('--out-dir',type=Path,required=True)
    a=p.parse_args()
    try:o=run(a.mode,a.index,a.out_dir,a.source);code=0 if o.get('passed') else 1
    except Exception as exc:
        o={'status':'single V2 worker exception','passed':False,'operator_version':OPERATOR_VERSION,'mode':a.mode,'index':a.index,'error_type':type(exc).__name__,'error':str(exc),'traceback':traceback.format_exc()};code=1
        a.out_dir.mkdir(parents=True,exist_ok=True)
        (a.out_dir/f'term_{a.mode}_{a.index}.json').write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps(o,indent=2,sort_keys=True));return code

if __name__=='__main__':raise SystemExit(main())
