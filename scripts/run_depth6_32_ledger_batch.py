#!/usr/bin/env python3
"""Run one deterministic depth-6 [3,2] KEYMETA-first shard from frozen records only."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def _sha256_file(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def run_batch(records, *, shard: int, nshards: int, out_dir: Path,
              keymeta_backend, selector_backend, master_backend, basis_backend,
              tol: float=1e-8) -> dict:
    shard=int(shard); nshards=int(nshards)
    if nshards<=0 or not 0<=shard<nshards:
        raise RuntimeError('invalid shard/nshards')
    seen=set()
    chosen=[]
    for rec in records:
        oid=int(rec['orbit_id'])
        if oid in seen:
            raise RuntimeError(f'duplicate orbit {oid}')
        seen.add(oid)
        if int(rec.get('source_shard',-1))==shard:
            chosen.append(rec)
    if not chosen:
        raise RuntimeError(f'no records assigned to shard {shard}')

    out_dir=Path(out_dir)
    block_dir=out_dir/'blocks'
    block_dir.mkdir(parents=True,exist_ok=True)
    manifest=[]
    columns=0
    for rec in sorted(chosen,key=lambda r:int(r['orbit_id'])):
        oid=int(rec['orbit_id'])
        m=int(rec['m'])
        keymeta_path=block_dir/f'orbit-{oid:04d}.keymeta.pkl.gz'
        meta=keymeta_backend.compute_and_persist_block(
            rec,selector_backend,master_backend,basis_backend,
            keymeta_path=keymeta_path,raw_path=None,tol=tol,
        )
        if int(meta['orbit_id'])!=oid or int(meta['m'])!=m:
            raise RuntimeError(f'orbit {oid}: persisted KEYMETA identity mismatch')
        manifest.append({
            'orbit_id':oid,
            'm':m,
            'keymeta_file':keymeta_path.name,
            'keymeta_sha256':_sha256_file(keymeta_path),
            'q_key_count':int(meta['q_key_count']),
            'rows':int(meta['rows']),
        })
        columns+=m

    payload={
        'schema_version':1,
        'kind':'BQG_DEPTH6_32_KEYMETA_FIRST_SHARD_MANIFEST',
        'status':'PASS_KEYMETA_FIRST_SHARD_COMPLETE',
        'shard':shard,
        'nshards':nshards,
        'completed_blocks':len(manifest),
        'completed_columns':columns,
        'records':manifest,
        'raw_matrices_persisted':False,
        'rank_certified':False,
        'numerical_closure_claimed':False,
        'claim_boundary':(
            'Complete KEYMETA persistence for one authorized frozen shard only. '
            'No global actual-q coverage, peeling, rank, kernel, or [3,2] closure claim.'
        ),
    }
    (out_dir/f'shard-{shard:03d}-manifest.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    return payload


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--ledger',type=Path,required=True)
    ap.add_argument('--gate',type=Path,required=True)
    ap.add_argument('--shard',type=int,required=True)
    ap.add_argument('--nshards',type=int,default=112)
    ap.add_argument('--out-dir',type=Path,required=True)
    ap.add_argument('--tol',type=float,default=1e-8)
    a=ap.parse_args()

    ledger=json.loads(a.ledger.read_text())
    gate=json.loads(a.gate.read_text())
    import depth6_32_replay_authorization as A
    auth=A.authorize_replay(ledger,gate)
    if not auth['numerical_replay_authorized']:
        raise RuntimeError('numerical replay not authorized')
    if int(a.nshards)!=int(auth['nshards']):
        raise RuntimeError('nshards mismatch with authorization')

    import bqg_depth6_generic_jucys_selector as selector_backend
    import bqg_depth6_generic_master_row as master_backend
    import depth6_32_frozen_m_basis as basis_backend
    import depth6_32_keymeta_first as keymeta_backend

    result=run_batch(
        ledger['records'],shard=a.shard,nshards=a.nshards,out_dir=a.out_dir,
        keymeta_backend=keymeta_backend,selector_backend=selector_backend,
        master_backend=master_backend,basis_backend=basis_backend,tol=a.tol,
    )
    result['authorization_kind']=auth['kind']
    result['identity_sha256']=auth['identity_sha256']
    print(json.dumps(result,sort_keys=True))
    return 0


if __name__=='__main__': raise SystemExit(main())
