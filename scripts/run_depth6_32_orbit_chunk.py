#!/usr/bin/env python3
"""Run one explicit recovery chunk of already-frozen depth-6 [3,2] orbit IDs."""
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


def _write_chunk_manifest(out_dir: Path, *, chunk: int, requested: list[int],
                          rows: list[dict], columns: int, status: str) -> dict:
    payload={
        'schema_version':1,
        'kind':'BQG_DEPTH6_32_KEYMETA_RECOVERY_CHUNK_MANIFEST',
        'status':status,
        'chunk':int(chunk),
        'orbit_ids':sorted(requested),
        'completed_blocks':len(rows),
        'completed_columns':int(columns),
        'records':list(rows),
        'raw_matrices_persisted':False,
        'rank_certified':False,
        'numerical_closure_claimed':False,
        'claim_boundary':(
            'Durable KEYMETA checkpoint for an explicit authorized frozen-orbit recovery chunk. '
            'A partial checkpoint certifies only the listed persisted blocks; it makes no global '
            'actual-q coverage, peeling, rank, kernel, or [3,2] closure claim.'),
    }
    path=Path(out_dir)/f'chunk-{int(chunk):04d}-manifest.json'
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\\n')
    tmp.replace(path)
    return payload


def run_chunk(records, orbit_ids, *, chunk: int, out_dir: Path,
              keymeta_backend, selector_backend, master_backend, basis_backend,
              tol: float=1e-8) -> dict:
    requested=[int(x) for x in orbit_ids]
    if not requested:
        raise RuntimeError('empty orbit chunk')
    if len(requested)!=len(set(requested)):
        raise RuntimeError('duplicate requested orbit id')

    by_id={}
    for rec in records:
        oid=int(rec['orbit_id'])
        if oid in by_id:
            raise RuntimeError(f'duplicate ledger orbit {oid}')
        by_id[oid]=rec
    unknown=set(requested).difference(by_id)
    if unknown:
        raise RuntimeError(f'unknown requested orbit ids: {sorted(unknown)[:10]}')

    out_dir=Path(out_dir)
    block_dir=out_dir/'blocks'
    block_dir.mkdir(parents=True,exist_ok=True)
    rows=[]
    columns=0
    # Persist a manifest checkpoint after every successful block so a runner
    # shutdown does not discard already completed numerical work.
    _write_chunk_manifest(out_dir,chunk=chunk,requested=requested,rows=rows,
                          columns=columns,status='PASS_KEYMETA_RECOVERY_CHUNK_PARTIAL')
    for oid in sorted(requested):
        rec=by_id[oid]
        m=int(rec['m'])
        keymeta_path=block_dir/f'orbit-{oid:04d}.keymeta.pkl.gz'
        meta=keymeta_backend.compute_and_persist_block(
            rec,selector_backend,master_backend,basis_backend,
            keymeta_path=keymeta_path,raw_path=None,tol=tol)
        if int(meta['orbit_id'])!=oid or int(meta['m'])!=m:
            raise RuntimeError(f'orbit {oid}: persisted KEYMETA identity mismatch')
        rows.append({
            'orbit_id':oid,
            'm':m,
            'keymeta_file':keymeta_path.name,
            'keymeta_sha256':_sha256_file(keymeta_path),
            'q_key_count':int(meta['q_key_count']),
            'rows':int(meta['rows']),
        })
        columns+=m
        _write_chunk_manifest(out_dir,chunk=chunk,requested=requested,rows=rows,
                              columns=columns,status='PASS_KEYMETA_RECOVERY_CHUNK_PARTIAL')

    return _write_chunk_manifest(
        out_dir,chunk=chunk,requested=requested,rows=rows,columns=columns,
        status='PASS_KEYMETA_RECOVERY_CHUNK_COMPLETE')

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--ledger',type=Path,required=True)
    ap.add_argument('--gate',type=Path,required=True)
    ap.add_argument('--orbit-ids-json',required=True,
                    help='JSON array of already-frozen orbit IDs')
    ap.add_argument('--chunk',type=int,required=True)
    ap.add_argument('--out-dir',type=Path,required=True)
    ap.add_argument('--tol',type=float,default=1e-8)
    a=ap.parse_args()

    ledger=json.loads(a.ledger.read_text())
    gate=json.loads(a.gate.read_text())
    import depth6_32_replay_authorization as A
    auth=A.authorize_replay(ledger,gate)
    if not auth['numerical_replay_authorized']:
        raise RuntimeError('numerical replay not authorized')
    orbit_ids=json.loads(a.orbit_ids_json)
    if not isinstance(orbit_ids,list):
        raise RuntimeError('orbit-ids-json must decode to a list')

    import bqg_depth6_generic_jucys_selector as selector_backend
    import bqg_depth6_generic_master_row as master_backend
    import depth6_32_frozen_m_basis as basis_backend
    import depth6_32_keymeta_first as keymeta_backend

    result=run_chunk(
        ledger['records'],orbit_ids,chunk=a.chunk,out_dir=a.out_dir,
        keymeta_backend=keymeta_backend,selector_backend=selector_backend,
        master_backend=master_backend,basis_backend=basis_backend,tol=a.tol)
    result['authorization_kind']=auth['kind']
    result['identity_sha256']=auth['identity_sha256']
    print(json.dumps(result,sort_keys=True))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
