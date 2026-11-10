#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,tarfile
from pathlib import Path
import numpy as np
try:
    from .rank_matrix_guard import validate_full_column_matrix, validate_shard_coverage, validate_index_coverage, validate_index_value, validate_positive_integer, parse_indexed_map_filename
    from .safe_archive import safe_extract_tar
    from .safe_pickle import safe_load_numpy_pickle
except ImportError:
    from rank_matrix_guard import validate_full_column_matrix, validate_shard_coverage, validate_index_coverage, validate_index_value, validate_positive_integer, parse_indexed_map_filename, validate_index_value, parse_indexed_map_filename
    from safe_archive import safe_extract_tar
    from safe_pickle import safe_load_numpy_pickle
import itertools
import bqg_depth6_41 as Z
import bqg_safe_boundary_master_gate as R

def gauss_ok(sp):
    return all(R.allowed_k2_t(*R.local_spins(sp,v)) for v in R.VERT)

def geometric_outs(sp):
    sp=tuple(sp);J=max(sp)+1;out=set()
    for a,b in itertools.combinations(R.NEIG[0],2):
        es=[R.EIDX[tuple(sorted((0,a)))],R.EIDX[tuple(sorted((a,b)))],R.EIDX[tuple(sorted((b,0)))]]
        for ds in itertools.product((-1,1),repeat=3):
            z=list(sp);ok=True
            for e,dd in zip(es,ds):
                z[e]+=dd
                if z[e]<0 or z[e]>J:ok=False;break
            if ok:
                z=tuple(z)
                if gauss_ok(z):out.add(Z.canon4(z)[0])
    return out

def extract_artifacts(root,work,expected):
    root=Path(root);work=Path(work);work.mkdir(parents=True,exist_ok=True)
    ts=sorted(root.rglob('bqg-s4sign-shard-*.tar.gz'))
    if len(ts)!=expected:raise SystemExit(f'expected {expected} tarballs, got {len(ts)}')
    ds=[]
    for t in ts:
        d=work/t.stem.replace('.tar','');d.mkdir(exist_ok=True)
        with tarfile.open(t, 'r:gz') as f:safe_extract_tar(f, d)
        ds.append(d)
    return ds

def load_index(ds,expected):
    summaries=[];records={};paths={}
    for d in ds:
        ss=list(d.rglob('summary.json'))
        if len(ss)!=1:raise SystemExit(f'bad summary count {d} {len(ss)}')
        s=json.load(open(ss[0]));summaries.append(s)
        if (s['shell_states'],s['s4sign_blocks'],s['s4sign_dim'])!=(264962,11956,130112):raise SystemExit(('ledger',s))
        if s['failed_blocks']!=0 or s['ok_blocks']!=s['assigned_blocks']:raise SystemExit(('incomplete shard',s))
        for r in s['records']:
            i=validate_index_value(r.get('i'), "S4-sign record index")
            validate_positive_integer(r.get('d'), f"S4-sign block {i} dimension")
            if i in records:raise SystemExit(('duplicate record',i))
            records[i]=r
        for p in d.rglob('b*.pkl'):
            i=parse_indexed_map_filename(p.stem, "S4-sign map filename")
            if i in paths:raise SystemExit(('duplicate map',i))
            paths[i]=p
    validate_shard_coverage(
        [s.get('shard') for s in summaries], expected, "S4-sign aggregate"
    )
    validate_index_coverage(
        records.keys(), paths.keys(), 11956, "S4-sign aggregate"
    )
    if sum(int(r['d']) for r in records.values())!=130112:raise SystemExit('column total mismatch')
    return records,paths

def certify(records,paths,outfile):
    remaining=set(records);support={};occ={}
    for i,r in records.items():
        ss=geometric_outs(tuple(r['qin']));support[i]=ss
        for q in ss:occ.setdefault(q,set()).add(i)
    rounds=[];cert={}
    while True:
        ready=[]
        for i in sorted(remaining):
            uq=[q for q in support[i] if occ[q]=={i}]
            if not uq:continue
            with open(paths[i],'rb') as f:block=safe_load_numpy_pickle(f)
            mats=[block[q] for q in uq if q in block]
            del block
            if not mats:continue
            M=np.vstack(mats);d=int(records[i]['d']);M=validate_full_column_matrix(M,d,f"block {i}")
            if M.shape[0]<d:continue
            s=np.linalg.svd(M,compute_uv=False);smax=float(s[0]) if len(s) else 0.
            tol=max(M.shape)*np.finfo(float).eps*max(smax,1.)*100
            rank=int((s>tol).sum());sigma=float(s[-1]) if len(s)>=d else 0.
            if rank==d and sigma>max(tol,1e-10):ready.append((i,d,len(uq),M.shape[0],sigma,tol))
        if not ready:break
        rounds.append({'round':len(rounds),'blocks':len(ready),'columns':sum(x[1] for x in ready)})
        for i,d,nq,nrows,sigma,tol in ready:
            cert[str(i)]={'d':d,'unique_q':nq,'rows':nrows,'sigma_min':sigma,'rank_tol':tol}
        for i,*_ in ready:remaining.remove(i)
        for i,*_ in ready:
            for q in support[i]:occ[q].discard(i)
        print('ROUND',rounds[-1],'remaining',len(remaining),flush=True)
    res={'status':'PASS' if not remaining else 'RESIDUAL_CORE','scope':'depth6 S4-sign input = [1^5] + [2,1,1,1]','target_blocks':11956,'target_columns':130112,'certified_blocks':len(cert),'certified_columns':sum(v['d'] for v in cert.values()),'remaining_blocks':len(remaining),'remaining_columns':sum(int(records[i]['d']) for i in remaining),'rounds':rounds,'global_sigma_min':min((v['sigma_min'] for v in cert.values()),default=None),'remaining_ids':sorted(remaining),'proof_logic':'Conservative geometric-support occupancy; each removed S4-sign block has full-column-rank numeric H0 image on rows unique at its removal step.'}
    Path(outfile).write_text(json.dumps(res,indent=2)+'\n')
    print(json.dumps({k:res[k] for k in ['status','certified_blocks','certified_columns','remaining_blocks','remaining_columns','global_sigma_min']},indent=2),flush=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--artifacts',required=True);ap.add_argument('--expected-shards',type=int,required=True);ap.add_argument('--work',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
    ds=extract_artifacts(a.artifacts,a.work,a.expected_shards);rec,paths=load_index(ds,a.expected_shards);certify(rec,paths,a.out)
if __name__=='__main__':main()
