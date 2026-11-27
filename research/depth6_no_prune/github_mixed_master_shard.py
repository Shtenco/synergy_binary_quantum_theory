#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,math,os,pickle,time,multiprocessing as mp
from pathlib import Path
import sys
import numpy as np

SOURCE_DIR = Path(__file__).resolve().parents[1] / 'depth6_recovered_source'
sys.path.insert(1, str(SOURCE_DIR))
import bqg_radial_return_reachability_gate as Reach
import bqg_depth6_41 as Z
import bqg_depth6_generic_jucys_selector as J
import bqg_depth6_generic_master_row as M
import bqg_safe_boundary_master_gate as R
import bqg_covariant_volume_recoupling as Rec

TARGET={
 '32':(130903,2755),
 '311':(153455,2719),
 '221':(130503,2749),
 '2111':(103318,2712),
}

def shell(depth=6):
    cur={(1,)*10}
    for k in range(depth):
        nxt=set()
        for v in Reach.V:nxt |= Reach.apply_H(cur,v)
        _,cur=Reach.summarize(nxt)
        print('SHELL',k+1,len(cur),flush=True)
    return tuple(sorted(cur))

def cost_balanced(active,nshards):
    # active entries: (orbit_id, rep, m, coord_dim, proxy)
    bins=[[] for _ in range(nshards)];load=[0.0]*nshards
    for z in sorted(active,key=lambda x:(x[4],x[2],x[0]),reverse=True):
        j=min(range(nshards),key=lambda b:(load[b],b));bins[j].append(z);load[j]+=z[4]
    return bins,load

_RECS=None;_KEY=None;_OUT=None

def init_worker(recs,key,outdir):
    global _RECS,_KEY,_OUT
    _RECS=recs;_KEY=key;_OUT=Path(outdir)

def work(pos):
    oid,rep,m,d,proxy=_RECS[pos];t=time.time()
    try:
        W=J.fast_selector(tuple(rep),_KEY)
        if W.shape!=(d,m):
            # d here is coordinate dimension
            if W.shape[1]!=m:raise RuntimeError(('selector shape',oid,W.shape,d,m))
        mp0,_=M.master_map(tuple(rep),_KEY,W)
        p=_OUT/f'b{oid:05d}.pkl';pickle.dump(mp0,open(p,'wb'),protocol=5)
        return {'orbit_index':oid,'rep':list(rep),'m':m,'coord_dim':int(W.shape[0]),'row_blocks':len(mp0),
                'rows':int(sum(A.shape[0] for A in mp0.values())),'sec':time.time()-t,'ok':True}
    except Exception as e:
        return {'orbit_index':oid,'rep':list(rep),'m':m,'coord_dim':d,'sec':time.time()-t,'ok':False,'error':repr(e)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--irrep',choices=TARGET,required=True);ap.add_argument('--shard',type=int,required=True);ap.add_argument('--shards',type=int,default=128);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--out',type=Path,default=Path('artifact'))
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True);t0=time.time()
    sp=shell(6)
    if len(sp)!=264962:raise SystemExit(f'bad shell {len(sp)}')
    reps=sorted({Z.canon5(tuple(s))[0] for s in sp})
    if len(reps)!=2757:raise SystemExit(f'bad S5 orbit count {len(reps)}')
    seed=J.CFG[a.irrep][0];active=[];total=0
    for oid,rep in enumerate(reps):
        m=J.exact_mult(rep,a.irrep)
        if not m:continue
        d=J.layout(rep,seed)[2]
        # selector/H cost grows with coordinate rows and multiplicity; add floor so tiny blocks still count.
        proxy=float(max(1,d)*max(1,m))
        active.append((oid,rep,m,d,proxy));total+=m
    target,ac=TARGET[a.irrep]
    if total!=target or len(active)!=ac:raise SystemExit(f'bad ledger {a.irrep} {len(active)} {total}')
    bins,loads=cost_balanced(active,a.shards);mine=bins[a.shard]
    print('LEDGER',a.irrep,len(active),total,'SHARD',a.shard,'N',len(mine),'COLS',sum(z[2] for z in mine),'PROXY',loads[a.shard],flush=True)
    # Prewarm local volume blocks for input S5 orbit members. Runtime output blocks extend the cache as needed.
    lsset=set()
    for oid,rep,m,d,proxy in mine:
        for p in Z.S5:
            q=Z.S.permute_spins(rep,p) if hasattr(Z,'S') else rep
            for v in range(5):lsset.add(tuple(R.local_spins(q,v)))
    tw=time.time();nv=0
    for ls in lsset:
        for JJ in Rec.allowed_total_J(ls):Rec.reduced_volume_block(ls,JJ);nv+=1
    print('PREWARM',len(lsset),nv,'sec',time.time()-tw,flush=True)
    ctx=mp.get_context('fork');records=[]
    with ctx.Pool(a.workers,initializer=init_worker,initargs=(mine,a.irrep,a.out)) as pool:
        for n,z in enumerate(pool.imap_unordered(work,range(len(mine)),chunksize=1),1):
            records.append(z)
            if n%5==0 or n==len(mine):print('DONE',n,'/',len(mine),'ok',sum(x['ok'] for x in records),'sec',round(time.time()-t0,1),flush=True)
    records.sort(key=lambda x:x['orbit_index'])
    summary={'irrep':a.irrep,'shard':a.shard,'shards':a.shards,'shell_states':len(sp),'s5_orbits':len(reps),'active_blocks':len(active),'target_columns':total,'assigned_blocks':len(mine),'assigned_columns':sum(z[2] for z in mine),'ok_blocks':sum(x['ok'] for x in records),'failed_blocks':sum(not x['ok'] for x in records),'elapsed_sec':time.time()-t0,'records':records}
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    if summary['failed_blocks']:raise SystemExit('failed blocks')
if __name__=='__main__':main()
