#!/usr/bin/env python3
"""v1.1: one REAL KKV auxiliary trace path, source all-jhalf, exact sine HE.
A path is NOT a complete ordered triple. Eight path outputs need coherent sum.
Do not equate timeout with zero; do not equate this raw result with physical HL.
"""
import argparse, itertools, json, math, sys, time, hashlib, base64, zlib
from pathlib import Path
HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_lorentzian_sine_ordered_triple_gate as SINE
import peter_weyl_lorentzian_ordered_triple_gate as RAW
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

def norm(state):
    return math.sqrt(sum(abs(z)**2 for z in state.values()))

def oriented_24(neighbors):
    result=[]
    for face,omit in enumerate(neighbors):
        base=tuple(x for x in neighbors if x!=omit)
        for perm in itertools.permutations(base):
            ids=[base.index(x) for x in perm]
            inversions=sum(ids[i]>ids[j] for i in range(3) for j in range(i+1,3))
            result.append(((-1)**(face+inversions),perm))
    if len(result)!=24: raise AssertionError('24 terms required')
    return result

def compute(source_v=0,index=18,indices=(0,0,0),Jmax2=7,stop_after='CK1',emit_checkpoint=False):
    ZVM.patch_and_clear()
    if not (0 <= index <24):raise ValueError('term index outside 0..23')
    if Jmax2 < 7: raise ValueError('Jmax2 must be >=7 for all-jhalf')
    if len(indices)!=3 or any(x not in (0,1) for x in indices):
        raise ValueError('auxiliary trace indices must be 0 or 1')
    i,j,k=indices
    sign,(a,b,c)=oriented_24(PW.NEIG[source_v])[index]
    gauss_key=PW.basis_full_jhalf()[0]
    cov=CV.gauss_to_covariant({gauss_key:1+0j},source_v)
    record={'source_v':source_v,'term':index,'orientation':sign,
       'ordered_neighbors':[a,b,c],'auxiliary_path':[i,j,k],
       'stop_after':stop_after,'Jmax2':Jmax2,'steps':[],
       'physical_HL_computed':False,'complete_aux_trace':False}
    def stage(name, out, dt, diag=None):
        res={'stage':name,'support':len(out),'norm':norm(out),'elapsed_s':dt}
        if diag is not None:res['diagnostics']=diag
        record['steps'].append(res)
        print('REAL_PATH_STEP',json.dumps(res),flush=True)
    def execute():
        t=time.perf_counter()
        q,leak=RAW.COMP.C_volume_component(cov,source_v,c,k,i,Jmax2)
        stage('CV',q,time.perf_counter()-t,{'basis_leakage':float(leak)})
        if not q or stop_after=='CV':return
        t=time.perf_counter()
        out,diag=RAW.KCOMP.C_K_component(q,source_v,b,j,k,Jmax2)
        stage('CK1',out,time.perf_counter()-t,diag)
        # Preserve exact complex coefficients (NOT just their norms) and measure
        # final-scalar-eligible J=0,1 before the expensive second C(K).
        byJ={}
        for key,amp in out.items():
            J2=int(key[2])
            byJ.setdefault(str(J2),[0,0.0])
            byJ[str(J2)][0]+=1
            byJ[str(J2)][1]+=abs(amp)**2
        scalar={key:amp for key,amp in out.items() if key[2] in (0,2)}
        stat={'counts_norm2_by_J2':byJ,'J01_support':len(scalar),
              'J01_norm':norm(scalar),'total_support':len(out),'total_norm':norm(out)}
        record['CK1_J_STAT']=stat
        print('CK1_J_STAT',json.dumps(stat,sort_keys=True),flush=True)
        if emit_checkpoint:
            payload={'schema':'BQG_CK1_EXACT_V1','source_v':source_v,'term':index,
              'indices':[i,j,k],'Jmax2':Jmax2,'ordered_neighbors':[a,b,c],
              'operator':'RAW.KCOMP.C_K_component after RAW.COMP.C_volume_component in SINE.run',
              'state':[[repr(key),float(z.real),float(z.imag)]
                       for key,z in sorted(out.items(),key=lambda kv:repr(kv[0]))]}
            raw=json.dumps(payload,separators=(',',':'),sort_keys=True).encode('utf-8')
            compressed=base64.b64encode(zlib.compress(raw,9)).decode('ascii')
            record['CK1_checkpoint_sha256']=hashlib.sha256(raw).hexdigest()
            print('CK1_CHECKPOINT_BASE64',compressed,flush=True)
            print('CK1_CHECKPOINT_SHA256',record['CK1_checkpoint_sha256'],flush=True)
        if not out or stop_after=='CK1':return
        out={key:amp for key,amp in out.items() if key[2] in (0,2)}
        stage('SCALAR_ELIGIBLE_J01',out,0)
        if not out:return
        t=time.perf_counter()
        out,diag=RAW.KCOMP.C_K_component(out,source_v,a,i,j,Jmax2)
        stage('CK2',out,time.perf_counter()-t,diag)
        record['sparse_states']=[[repr(key),float(z.real),float(z.imag)]
                                 for key,z in sorted(out.items(),key=lambda kv:repr(kv[0]))]
    SINE.run(source_v,executor=execute)
    record['completed_probed_path']=True
    return record

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=int,default=0);p.add_argument('--term',type=int,default=18)
    p.add_argument('--i',type=int,default=0);p.add_argument('--j',type=int,default=0);p.add_argument('--k',type=int,default=0)
    p.add_argument('--stop-after',choices=['CV','CK1','CK2'],default='CK1')
    p.add_argument('--Jmax2',type=int,default=7);p.add_argument('--output',type=Path)
    p.add_argument('--emit-checkpoint',action='store_true')
    a=p.parse_args()
    data=compute(a.source,a.term,(a.i,a.j,a.k),a.Jmax2,a.stop_after,a.emit_checkpoint)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print('REAL_PATH_COMPLETE',json.dumps({k:v for k,v in data.items() if k!='sparse_states'}),flush=True)
if __name__=='__main__':main()
