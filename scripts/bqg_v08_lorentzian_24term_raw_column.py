#!/usr/bin/env python3
"""Actual 24-term raw Lorentzian K-K-V column; not yet Hermitian HL or physical kernel.
All operator amplitudes come from existing Peter-Weyl SINE + RAW code.
Outputs are coherent complex sums, never norms of individual terms.
"""
import argparse,itertools,json,math
from pathlib import Path
def sign(base,p):
    ix=[base.index(q) for q in p]
    return -1 if sum(ix[i]>ix[j] for i in range(3) for j in range(i+1,3))%2 else 1
def ordered(neigh):
    out=[]
    for r,omit in enumerate(neigh):
        tri=tuple(x for x in neigh if x!=omit)
        for p in itertools.permutations(tri):
            out.append(((-1)**r*sign(tri,p),p))
    assert len(out)==24
    return out
def add(target,source,c=1):
    for k,v in source.items():target[k]=target.get(k,0j)+c*v
def nrm(s):return math.sqrt(sum(abs(v)**2 for v in s.values()))
def raw(v,start,count,Jmax2,index):
    import k5_peter_weyl_safe_hda_column as PW
    import peter_weyl_lorentzian_sine_ordered_triple_gate as SINE
    import peter_weyl_lorentzian_ordered_triple_gate as RAW
    import peter_weyl_covariant_volume_leg_gate as CV
    import peter_weyl_zeroaware_volume_migration_experiment as ZVM
    ZVM.patch_and_clear()
    ini=PW.basis_full_jhalf()[index]
    psi=CV.gauss_to_covariant({ini:1+0j},v)
    terms=ordered(PW.NEIG[v])
    total={};row=[]
    def executor():
        for k in range(start,min(start+count,24)):
            eps,(a,b,c)=terms[k]; tmp={}
            max_leak=0.
            stage={'after_CV':0,'after_CK1':0,'after_J02':0,'after_CK2':0,'CV_nonzero_aux_paths':0,'CK1_nonzero_aux_paths':0,'CK2_nonzero_aux_paths':0}
            for i,j,l in itertools.product(range(2),repeat=3):
                q,leak=RAW.COMP.C_volume_component(psi,v,c,l,i,Jmax2)
                max_leak=max(max_leak,float(leak))
                if q:
                    stage['after_CV']+=len(q); stage['CV_nonzero_aux_paths']+=1
                if not q:continue
                q,d=RAW.KCOMP.C_K_component(q,v,b,j,l,Jmax2)
                max_leak=max(max_leak,float(d['outer_complete_basis_leakage']),
                    float(d['internal_volume_sector_leakage']))
                if q:
                    stage['after_CK1']+=len(q);stage['CK1_nonzero_aux_paths']+=1
                q={key:z for key,z in q.items() if key[2] in (0,2)}
                stage['after_J02']+=len(q)
                if not q:continue
                q,d=RAW.KCOMP.C_K_component(q,v,a,i,j,Jmax2)
                max_leak=max(max_leak,float(d['outer_complete_basis_leakage']),
                    float(d['internal_volume_sector_leakage']))
                if q:
                    stage['after_CK2']+=len(q);stage['CK2_nonzero_aux_paths']+=1
                add(tmp,q)
            add(total,tmp,eps)
            row.append({'index':k,'sign':eps,'edges':[a,b,c],'support':len(tmp),
                'norm':nrm(tmp),'max_basis_leak':max_leak,'stage_counts':stage})
        return None
    SINE.run(v,executor=executor)
    return dict(vertex=v,input_index=index,Jmax2=Jmax2,start=start,count=len(row),
      terms=row,norm_raw=nrm(total),scalar_norm=nrm({k:z for k,z in total.items() if k[2]==0}),
      non_scalar_norm=nrm({k:z for k,z in total.items() if k[2]!=0}),
      support=len(total),
      states=[[repr(k),z.real,z.imag] for k,z in sorted(total.items(),key=lambda p:repr(p[0])) if abs(z)>1e-13],
      scope='raw 24-term amplitude only; NO Hermitian completion; NO full HL or full physical master')
def main():
    p=argparse.ArgumentParser()
    p.add_argument('--v',type=int,default=0);p.add_argument('--input',type=int,default=0)
    p.add_argument('--Jmax2',type=int,default=7);p.add_argument('--start',type=int,default=0)
    p.add_argument('--count',type=int,default=1);p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.Jmax2<7: p.error('all-jhalf input needs Jmax2 >=7')
    if not(0<=a.start<24 and 1<=a.count<=24-a.start):p.error('bad shard range')
    result=raw(a.v,a.start,a.count,a.Jmax2,a.input)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(result,indent=2)+'\n')
    # JSON records full sparse states to local file; GitHub log stays bounded.
    print(json.dumps({k:v for k,v in result.items() if k!='states'},indent=2))
if __name__=='__main__':main()
