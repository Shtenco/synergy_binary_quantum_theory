#!/usr/bin/env python3
"""Independent C(V) consistency audit before Lorentzian ordered triples.

Compare the Gauss-input full covariant C(V) composition engine with the independent
reference finite-leg generator on the SAME gauge input and same cutoff. No mixing
with H_L and no inference of physical kernel.
"""
import json,argparse
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_covariant_composition_gate as C
import peter_weyl_zeroaware_volume_migration_experiment as ZVM
def norm2(s):return sum(abs(v)**2 for v in s.values())
def nrm(s):return norm2(s)**0.5
def difference(a,b):
    return sum(abs(a.get(k,0)-b.get(k,0))**2 for k in set(a)|set(b))**0.5
def run(v=0,index=0):
    ZVM.patch_and_clear()
    input_key=PW.basis_full_jhalf()[index]
    psi=CV.gauss_to_covariant({input_key:1+0j},v)
    records=[]
    for target in PW.NEIG[v]:
        for cutoff in (3,7):
            for i in range(2):
                for j in range(2):
                    got,leak=C.C_volume_component(psi,v,target,i,j,cutoff)
                    row={"v":v,"target":target,"cutoff":cutoff,"i":i,"j":j,
                         "norm":nrm(got),"support":len(got),"leak":float(leak)}
                    if cutoff==3:
                        # Reference on original Gauss input, same source/target.
                        direct=dict(C.KG.local_volume_column(input_key,v))
                        if i==j:
                            base=CV.gauss_to_covariant(direct,v)
                        else:base={}
                        # independent 2-hit Covariant V from CV module.
                        hVh,_=CV.inverse_then_forward(input_key,v,target,i,j,cutoff,True)
                        ref=dict(base)
                        C.add(ref,hVh,-1)
                        row["independent_reference_norm"]=nrm(ref)
                        row["independent_difference"]=difference(got,ref)
                    records.append(row)
    # A true acceptance gate: compare independent volume constructions and both cutoffs.
    errs=[x["independent_difference"] for x in records if "independent_difference" in x]
    assert max(errs)<1e-10, ("C(V) independent implementations disagree",max(errs))
    for target in PW.NEIG[v]:
        for i in range(2):
            for j in range(2):
                a=next(x for x in records if x["target"]==target and x["cutoff"]==3 and x["i"]==i and x["j"]==j)
                b=next(x for x in records if x["target"]==target and x["cutoff"]==7 and x["i"]==i and x["j"]==j)
                assert abs(a["norm"]-b["norm"])<1e-11,("C(V) cutoff mismatch",target,i,j)
    if v==0 and index==0:
        assert all(x["norm"]>0.1 for x in records if x["target"] in (1,2,3))
        assert all(x["norm"]<1e-13 for x in records if x["target"]==4)
    print(json.dumps({"passed":True,"max_independent_difference":max(errs),
      "input":"all_jhalf","node":v,"input_index":index,
       "records":records,"note":"This is a comparison of real C(V) implementations, not an HL amplitude."},indent=2))
if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--v",type=int,default=0);ap.add_argument("--input",type=int,default=0)
    a=ap.parse_args();run(a.v,a.input)
