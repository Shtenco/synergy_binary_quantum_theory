#!/usr/bin/env python3
"""V24. Spectral lower bound for a NAIVE unitary matching of full local
S4-twirled Euclidean multiplicity blocks. NOT full BQG physical RG.
Published source: BQG_SERIAL_MULTIPLICITY_MASTER_RESULT.md
Git blob: 08f4de387cb7115b87121dd029de0917c56ed7b6
"""
import json
import numpy as np
from scipy.optimize import linprog

SCALES=[
 ("2",[5.067680982954,5.807713077901]),
 ("5/2",[5.769201915254,22.499018693032]),
 ("3",[9.001223358412,10.264684071156]),
 ("7/2",[7.140843485222,9.459669912782,36.381268968404]),
 ("4",[7.819342874141,14.040428696350,16.140265919741]),
]
def solve(a,b):
    """min_{s>=0,U unitary} || B U-s U A ||_2 for Hermitian equal-m blocks.
    Equality with min_s max_i |b_i-s*a_i| follows from Weyl inequalities."""
    a=np.asarray(a,float);b=np.asarray(b,float)
    if a.shape!=b.shape:raise ValueError("No square unitary across different m")
    if (a<=0).any() or (b<=0).any():raise ValueError("positive spectra expected")
    constraints=np.array([[float(x),-1] for x in a]+[[-float(x),-1] for x in a])
    opt=linprog([0,1],A_ub=constraints,b_ub=np.r_[b,-b],
                bounds=[(0,None),(0,None)],method="highs")
    if not opt.success:raise RuntimeError(opt.message)
    s=float(opt.x[0]);e=float(max(abs(s*a-b)))
    if abs(e-opt.x[1])>1e-8:raise RuntimeError("LP not attained")
    return s,e

def compute():
    output=[]
    for i in range(len(SCALES)-1):
        j,a=SCALES[i];next_j,b=SCALES[i+1]
        if len(a)!=len(b):
            output.append({"from":j,"to":next_j,"status":"UNRESOLVED_MULTIPLICITY_CHANGE",
                           "need":"Real microscopic ancilla-aware isometry; not a square unitary"})
            continue
        s,e=solve(a,b)
        gamma_a=float(a[1]-a[0]);gamma_b=float(b[1]-b[0])
        output.append({"from":j,"to":next_j,"m":len(a),
                       "best_common_positive_scale":s,
                       "unavoidable_full_block_spectral_error":e,
                       "unscaled_error":float(np.max(abs(np.array(a)-np.array(b)))),
                       "gamma_selection_from":gamma_a,"gamma_selection_to":gamma_b,
                       "is_physical_master_Rj":False})
    assert abs(output[0]["unavoidable_full_block_spectral_error"]-7.403132181051562)<1e-9
    assert abs(output[1]["unavoidable_full_block_spectral_error"]-5.069284677800814)<1e-9
    assert abs(output[3]["unavoidable_full_block_spectral_error"]-7.812384162910653)<1e-9
    return {"source":"Published local Euclidean S4-twirled multiplicity spectra",
            "original_evidence_blob":"08f4de387cb7115b87121dd029de0917c56ed7b6",
            "results":output,
            "full_BQG_RG_certified":False,
            "nonzero_full_physical_projector_certified":False,
            "physical_TT_and_EH_certified":False}
if __name__=="__main__":print(json.dumps(compute(),indent=2,ensure_ascii=False))
