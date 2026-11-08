#!/usr/bin/env python3
"""V24 minimal independent exact verifier. Scope: printed finite K5/W5 matrix,
plus separately archived PL16 E-only spectra. No new Lorentzian amplitudes.
"""
import json,math,hashlib
from fractions import Fraction as F
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
HASH="950423420850c887f8facb747be66a44cc4ed007655c6b04aac8ce21b347dbc9"
BLOB="0231c54d2e9a7fc9401995c999d3989f3f1f7648"
BRACKETS=[("32.41356","32.41357"),("32.81846","32.81848"),
          ("33.62288","33.62290"),("34.12068","34.12070"),
          ("36.90496","36.90498")]

def source():
    p=ROOT/"SOURCE_SNAPSHOT.json"
    if hashlib.sha256(p.read_bytes()).hexdigest()!=HASH:
        raise ValueError("SOURCE_HASH_MISMATCH")
    d=json.loads(p.read_text(),parse_float=str)
    if d["source_git_blob_sha"]!=BLOB:raise ValueError("SOURCE_BLOB_MISMATCH")
    return d

def complex_fraction(x):return F(x[0]),F(x[1])
def mul(x,y):return x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def conj(x):return x[0],-x[1]
def add(x,y):return x[0]+y[0],x[1]+y[1]
def sub(x,y):return x[0]-y[0],x[1]-y[1]

def inertia(a,shift):
    """Exact rational complex Hermitian LDL* inertia at non-eigenvalue shift."""
    n=len(a);L=[[(F(0),F(0)) for _ in range(n)] for _ in range(n)];D=[]
    for j in range(n):
        q=a[j][j][0]-F(shift)
        for k in range(j):q-=D[k]*(L[j][k][0]**2+L[j][k][1]**2)
        if q==0:raise ValueError("ZERO_PIVOT")
        D.append(q);L[j][j]=(F(1),F(0))
        for i in range(j+1,n):
            v=a[i][j]
            for k in range(j):v=sub(v,mul(mul(L[i][k],(D[k],F(0))),conj(L[j][k])))
            L[i][j]=(v[0]/q,v[1]/q)
    return sum(x<0 for x in D),D

def decimal_floor_sqrt(x,places=12):
    n=10**places
    q=math.isqrt((x.numerator*n*n)//x.denominator)
    assert F(q*q,n*n)<=x<F((q+1)**2,n*n)
    return f"{q//n}.{q%n:0{places}d}"

def verify():
    src=source()
    G=[[complex_fraction(z) for z in row] for row in src["complex_gram_re_im"]]
    assert len(G)==5 and all(len(row)==5 for row in G)
    assert all(G[i][j]==conj(G[j][i]) for i in range(5) for j in range(5))
    neg,pivots=inertia(G,"0")
    assert neg==0 and all(p>0 for p in pivots)
    intervals=[]
    for i,(lo,hi) in enumerate(BRACKETS):
        assert inertia(G,lo)[0]==i and inertia(G,hi)[0]==i+1
        intervals.append({"i":i+1,"lower_strict":lo,"upper_strict":hi})
    bounds=[]
    for k in range(1,6):
        lo=[F(z[0]) for z in BRACKETS]
        bounds.append({"k":k,"Q_operator_norm_strictly_above":
                       decimal_floor_sqrt(F(4,9)*lo[k-1]),
                       "Q_F_norm_strictly_above":
                       decimal_floor_sqrt(F(4,9)*sum(lo[:k]))})
    v12=json.loads((ROOT/"v12_results.json").read_text())
    v13=json.loads((ROOT/"v13_real_results.json").read_text())
    w17=v12["full_vs_compressed"]
    w43=v13["W43"]
    assert (w17["compressed_E0_rank"],w17["full_E0_map_rank"])==(2,17)
    assert (w43["even_dim"],w43["odd_dim"],w43["rank_compressed_E0"])==(27,16,32)
    pairerr=max(abs(x+y) for x,y in zip(w43["positive_eigenvalues_ascending"],
                                     reversed(w43["negative_eigenvalues_ascending"])))
    C=v13["C3"]
    assert C["measured_commutator_rank"]==10
    smin=C["minimum_nonzero_commutator_singular"]
    Fnorm=float(src["commutator_F_norm"])
    smax=C["maximum_commutator_singular"]
    hda=[]
    for r in range(10):
        floor=(4/9)*math.sqrt(max((10-r)*smin**2,Fnorm**2-r*smax**2,0))
        hda.append({"rank_at_most":r,"conditional_F_residual_at_least":floor})
    result={
      "scope":"W5 exact printed Gram, W17/W43 and HDA archived PL16 floats",
      "W5_eigen_brackets":intervals,"W5_kernel_nullity_Q_necessary_bounds":bounds,
      "W17_compressed_zero_count":15,"W17_full_E0_kernel_dim":0,
      "W17_corrected_single_node_gap":(2/3)*w17["full_E0_map_min_singular"],
      "W43_chiral_index_min_zero":11,"W43_chiral_zero_count":w43["nullity_compressed_E0"],
      "W43_nonzero_plus_minus_eigen_pairing_error":pairerr,
      "HDA_ten_probe_rank":10,"HDA_conditional_F_residual_bounds":hda,
      "HDA_six_distance_two_nonEE_required_F":(4/9)*C["one_hop_difference_conditional_gate"]["euclidean_commutator_combined_norm"],
      "HDA_six_distance_two_nonEE_required_rank":6,
      "new_true_SU2_V4_amplitudes":0,"physical_HDA_proved":False,
      "physical_continuum_proved":False,"interacting_TT_proved":False}
    return result

if __name__=="__main__":
    d=verify()
    (ROOT/"V24_REPLAY.json").write_text(json.dumps(d,indent=2)+"\n")
    print(json.dumps(d,indent=2))
