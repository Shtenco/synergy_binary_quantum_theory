#!/usr/bin/env python3
"""BQG V24 finite *compressed* chiral-index check; not a physical Lorentzian result."""
import json
import numpy as np
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def chiral_index(n_even,n_odd,rank_b):
    if any(type(x) is not int or x<0 for x in (n_even,n_odd,rank_b)):
        raise ValueError('invalid dimensions')
    if rank_b>min(n_even,n_odd):raise ValueError('invalid block rank')
    nullity=n_even+n_odd-2*rank_b
    return {'even':n_even,'odd':n_odd,'rank_B':rank_b,
            'nullity':nullity,'protected_index_min':abs(n_even-n_odd),
            'additional_unprotected_zeros':nullity-abs(n_even-n_odd)}

def verify():
    v12=json.loads((ROOT/'v12_results.json').read_text())
    v13=json.loads((ROOT/'v13_real_results.json').read_text())
    w17=v12['full_vs_compressed']
    w43=v13['W43']
    i17=chiral_index(1,16,1)
    i43=chiral_index(w43['even_dim'],w43['odd_dim'],w43['rank_compressed_E0']//2)
    assert i17['nullity']==w17['compressed_E0_nullity']==15
    assert i43['nullity']==w43['nullity_compressed_E0']==11
    assert i17['additional_unprotected_zeros']==i43['additional_unprotected_zeros']==0
    assert w17['full_E0_map_rank']==17 and w17['full_E0_nullity']==0
    p=np.array(w43['positive_eigenvalues_ascending'])
    n=np.array(w43['negative_eigenvalues_ascending'])
    assert len(p)==len(n)==16
    pairing=float(np.max(np.abs(p+n[::-1])))
    assert pairing<1e-10
    return {'W17':i17,'W43':i43,'W43_numerical_plusminus_pairing':pairing,
      'scope':'compressed finite E-only chiral operators only; NOT the full H_L',
      'mathematical_identity':'nullity([[0,B],[B*,0]]) = n_even + n_odd - 2 rank(B)',
      'full_Lorentzian_physical_kernel_proved':False}
if __name__=='__main__':
    out=verify()
    print(json.dumps(out,indent=2))
