#!/usr/bin/env python3
"""Finite PL16 E-only output-habitat dimension/conditioning theorem (V26).

No Lorentzian S/R, independent D, physical projector or full off-shell HDA.
The theorem assumes P is orthogonal and the archived ranks identify exact
finite maps with a common output Hilbert inner product.
"""
import json, math
from pathlib import Path
BASE=Path(__file__).resolve().parent

def replay():
    d=json.loads((BASE/'v12_results.json').read_text())
    e=json.loads((BASE/'v13_real_results.json').read_text())
    w=d['full_vs_compressed']; c=d['C3'];w43=d['Krylov_W43']; c13=e['C3']
    n=w['W17_dimension']; rP=w['compressed_E0_rank'];rE=w['full_E0_map_rank']
    nr=n-rP
    assert (n,rP,rE,nr)==(17,2,17,15)
    assert c['lost_commutator_fraction_on_W17']==1.0
    assert c13['measured_commutator_rank']==10
    assert w43['total_dimension']==43 and w43['seed_Cr_output_retention_fraction']==1.0
    s=w['full_E0_map_min_singular'];sC=c13['minimum_nonzero_commutator_singular']
    return {
        'scope':'published PL16 E-only finite maps; interpreted as exact matrix ranks',
        'W17_dim':n,'E0_compressed_rank':rP,'E0_full_rank':rE,
        'W17_compressed_kernel_dim':nr,
        'outside_W17_E0_on_compressed_kernel_rank_at_least':nr,
        'minimal_habitat_dim_for_W17_and_all_E0_on_compressed_kernel':n+nr,
        'outside_E0_on_kernel_Frobenius_lower':math.sqrt(nr)*s,
        'outside_scaled_E0_on_kernel_Frobenius_lower':(2/3)*math.sqrt(nr)*s,
        'EE_ten_commutator_rank_outside_W17':c13['measured_commutator_rank'],
        'minimal_habitat_dim_for_W17_and_ten_commutators':n+c13['measured_commutator_rank'],
        'minimal_combined_one_hit_habitat_dim_at_least':n+max(nr,c13['measured_commutator_rank']),
        'W43_existing_dimension':w43['total_dimension'],
        'W43_excess_above_minimal_dimension':w43['total_dimension']-(n+max(nr,c13['measured_commutator_rank'])),
        'EE_commutator_rank_stable_if_operator_error_less_than':sC,
        'EE_commutator_rank_stable_if_each_column_error_less_than':sC/math.sqrt(c13['measured_commutator_rank']),
        'not_proven':'W43 invariance under full H, 32 V4, D, physical kernel or continuum',
        'full_HDA_proved':False,
        'fresh_SU2_amplitudes':0,
    }
if __name__=='__main__':print(json.dumps(replay(),indent=2))
