#!/usr/bin/env python3
"""BQG V25: finite PL16/W17 E-only geometry, NOT full Lorentzian HDA."""
import json, math
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def run():
    v12=json.loads((ROOT/"v12_results.json").read_text())
    v13=json.loads((ROOT/"v13_real_results.json").read_text())
    w=v12["full_vs_compressed"]
    c=v12["C3"]
    c13=v13["C3"]
    assert (w["compressed_E0_rank"],w["compressed_E0_nullity"],w["full_E0_map_rank"])==(2,15,17)
    assert c["lost_commutator_fraction_on_W17"]==1.0
    assert c13["measured_commutator_rank"]==10
    sig=w["full_E0_map_min_singular"]
    norm=c["norm_total_C3"]
    result={
      "PL16_W17": {
        "nullity_PEP":15,
        "rank_full_EP":17,
        "offhabitat_rank_at_least":15,
        "offhabitat_Frobenius_lower":math.sqrt(15)*sig,
        "offhabitat_scaled_Frobenius_lower":(2/3)*math.sqrt(15)*sig,
        "offhabitat_scaled_operator_lower":(2/3)*sig,
        "theorem":"For x in ker(PEP) with Px=x, E P x=(I-P)E P x."
      },
      "PL16_HDA_ten_probes":{
        "EE_commutator_Frobenius_norm":norm,
        "EE_commutator_orthogonal_to_W17":True,
        "if_all_nonEE_and_independent_D_compensation_inside_W17_residual_F_lower":(4/9)*norm,
        "if_compensation_inside_W17_residual_operator_lower":(4/9)*c13["maximum_commutator_singular"],
        "theorem":"If PC=0 and PK=K, then ||(4/9)C+K||_F^2=(4/9)^2||C||_F^2+||K||_F^2."
      },
      "full_Lorentzian_HDA_proved":False,
      "fresh_physical_SU2_amplitudes":0
    }
    return result

if __name__=="__main__": print(json.dumps(run(),indent=2,ensure_ascii=False))
