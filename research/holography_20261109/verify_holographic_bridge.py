#!/usr/bin/env python3
"""Exact finite holography-inspired bridges for BQG chosen mathematical carriers.
NOT a proof of physical holographic duality or new SU(2) amplitudes.
"""
from itertools import combinations, product
from math import comb, log, sqrt
import json
import numpy as np
import sympy as sp

def run():
    labels=list(product(range(2),repeat=2))
    columns=[(a,b,a^b) for a,b in labels]
    assert len(set(columns))==4
    recover={f"{i}{j}":len({(c[i],c[j]) for c in columns})==4
             for i,j in combinations(range(3),2)}
    assert all(recover.values())

    V=sp.zeros(8,4)
    for col,c in enumerate(columns):
        V[4*c[0]+2*c[1]+c[2],col]=1
    assert V.T*V==sp.eye(4)
    O=sp.Matrix([[1,2,0,0],[2,3,0,0],[0,0,5,7],[0,0,7,11]])
    assert V.T*(V*O*V.T)*V==O
    # This parity map is a classical 2-of-3 recovery, NOT a general
    # quantum code correcting arbitrary erasure of one boundary qubit.
    assert columns[0][0]!=columns[2][0]

    walsh=[sp.Matrix([(-1)**a,(-1)**b,(-1)**(a+b)])/sp.sqrt(3)
           for a,b in labels]
    assert sum(walsh,sp.zeros(3,1))==sp.zeros(3,1)
    for i,j in combinations(range(4),2):
        assert sp.simplify(walsh[i].dot(walsh[j]))==-sp.Rational(1,3)
    assert sum((w*w.T for w in walsh),sp.zeros(3,3))==sp.Rational(4,3)*sp.eye(3)

    nodes=list(product(range(2),repeat=4))
    edges=[(u,v) for u,v in combinations(nodes,2)
           if sum(x!=y for x,y in zip(u,v))==1]
    cut=[(u,v) for u,v in edges if u[0]!=v[0]]
    assert len(nodes)==16 and len(edges)==32 and len(cut)==8
    rank_bound=2**len(cut)
    entropy_bound_bits=len(cut)
    # A product of m Bell pairs is a witness saturating the cut rank
    # and the m-bit entanglement entropy for the chosen abstract network.
    for m in range(1,9):
        d=2**m
        x=np.eye(d)/sqrt(d)
        rho=x@x.T
        assert np.allclose(rho,np.eye(d)/d,atol=1e-15)
        assert abs(float(np.trace(rho@rho))-1/d)<1e-13

    assert comb(8,4)==70
    return {
        "kind":"finite_math_bridge_not_physical_holography",
        "walsh_pairwise_dot":"-1/3",
        "walsh_frame_operator":"4/3 I_3",
        "bulk_classical_states":4,
        "boundary_parity_codewords":["".join(map(str,c)) for c in columns],
        "any_two_classical_flags_recover_label":recover,
        "boundary_encoding_isometry_rank":int(V.rank()),
        "operator_reconstruction_on_global_code_subspace":True,
        "quantum_erasure_correction_for_any_two_qubit_logical_state":False,
        "wall_microstate_count_with_4_of_8_lit":70,
        "wall_missing_information_bits":log(70,2),
        "Q4_vertices":16,"Q4_edges":32,"Q4_coordinate_cut_edges":8,
        "Q4_bond_dimension_assumed":2,
        "Q4_schmidt_rank_bound":rank_bound,
        "Q4_von_neumann_entropy_bound_bits":entropy_bound_bits,
        "Q4_bell_pair_saturation_pure_state":True,
        "original_BQG_S_U_2_amplitudes_computed":False,
        "physical_Bekenstein_Hawking_law_derived_from_BQG":False,
        "AdS_CFT_duality_established_for_BQG":False,
        "physical_quantum_gravity_proved":False,
    }

if __name__=="__main__":
    print(json.dumps(run(),ensure_ascii=False,indent=2))
