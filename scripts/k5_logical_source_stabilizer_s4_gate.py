#!/usr/bin/env python3
"""Exact source-node stabilizer S4 representation on the 32D K5 logical carrier.

A global vertex permutation fixing source node 0 induces (i) a permutation of
the four spectator logical factors and (ii) a local four-leg recoupling matrix
at every node.  This gate builds the resulting 32x32 unitaries explicitly from
the already tested one-cell logical S4 representation, checks the group law,
and decomposes the 32D representation into S4 irreps.

The result does not assume Lorentzian covariance.  It only supplies the exact
domain representation needed for a future held-out test that H_L^dag H_L lies
in its commutant before symmetry is used to reduce expensive columns.
"""
from __future__ import annotations
import argparse, itertools, json, math, sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW
import logical_s4_twirl_gate as LS4


def compose(g,h):
    # Apply h first, then g.
    return tuple(g[h[v]] for v in range(5))


def invperm(g):
    out=[0]*5
    for i,j in enumerate(g): out[j]=i
    return tuple(out)


def local_slot_permutation(g,v):
    """New slot at g(v) takes the old slot whose neighbour maps there."""
    gv=g[v]; gi=invperm(g); q=[]
    for new_neighbor in PW.NEIG[gv]:
        old_neighbor=gi[new_neighbor]
        q.append(PW.NEIG[v].index(old_neighbor))
    return tuple(q)


def global_U(g, local_basis):
    local={v:LS4.logical_representation(local_slot_permutation(g,v),local_basis) for v in range(5)}
    bits=list(itertools.product(range(2),repeat=5)); idx={b:i for i,b in enumerate(bits)}
    U=np.zeros((32,32),complex)
    for kin in bits:
        col=idx[kin]
        for kout in bits:
            amp=1+0j
            for v in range(5):
                amp*=local[v][kout[g[v]],kin[v]]
                if abs(amp)<1e-15: break
            if abs(amp)>1e-15: U[idx[kout],col]=amp
    return U


def cycle_type_on_neighbors(g):
    seen=set(); lengths=[]
    for x in (1,2,3,4):
        if x in seen: continue
        y=x;n=0
        while y not in seen:
            seen.add(y);n+=1;y=g[y]
        lengths.append(n)
    return tuple(sorted(lengths,reverse=True))

CLASSES={(1,1,1,1):("1",1),(2,1,1):("2",6),(2,2):("22",3),(3,1):("3",8),(4,):("4",6)}
IRREPS={
    "[4]":{"1":1,"2":1,"22":1,"3":1,"4":1},
    "[31]":{"1":3,"2":1,"22":-1,"3":0,"4":-1},
    "[22]":{"1":2,"2":0,"22":2,"3":-1,"4":0},
    "[211]":{"1":3,"2":-1,"22":-1,"3":0,"4":1},
    "[1111]":{"1":1,"2":-1,"22":1,"3":1,"4":-1},
}


def run():
    lb=LS4.singlet_basis()
    perms=[]
    for p in itertools.permutations((1,2,3,4)):
        g=(0,)+tuple(p); perms.append(g)
    reps={g:global_U(g,lb) for g in perms}
    unit=max(float(np.linalg.norm(U.conj().T@U-np.eye(32))) for U in reps.values())

    # Check both composition conventions; one must realize the group exactly.
    err_direct=0.;err_reverse=0.
    for g in perms:
        for h in perms:
            err_direct=max(err_direct,float(np.linalg.norm(reps[g]@reps[h]-reps[compose(g,h)])))
            err_reverse=max(err_reverse,float(np.linalg.norm(reps[g]@reps[h]-reps[compose(h,g)])))
    convention="homomorphism" if err_direct<=err_reverse else "anti_homomorphism"
    group_err=min(err_direct,err_reverse)

    chars={name:[] for name,_ in [v for v in CLASSES.values()]}
    class_values={}
    for ctype,(name,size) in CLASSES.items():
        vals=[np.trace(reps[g]) for g in perms if cycle_type_on_neighbors(g)==ctype]
        class_values[name]={"size":size,"count":len(vals),"character_real":float(np.mean(vals).real),"character_imag":float(np.mean(vals).imag),"spread":float(max(abs(z-np.mean(vals)) for z in vals))}

    mult={}
    for ir,chi in IRREPS.items():
        s=0j
        for name,row in class_values.items(): s+=row["size"]*np.conj(chi[name])*complex(row["character_real"],row["character_imag"])
        mult[ir]=int(round((s/24).real))
    dim_check=sum(mult[k]*IRREPS[k]["1"] for k in mult)
    comm_dim=sum(x*x for x in mult.values())
    char_comm=(sum(row["size"]*(row["character_real"]**2+row["character_imag"]**2) for row in class_values.values())/24)

    checks={
        "24_stabilizer_elements":len(perms)==24,
        "global_representation_unitary":unit<2e-11,
        "group_law":group_err<2e-10,
        "characters_constant_on_classes":max(row["spread"] for row in class_values.values())<2e-10,
        "irrep_multiplicities_nonnegative":all(x>=0 for x in mult.values()),
        "irrep_dimensions_sum_to_32":dim_check==32,
        "commutant_dimension_character_identity":abs(comm_dim-char_comm)<2e-8,
    }
    return {
        "status":"exact K5 source-node stabilizer S4 representation on 32D logical carrier",
        "passed":bool(all(checks.values())),"source_node":0,"logical_dimension":32,"group_order":24,
        "max_unitarity_error":unit,"group_convention":convention,"group_law_error":group_err,
        "class_characters":class_values,"irrep_multiplicities":mult,"dimension_check":dim_check,
        "commutant_complex_dimension":comm_dim,"commutant_dimension_from_character":float(char_comm),
        "checks":checks,
        "scientific_use":"Before any Lorentzian column reduction, directly verify that the measured source-node Gram M_L=H_L^dag H_L commutes with every U_g. If that held-out covariance passes, M_L is restricted to this commutant; symmetry may then reduce the number of independent microscopic columns/parameters.",
        "claim_boundary":"Representation theory only; this gate does not assume or prove H_L covariance, does not reconstruct M_L, and does not alter the preregistered Lorentzian operator."
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path);a=ap.parse_args();out=run();txt=json.dumps(out,indent=2);print(txt)
    if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(txt+'\n',encoding='utf-8')
    return 0 if out['passed'] else 1

if __name__=='__main__':raise SystemExit(main())
