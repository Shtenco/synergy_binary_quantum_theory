#!/usr/bin/env python3
"""First real Peter-Weyl Lorentzian ordered triple amplitude gate.

This is the smallest physical-amplitude object beyond the C(V)/C(K)
prerequisites.  On the frozen all-j=1/2 Gauss input at source node v, evaluate

    T_abc = Tr_aux[ C_a(K) C_b(K) C_c(V) ]

with genuine state-to-state Peter-Weyl operators and right-to-left action

    C_c(V) -> C_b(K) -> C_a(K).

The auxiliary trace is

    sum_{i,j,k=0,1} C_a(K)_{ij} C_b(K)_{jk} C_c(V)_{ki}.

No precomputed columns are multiplied.  The source total-J label is retained
between legs.  Jmax=7/2 is the preregistered exact single-H_L support wall.

This gate deliberately computes ONE ordered triple only.  It is therefore not
yet the epsilon-oriented Lorentzian node Hamiltonian.  Its purpose is to test
whether the real three-leg amplitude is nonzero, closes back to an SU(2) scalar
after the auxiliary trace, respects the cutoff and has no complete-basis charge
leakage.

Overall constants, beta and the later Hermitian Lorentzian ordering are absent;
this is the raw structural core amplitude.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0,str(HERE))

import k5_peter_weyl_safe_hda_column as PW
import bqg_pruning_audit as PRUNE
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_covariant_composition_gate as COMP
import peter_weyl_covariant_K_composition_gate as KCOMP

JMAX2=7
TOL=1e-11


def add(dst,src,scale=1.0,tol=TOL):
    for key,amp in src.items():
        prior=dst.get(key,0j)
        incoming=scale*amp
        z=prior+incoming
        if abs(z)>tol:
            dst[key]=z
        else:
            PRUNE.record_merge("lorentzian_ordered_triple.add",key,prior,incoming,z,tol)
            if key in dst:
                del dst[key]


def norm2(state):
    return float(sum(abs(a)**2 for a in state.values()))


def weight_by_source_J(state):
    out={}
    for key,amp in state.items():
        J2=key[2]
        out[J2]=out.get(J2,0.0)+abs(amp)**2
    return {str(J2/2):float(x) for J2,x in sorted(out.items())}


def max_spin(state):
    return max((max(key[0]) for key in state),default=0)/2


def update_diag(dst,src):
    for name,val in src.items():
        if isinstance(val,(int,float)):
            dst[name]=max(dst.get(name,0.0),float(val))


def run(source_v=0,a=None,b=None,c=None):
    # Apply the preregistered zero-aware exact-nullspace convention before any
    # physical amplitude is generated.  This changes no frozen H_E threshold.
    import peter_weyl_zeroaware_volume_migration_experiment as ZVM
    ZVM.patch_and_clear()

    neigh=PW.NEIG[source_v]
    if a is None or b is None or c is None:
        a,b,c=neigh[:3]
    if len({a,b,c})!=3 or any(x not in neigh for x in (a,b,c)):
        raise ValueError('a,b,c must be three distinct neighbors of source_v')

    initial=PW.basis_full_jhalf()[0]
    psi=CV.gauss_to_covariant({initial:1+0j},source_v)

    total={}
    aux_path_rows=[]
    diagmax={
        'CV_complete_basis_leakage':0.0,
        'CK_outer_complete_basis_leakage':0.0,
        'CK_internal_volume_sector_leakage':0.0,
        'CK_complete_charge_basis_leakage':0.0,
    }

    # (ABC)_ii=sum_jk A_ij B_jk C_ki, action read right-to-left.
    for i,j,k in itertools.product(range(2),repeat=3):
        s1,leakV=COMP.C_volume_component(
            psi,source_v,c,k,i,JMAX2
        )
        diagmax['CV_complete_basis_leakage']=max(
            diagmax['CV_complete_basis_leakage'],float(leakV)
        )
        if not s1:
            aux_path_rows.append({
                'indices':[i,j,k],'after_CV_support':0,
                'after_first_CK_support':0,'after_second_CK_support':0,
                'final_path_norm':0.0,
            })
            continue

        # Exact scalar-channel pruning rule: after the second leg a final rank
        # 0/1 leg can return to J=0 only from J=0 or J=1.  Here the second leg
        # is C_b(K), so J=2 is removed only AFTER its complete action.
        s2,d2=KCOMP.C_K_component(
            s1,source_v,b,j,k,JMAX2
        )
        update_diag(diagmax,{
            'CK_outer_complete_basis_leakage':d2['outer_complete_basis_leakage'],
            'CK_internal_volume_sector_leakage':d2['internal_volume_sector_leakage'],
            'CK_complete_charge_basis_leakage':d2['complete_charge_basis_leakage'],
        })
        s2_scalar_relevant={key:amp for key,amp in s2.items() if key[2] in (0,2)}

        if not s2_scalar_relevant:
            aux_path_rows.append({
                'indices':[i,j,k],'after_CV_support':len(s1),
                'after_first_CK_support':len(s2),
                'after_first_CK_scalar_relevant_support':0,
                'after_second_CK_support':0,'final_path_norm':0.0,
            })
            continue

        s3,d3=KCOMP.C_K_component(
            s2_scalar_relevant,source_v,a,i,j,JMAX2
        )
        update_diag(diagmax,{
            'CK_outer_complete_basis_leakage':d3['outer_complete_basis_leakage'],
            'CK_internal_volume_sector_leakage':d3['internal_volume_sector_leakage'],
            'CK_complete_charge_basis_leakage':d3['complete_charge_basis_leakage'],
        })
        add(total,s3)
        aux_path_rows.append({
            'indices':[i,j,k],
            'after_CV_support':len(s1),
            'after_first_CK_support':len(s2),
            'after_first_CK_scalar_relevant_support':len(s2_scalar_relevant),
            'after_second_CK_support':len(s3),
            'final_path_norm':math.sqrt(norm2(s3)),
        })

    weights=weight_by_source_J(total)
    total2=sum(weights.values())
    scalar2=sum(v for J,v in ((float(k),x) for k,x in weights.items()) if abs(J)<1e-15)
    nonscalar2=max(0.0,total2-scalar2)
    nonscalar_fraction=nonscalar2/max(total2,1e-30)
    outnorm=math.sqrt(norm2(total))
    mspin=max_spin(total)

    passed=(
        len(total)>0
        and outnorm>1e-10
        and nonscalar_fraction<1e-8
        and diagmax['CV_complete_basis_leakage']<1e-9
        and diagmax['CK_outer_complete_basis_leakage']<1e-9
        and diagmax['CK_internal_volume_sector_leakage']<1e-9
        and diagmax['CK_complete_charge_basis_leakage']<1e-9
        and mspin<=3.5+1e-12
    )

    ranked=sorted(total.items(),key=lambda kv:abs(kv[1]),reverse=True)[:12]
    return {
        'status':'first real Peter-Weyl ordered Lorentzian K-K-V triple',
        'passed':bool(passed),
        'source_node':source_v,
        'ordered_edges':[a,b,c],
        'definition':'Tr_aux[C_a(K) C_b(K) C_c(V)] with right-to-left state action',
        'Jmax':JMAX2/2,
        'auxiliary_index_paths':8,
        'output_support':len(total),
        'output_norm':outnorm,
        'output_weight_by_source_J':weights,
        'final_nonscalar_J_weight_fraction':nonscalar_fraction,
        'max_spin_reached':mspin,
        'max_diagnostics':diagmax,
        'path_diagnostics':aux_path_rows,
        'largest_output_amplitudes':[
            {'abs_amp':abs(amp),'amp':[amp.real,amp.imag],'key':repr(key)}
            for key,amp in ranked
        ],
        'exact_pruning_used':'After the middle C(K), J=2 is discarded because the last rank-(0+1) C(K) cannot couple J=2 to final J=0.',
        'beta_note':'No beta or (1+beta^2) factor is inserted; this is the raw structural K-K-V core.',
        'hermiticity_note':'K=[V,H_E] is anti-Hermitian in the present convention. This ordered raw triple is not yet the final Hermitian Lorentzian Hamiltonian ordering.',
        'next_use':'If green, compute the six signed permutations on one tetrahedral face, then the full four-face epsilon sum before introducing the fixed Lorentzian prefactor.',
        'scope_note':'One real ordered triple only; not yet the epsilon-oriented H_L node operator and not an HDA closure claim.',
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--v',type=int,default=0)
    ap.add_argument('--a',type=int)
    ap.add_argument('--b',type=int)
    ap.add_argument('--c',type=int)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    out=run(args.v,args.a,args.b,args.c)
    text=json.dumps(out,indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text+'\n',encoding='utf-8')
    return 0 if out['passed'] else 1


if __name__=='__main__':
    raise SystemExit(main())
