#!/usr/bin/env python3
"""Matrix-covariant Lorentzian K leg C_e(K)=h_e[h_e^-1,K_v].

This is the next genuine Thiemann factor after the verified C_e(V) and
K_v=[V_v,H_E,v].  The outer inverse holonomy creates one spectator fundamental
charge at each endpoint.  A complete gauge-invariant H_E move must preserve the
endpoint representation J=1/2, but an individual fixed-index primitive branch
need not do so.  Therefore NO primitive branch is projected early to J=1/2.

Every completed primitive is decomposed in the COMPLETE direct sum of total-J
recoupling sectors at the two charged endpoints.  All holonomy-index,
orientation and adjoint contributions are then summed.  Only on that full H_E
state is charge conservation tested:

    H_E : J=1/2 -> J=1/2.

Internal volume insertions are likewise evaluated without representation
truncation,

    V = direct_sum_J sqrt(|Q_J|),
    Q_J=P_J[J1.(J2xJ3)]P_J.

The charged extrinsic-curvature operator is K_v=[V_v,H_E,v].  The final forward
holonomy closes the target spectator and leaves the source matrix-covariant
J=0 plus J=1 content appropriate to a fundamental 2x2 operator.  Any surviving
J>1 output is a hard covariance failure.

No beta or HDA normalization is fitted.  This gate builds one exact covariant K
factor; the traced Lorentzian triple remains the next gate.
"""
from __future__ import annotations

import argparse
import functools
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import k5_peter_weyl_safe_hda_column as PW
import bqg_pruning_audit as PRUNE
import charged_intertwiner_recoupling_gate as CH
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_lorentzian_K_block_gate as KG


def add_state(dst, src, scale=1.0, tol=1e-11):
    for key, amp in src.items():
        prior = dst.get(key, 0j)
        incoming = scale * amp
        z = prior + incoming
        PRUNE.record_merge("covariant_K_leg.add_state", key, prior, incoming, z, tol)
        if PRUNE.keep(z,tol):
            dst[key] = z
        elif key in dst:
            del dst[key]


def label_gauss(K):
    return ("G", int(K))


def label_charge(J2, M2, K12, K34):
    return ("C", int(J2), int(M2), int(K12), int(K34))


def tensor_for_label(v, spins_local, label):
    if label[0] == "G":
        return PW.oriented_intertwiner(v, spins_local, label[1])
    if label[0] == "C":
        _, J2, M2, K12, K34 = label
        T = CH.charged_tensor(tuple(spins_local), K12, K34, J2, M2)
        return CV.orient_local(T, spins_local, v)
    raise ValueError(label)


def branch_from_key(key, amp=1+0j):
    spins, labels = key
    tensors = tuple(
        tensor_for_label(v, PW.local_spins(spins, v), labels[v])
        for v in PW.VERT
    )
    return spins, tensors, amp


def project_branch_complete_charges(branch, charged_nodes=(0, 1), tol=1e-11):
    spins, tensors, amp = branch
    opts_by_node = []
    max_basis_leak = 0.0
    for v in PW.VERT:
        ls = PW.local_spins(spins, v)
        X = tensors[v]
        opts = []
        recon = np.zeros_like(X)
        if v in charged_nodes:
            Xu = CV.unorient_local(X, ls, v)
            recon_u = np.zeros_like(Xu)
            for J2 in CV.all_total_J2(ls):
                rec = CH.allowed_charged_labels(tuple(ls), J2)
                for M2 in PW.m2vals_t(J2):
                    for K12, K34 in rec:
                        B = CH.charged_tensor(tuple(ls), K12, K34, J2, M2)
                        c = np.vdot(B, Xu)
                        if PRUNE.keep(c,1e-13):
                            opts.append((label_charge(J2, M2, K12, K34), c))
                            recon_u += c * B
                        else:
                            PRUNE.record_threshold("covariant_K_leg.charged_projection.1e-13",(spins,v,J2,M2,K12,K34),c,1e-13)
            recon = CV.orient_local(recon_u, ls, v)
        else:
            for K in PW.allowed_k2_t(*ls):
                B = PW.oriented_intertwiner(v, ls, K)
                c = np.vdot(B, X)
                if abs(c) > 1e-13:
                    opts.append((label_gauss(K), c))
                    recon += c * B
                else:
                    PRUNE.record_threshold("covariant_K_leg.gauss_projection.1e-13",(spins,v,K),c,1e-13)
        nrm = float(np.linalg.norm(X))
        leak = float(np.linalg.norm(X - recon) / max(nrm, 1e-30))
        max_basis_leak = max(max_basis_leak, leak)
        if not opts:
            return {}, max_basis_leak
        opts_by_node.append(tuple(opts))

    out = {}
    for choice in itertools.product(*opts_by_node):
        val = amp
        labels = []
        for lab, c in choice:
            labels.append(lab); val *= c
        key = (spins, tuple(labels))
        if PRUNE.keep(val,tol):
            out[key] = out.get(key, 0j) + val
        else:
            PRUNE.record_threshold("covariant_K_leg.complete_branch_projection",key,val,tol)
    for key, value in out.items():
        if abs(value) <= tol:
            PRUNE.record_threshold("covariant_K_leg.complete_branch_output",key,value,tol)
    return {k:a for k,a in out.items() if PRUNE.keep(a,tol)}, max_basis_leak


@functools.lru_cache(None)
def canonical_volume_block_general(spins_local, J2):
    spins_local = tuple(spins_local)
    qblocks = []
    for M2 in PW.m2vals_t(J2):
        _, Qb, _, _ = CH.q_block(spins_local, J2, M2)
        qblocks.append(Qb)
    Q = sum(qblocks) / len(qblocks)
    Q = 0.5 * (Q + Q.conj().T)
    return CH.canonical_volume_block(Q)


def apply_volume_allJ_oriented(T, spins_local, v):
    X = CV.unorient_local(T, spins_local, v)
    Y = np.zeros_like(X)
    recon = np.zeros_like(X)
    for J2 in CV.all_total_J2(spins_local):
        rec = CH.allowed_charged_labels(tuple(spins_local), J2)
        if not rec:
            continue
        Vb = canonical_volume_block_general(tuple(spins_local), J2)
        for M2 in PW.m2vals_t(J2):
            basis = [CH.charged_tensor(tuple(spins_local), a, b, J2, M2) for a,b in rec]
            coeff = np.asarray([np.vdot(B, X) for B in basis], complex)
            for c,B in zip(coeff,basis): recon += c*B
            outc = Vb @ coeff
            for c,B in zip(outc,basis): Y += c*B
    leak = float(np.linalg.norm(X-recon)/max(np.linalg.norm(X),1e-30))
    return CV.orient_local(Y, spins_local, v), leak


def apply_volume_complete_state(state, v):
    out = {}
    for key, amp in state.items():
        spins, labels = key
        lab = labels[v]
        if lab[0] != "C": raise ValueError("charged endpoint expected")
        _, J2, M2, K12, K34 = lab
        ls = PW.local_spins(spins, v)
        rec = CH.allowed_charged_labels(tuple(ls), J2)
        idx = rec.index((K12,K34))
        Vb = canonical_volume_block_general(tuple(ls), J2)
        for row,(A,B) in enumerate(rec):
            c = Vb[row,idx]
            labs = list(labels); labs[v] = label_charge(J2,M2,A,B)
            ko = (spins,tuple(labs))
            if PRUNE.keep(amp*c,1e-13):
                out[ko] = out.get(ko,0j)+amp*c
            else:
                PRUNE.record_threshold("covariant_K_leg.volume_matrix_element.1e-13",ko,amp*c,1e-13)
    for key, value in out.items():
        if abs(value)<=1e-11:
            PRUNE.record_threshold("covariant_K_leg.volume_state_output.1e-11",key,value,1e-11)
    return {k:a for k,a in out.items() if PRUNE.keep(a,1e-11)}


def apply_sequence_to_branch(branch, seq, source_v, Jmax2):
    branches=[branch]; max_v_leak=0.0
    for op in seq:
        if op[0]=="V":
            nb=[]
            for spins,tensors,amp in branches:
                t=list(tensors); ls=PW.local_spins(spins,source_v)
                t[source_v],leak=apply_volume_allJ_oriented(t[source_v],ls,source_v)
                max_v_leak=max(max_v_leak,leak)
                nb.append((spins,tuple(t),amp))
            branches=nb
        else:
            nb=[]
            for br in branches: nb.extend(PW.apply_path_branch(br,op[1],op[2],op[3],Jmax2))
            branches=nb
        if not branches: break
    return branches,max_v_leak


def apply_HE_complete_key(key, source_v, Jmax2, charged_nodes=(0,1)):
    base=branch_from_key(key); out={}; max_v_leak=0.0; max_basis_leak=0.0
    for sign,spec in PW.oriented_specs(source_v):
        v,a,b,c=spec
        for adj in (False,True):
            pref=0.5*sign
            for coef,seq0 in PW.T_sequences(v,a,b,c):
                seq=PW.adjoint_sequence(seq0) if adj else seq0
                branches,vleak=apply_sequence_to_branch(base,seq,source_v,Jmax2)
                max_v_leak=max(max_v_leak,vleak)
                for br in branches:
                    projected,bleak=project_branch_complete_charges(br,charged_nodes)
                    max_basis_leak=max(max_basis_leak,bleak)
                    add_state(out,projected,pref*coef)
    return out,max_v_leak,max_basis_leak


@functools.lru_cache(None)
def HE_complete_cached(key,source_v,Jmax2):
    out,vleak,bleak=apply_HE_complete_key(key,source_v,Jmax2)
    return tuple(out.items()),vleak,bleak


def apply_HE_complete_state(state,source_v,Jmax2):
    out={}; max_v=0.0; max_b=0.0
    for key,amp in state.items():
        items,vleak,bleak=HE_complete_cached(key,source_v,Jmax2)
        max_v=max(max_v,vleak); max_b=max(max_b,bleak)
        for ko,c in items: out[ko]=out.get(ko,0j)+amp*c
    for key, value in out.items():
        if abs(value)<=1e-10:
            PRUNE.record_threshold("covariant_K_leg.HE_state_output.1e-10",key,value,1e-10)
    return {k:a for k,a in out.items() if PRUNE.keep(a,1e-10)},max_v,max_b


def charge_weights(state,node):
    w={}
    for (_,labels),amp in state.items():
        lab=labels[node]
        if lab[0]=="C": w[lab[1]]=w.get(lab[1],0.0)+abs(amp)**2
    return w


def wrong_charge_fraction(state,node,expected_J2=1):
    w=charge_weights(state,node); total=sum(w.values())
    wrong=sum(x for J,x in w.items() if J!=expected_J2)
    return wrong/max(total,1e-30),{str(J/2):float(x) for J,x in sorted(w.items())}


def apply_K_complete_state(state,source_v,target_v,Jmax2):
    HE,vleak1,bleak1=apply_HE_complete_state(state,source_v,Jmax2)
    he_wrong_v,he_wv=wrong_charge_fraction(HE,source_v)
    he_wrong_w,he_ww=wrong_charge_fraction(HE,target_v)
    VH=apply_volume_complete_state(HE,source_v)
    Vstate=apply_volume_complete_state(state,source_v)
    HV,vleak2,bleak2=apply_HE_complete_state(Vstate,source_v,Jmax2)
    out={}; add_state(out,VH,+1); add_state(out,HV,-1)
    k_wrong_v,k_wv=wrong_charge_fraction(out,source_v)
    k_wrong_w,k_ww=wrong_charge_fraction(out,target_v)
    return out,max(vleak1,vleak2),max(bleak1,bleak2),{
        "HE_wrong_source_fraction":he_wrong_v,
        "HE_wrong_target_fraction":he_wrong_w,
        "HE_source_weights":he_wv,"HE_target_weights":he_ww,
        "K_wrong_source_fraction":k_wrong_v,
        "K_wrong_target_fraction":k_wrong_w,
        "K_source_weights":k_wv,"K_target_weights":k_ww,
    }


def inverse_outer_complete(initial,v,w,k,j,Jmax2):
    branches=[]
    for br in [PW.initial_factorized_oriented(initial)]: branches.extend(PW.apply_hit_branch(br,w,v,k,j,Jmax2))
    out={}; max_leak=0.0
    for br in branches:
        projected,leak=project_branch_complete_charges(br,(v,w))
        max_leak=max(max_leak,leak); add_state(out,projected)
    wrong_v,wv=wrong_charge_fraction(out,v); wrong_w,ww=wrong_charge_fraction(out,w)
    return out,max_leak,max(wrong_v,wrong_w),{"source":wv,"target":ww}


def close_complete_state_covariantly(state,v,w,i,k,Jmax2):
    out={}
    for key,amp in state.items():
        br=branch_from_key(key,amp)
        for cb in PW.apply_hit_branch(br,v,w,i,k,Jmax2): add_state(out,CV.project_covariant_branches([cb],v))
    return out


def covariant_K_leg(initial,v,w,i,j,Jmax2):
    total={}; diagmax={
        "outer_complete_basis_leakage":0.0,"outer_wrong_charge_fraction":0.0,
        "internal_volume_sector_leakage":0.0,"complete_charge_basis_leakage":0.0,
        "HE_wrong_charge_fraction":0.0,"K_wrong_charge_fraction":0.0,
    }; sample={}
    for k in range(2):
        inv,oleak,owrong,_=inverse_outer_complete(initial,v,w,k,j,Jmax2)
        Kstate,vleak,bleak,diag=apply_K_complete_state(inv,v,w,Jmax2)
        diagmax["outer_complete_basis_leakage"]=max(diagmax["outer_complete_basis_leakage"],oleak)
        diagmax["outer_wrong_charge_fraction"]=max(diagmax["outer_wrong_charge_fraction"],owrong)
        diagmax["internal_volume_sector_leakage"]=max(diagmax["internal_volume_sector_leakage"],vleak)
        diagmax["complete_charge_basis_leakage"]=max(diagmax["complete_charge_basis_leakage"],bleak)
        diagmax["HE_wrong_charge_fraction"]=max(diagmax["HE_wrong_charge_fraction"],diag["HE_wrong_source_fraction"],diag["HE_wrong_target_fraction"])
        diagmax["K_wrong_charge_fraction"]=max(diagmax["K_wrong_charge_fraction"],diag["K_wrong_source_fraction"],diag["K_wrong_target_fraction"])
        sample=diag
        add_state(total,close_complete_state_covariantly(Kstate,v,w,i,k,Jmax2))
    return total,diagmax,sample


def covariant_norm2(state): return float(sum(abs(a)**2 for a in state.values()))
def matrix_norm(M): return math.sqrt(sum(covariant_norm2(s) for row in M for s in row))


def run(v=0,w=1):
    JMAX2=5; initial=PW.basis_full_jhalf()[0]
    Kgauss=KG.apply_K_local({initial:1+0j},v,JMAX2); Kcov=CV.gauss_to_covariant(Kgauss,v)
    C=[[{} for _ in range(2)] for _ in range(2)]
    diagmax={"outer_complete_basis_leakage":0.0,"outer_wrong_charge_fraction":0.0,"internal_volume_sector_leakage":0.0,"complete_charge_basis_leakage":0.0,"HE_wrong_charge_fraction":0.0,"K_wrong_charge_fraction":0.0}; sample={}
    for i in range(2):
        for j in range(2):
            hKh,diag,smp=covariant_K_leg(initial,v,w,i,j,JMAX2)
            for name in diagmax: diagmax[name]=max(diagmax[name],diag[name])
            sample=smp
            out={}
            if i==j: add_state(out,Kcov,+1)
            add_state(out,hKh,-1); C[i][j]=out
    weights=CV.weight_by_J(C); total_weight=sum(weights.values()); j1=weights.get("1.0",0.0)
    high=sum(x for j,x in weights.items() if float(j)>1.0+1e-15)/max(total_weight,1e-30)
    Cnorm=matrix_norm(C); supports=[[len(C[i][j]) for j in range(2)] for i in range(2)]
    max_spin=max((max(key[0]) for row in C for state in row for key in state),default=0)/2
    passed=(len(Kgauss)>0 and diagmax["outer_complete_basis_leakage"]<1e-10 and diagmax["outer_wrong_charge_fraction"]<1e-18 and diagmax["internal_volume_sector_leakage"]<1e-10 and diagmax["complete_charge_basis_leakage"]<1e-10 and diagmax["HE_wrong_charge_fraction"]<1e-18 and diagmax["K_wrong_charge_fraction"]<1e-18 and Cnorm>1e-10 and j1>1e-14 and high<1e-18 and max_spin<=2.5+1e-12)
    return {"status":"matrix-covariant Peter-Weyl Lorentzian K leg with full charge-sector summation","passed":bool(passed),"edge":[v,w],"input":"all ten links j=1/2; all five K=0","Jmax":2.5,"Gauss_K_support":len(Kgauss),"Gauss_K_norm":math.sqrt(PW.norm2_state(Kgauss)),**diagmax,"sample_full_sum_charge_weights":sample,"C_matrix_supports":supports,"C_matrix_Frobenius_covariant_state_norm":Cnorm,"C_weight_by_source_J":weights,"C_J1_weight":j1,"C_J_greater_than_1_weight_fraction":high,"max_spin_reached":max_spin,"definition":"C_ij(K_v)=delta_ij K_v-sum_k h_ik K_v h^-1_kj with K_v=[V_v,H_E,v]","projection_rule":"Primitive branches are decomposed in all charged J sectors; J=1/2 conservation is tested only after the complete H_E/K sum.","beta_note":"No beta coefficient is inserted; the later Lorentzian triple receives the frozen classical prefactor without HDA fitting.","next_use":"If this gate passes, assemble the first epsilon^{ijk}Tr[C_i(K)C_j(K)C_k(V)] structural Lorentzian column.","scope_note":"One covariant K factor only; the traced triple, full H_L and Lorentzian HDA remain open."}


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("--v",type=int,default=0); ap.add_argument("--w",type=int,default=1); ap.add_argument("--output",type=Path); a=ap.parse_args()
    out=run(a.v,a.w); text=json.dumps(out,indent=2); print(text)
    if a.output: a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text+"\n",encoding="utf-8")
    return 0 if out["passed"] else 1
if __name__=="__main__": raise SystemExit(main())
