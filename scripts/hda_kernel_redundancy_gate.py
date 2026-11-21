#!/usr/bin/env python3
"""Exact-kernel HDA redundancy theorem plus a near-zero counterexample.

For an exact common zero state of node constraints H_v, every commutator
[H_v,H_w] annihilates the state automatically. Therefore commutator-derived
D-like constraints do not further reduce the exact common kernel. Their role is
HDA/anomaly identification.

For approximate/near-zero sectors the implication needs a uniform control of
constraint operator norms (or a direct D residual test). A constructed Hermitian
sequence demonstrates why: M_H residual -> 0 while a commutator residual stays O(1)
when an operator norm diverges with the regulator.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np


def projector_from_master(M, rtol=1e-10):
    e, U = np.linalg.eigh(0.5*(M+M.conj().T))
    scale=max(float(np.max(np.abs(e))),1.0)
    z=np.abs(e)<=rtol*scale
    Q=U[:,z]
    return Q@Q.conj().T, e, int(z.sum())


def exact_control(seed=260905):
    rng=np.random.default_rng(seed)
    n=8; k=2
    X=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    U,_=np.linalg.qr(X)
    K=U[:,:k]; R=U[:,k:]
    constraints=[]
    for a in range(4):
        Y=rng.normal(size=(n-k,n-k))+1j*rng.normal(size=(n-k,n-k))
        A=0.5*(Y+Y.conj().T)+(2.0+0.3*a)*np.eye(n-k)
        constraints.append(R@A@R.conj().T)
    MH=sum(c.conj().T@c for c in constraints)
    comm=[]
    for i in range(len(constraints)):
        for j in range(i+1,len(constraints)):
            q=-1j*(constraints[i]@constraints[j]-constraints[j]@constraints[i])
            comm.append(0.5*(q+q.conj().T))
    MD=sum(q.conj().T@q for q in comm)
    P1,e1,r1=projector_from_master(MH)
    P2,e2,r2=projector_from_master(MH+MD)
    common_res=max(float(np.linalg.norm(c@P1)) for c in constraints)
    d_res=max(float(np.linalg.norm(q@P1)) for q in comm)
    return {
        "rank_exact_kernel_H_only":r1,
        "rank_exact_kernel_H_plus_commutator_D":r2,
        "projector_difference":float(np.linalg.norm(P1-P2)),
        "max_H_residual_on_kernel":common_res,
        "max_commutator_D_residual_on_kernel":d_res,
        "passed":bool(r1==k and r2==k and np.linalg.norm(P1-P2)<2e-9 and common_res<2e-9 and d_res<2e-9),
    }


def near_zero_negative_control():
    rows=[]
    good=True
    for eps in [1e-1,5e-2,2e-2,1e-2,5e-3]:
        # Hermitian constraints on e1,e2,e3. For psi=e1:
        # ||H1 psi||=eps, ||H2 psi||=eps, hence <M_H> = 2 eps^2 -> 0.
        # But H1 contains a 1/eps coupling e2<->e3, so [H1,H2] psi has O(1) norm.
        H1=np.array([[0,eps,0],[eps,0,1/eps],[0,1/eps,0]],dtype=float)
        H2=np.array([[0,0,eps],[0,0,0],[eps,0,0]],dtype=float)
        psi=np.array([1.,0.,0.])
        q=-1j*(H1@H2-H2@H1)
        mh=float(np.vdot(psi,(H1@H1+H2@H2)@psi).real)
        qres=float(np.linalg.norm(q@psi))
        opnorm=float(np.linalg.norm(H1,2))
        rows.append({"epsilon":eps,"H_master_expectation":mh,"commutator_residual":qres,"H1_operator_norm":opnorm})
    # Require the intended asymptotics: H-master falls by >100 while Q remains near one.
    good &= rows[-1]["H_master_expectation"] < rows[0]["H_master_expectation"]/100
    good &= min(r["commutator_residual"] for r in rows) > 0.9
    return {"passed":bool(good),"rows":rows,"lesson":"near-zero H master alone is insufficient if constraint norms are not uniformly controlled; require direct D/HDA residual convergence or a bounded-norm theorem"}


def run():
    ex=exact_control(); neg=near_zero_negative_control()
    return {
        "status":"exact HDA-kernel redundancy + near-zero fail-closed discriminator",
        "passed":bool(ex["passed"] and neg["passed"]),
        "exact_theorem":"if H_v psi=0 for every v, then [H_v,H_w] psi=0 for every pair; commutator-derived D cannot shrink the exact common kernel",
        "exact_control":ex,
        "near_zero_negative_control":neg,
        "production_rule":"exact-zero projector may use the complete H_v family with D as HDA validation; a refinement/near-zero projector must additionally demonstrate ||D P_low|| -> 0 (or an equivalent uniform bound), not merely lambda_low(M_H)->0"
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=Path('verification_results/HDA_KERNEL_REDUNDANCY.json'))
    a=ap.parse_args(); out=run(); a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8'); print(json.dumps(out,indent=2))
    return 0 if out['passed'] else 1

if __name__=='__main__': raise SystemExit(main())
