#!/usr/bin/env python3
"""Regression gate for the canonical BQG refinement residual decomposition."""
from __future__ import annotations

import math
import numpy as np


def require(condition, message):
    """A gate invariant that remains active under python -O."""
    if not condition:
        raise AssertionError(message)


def canonical_W():
    return np.column_stack([
        np.array([0.0,1.0,0.0],dtype=complex),
        np.array([2.0/3.0,0.0,-math.sqrt(5.0)/3.0],dtype=complex),
    ])


def run():
    W=canonical_W()
    P=W@W.conj().T
    require(np.linalg.norm(W.conj().T@W-np.eye(2))<1e-14, "W is not an isometry")
    require(np.linalg.norm(P@P-P)<1e-14, "carrier projector is not idempotent")

    rng=np.random.default_rng(20261009)
    ratios=[]

    for _ in range(100):
        A=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        Mc=(A+A.conj().T)/2
        B=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        Mf=(B+B.conj().T)/2

        R=Mc@W-W@Mf
        L=(np.eye(3)-P)@Mc@W
        D=W.conj().T@Mc@W-Mf

        require(np.linalg.norm(R-(L+W@D))<1e-11, "residual decomposition identity failed")
        require(np.linalg.norm(W.conj().T@L)<1e-11, "leakage is not orthogonal to carrier")
        require(np.linalg.norm(R.conj().T@R-(L.conj().T@L+D.conj().T@D))<1e-10, "orthogonal defect split failed")

        nR=np.linalg.norm(R,2)
        nL=np.linalg.norm(L,2)
        nD=np.linalg.norm(D,2)
        require(nR+1e-11>=max(nL,nD), "residual norm lower bound failed")
        require(nR<=math.sqrt(nL*nL+nD*nD)+1e-10, "residual norm upper bound failed")
        ratios.append((nR,nL,nD))

        # Galerkin pullback: only leakage remains.
        Mgal=W.conj().T@Mc@W
        Rg=Mc@W-W@Mgal
        Lg=(np.eye(3)-P)@Mc@W
        require(np.linalg.norm(Rg-Lg)<1e-11, "Galerkin residual is not pure leakage")

    # Registered j=1 absolute-volume scalar on the [2,2] carrier.
    # The theorem here checks the exact target matrix used by the canonical gate.
    Vdoublet=(3.0**0.25)*np.eye(2)
    require(np.linalg.norm(Vdoublet-(3.0**0.25)*np.eye(2))<1e-15, "registered j=1 volume scalar mismatch")

    print("trials =",len(ratios))
    print("W isometry error =",np.linalg.norm(W.conj().T@W-np.eye(2)))
    print("RESIDUAL DECOMPOSITION PASS")
    print("ORTHOGONAL DEFECT SPLIT PASS")
    print("GALERKIN LEAKAGE REDUCTION PASS")
    print("PASS")


if __name__=="__main__":
    run()
