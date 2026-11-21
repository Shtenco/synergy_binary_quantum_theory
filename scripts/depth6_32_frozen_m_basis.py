#!/usr/bin/env python3
"""Construct a numerical [3,2] coordinate basis from frozen identity only.

This module intentionally never calls exact_mult and never enumerates the
Depth-6 shell or S5 orbit list.  The irrep multiplicity m and coordinate
dimension are immutable inputs from the schema-7 assignment identity ledger.
The selector backend is used only as a projector/eigenspace coordinate
constructor.  Every returned basis is checked fail-closed by J4, J5, projector
and orthonormality residuals.

A PASS here is a basis-coordinate certificate, not structural closure, actual-q
coverage, numerical rank evidence, or [3,2] master closure.
"""
from __future__ import annotations

import numpy as np


def build_frozen_m_basis(rep, key: str, frozen_m: int, frozen_coord_dim: int,
                         backend, *, extra: int = 0, seednum: int = 12345,
                         tol: float = 1e-8):
    rep=tuple(int(x) for x in rep)
    m=int(frozen_m); expected_d=int(frozen_coord_dim)
    if m<=0:
        raise RuntimeError('frozen m must be positive')
    if expected_d<=0:
        raise RuntimeError('frozen coord_dim must be positive')
    if key not in backend.CFG:
        raise RuntimeError(f'unknown frozen irrep key {key}')

    seed,j4,s4,j5,s5=backend.CFG[key]
    _,_,d=backend.layout(rep,seed)
    if int(d)!=expected_d:
        raise RuntimeError(f'coord_dim mismatch: backend={d} frozen={expected_d}')
    if m>d:
        raise RuntimeError(f'frozen m exceeds coordinate dimension: m={m} d={d}')

    k=min(d,m+int(extra))
    salt={'32':1,'311':2,'221':3,'2111':4}[key]
    rng=np.random.default_rng(seednum + 17*sum((i+1)*x for i,x in enumerate(rep)) + 1000003*salt)
    C=rng.standard_normal((d,k))+1j*rng.standard_normal((d,k))
    Y=backend.apply_selector_projector(rep,key,C)
    Q,R=np.linalg.qr(Y,mode='reduced')
    s=np.linalg.svd(R,compute_uv=False)
    rank_tol=max(Y.shape)*np.finfo(float).eps*max(float(s[0]) if len(s) else 1.0,1.0)*100
    rank=int((s>rank_tol).sum())
    if rank<m:
        raise RuntimeError(f'frozen-m selector rank {rank} < frozen m {m}')
    Q=Q[:,:m]

    e4=float(np.linalg.norm(backend.apply_sum(rep,seed,Q,backend.J4_TERMS)-j4*Q)/max(np.linalg.norm(Q),1e-30))
    e5=float(np.linalg.norm(backend.apply_sum(rep,seed,Q,backend.J5_TERMS)-j5*Q)/max(np.linalg.norm(Q),1e-30))
    ep=float(np.linalg.norm(backend.apply_selector_projector(rep,key,Q)-Q)/max(np.linalg.norm(Q),1e-30))
    oo=float(np.linalg.norm(Q.conj().T@Q-np.eye(m)))
    max_residual=max(e4,e5,ep,oo)
    if max_residual>float(tol):
        raise RuntimeError(
            f'frozen-m basis residual exceeds tolerance: max={max_residual} tol={tol} '
            f'j4={e4} j5={e5} projector={ep} orthonormal={oo}'
        )

    cert={
        'schema_version':1,
        'kind':'BQG_DEPTH6_FROZEN_M_NUMERICAL_BASIS_WITNESS',
        'status':'PASS_FROZEN_M_NUMERICAL_BASIS_ONLY',
        'irrep_key':key,
        'rep':list(rep),
        'frozen_m':m,
        'frozen_coord_dim':expected_d,
        'basis_shape':[int(Q.shape[0]),int(Q.shape[1])],
        'numerical_rank':rank,
        'rank_tolerance':float(rank_tol),
        'j4_residual':e4,
        'j5_residual':e5,
        'projector_residual':ep,
        'orthonormality_residual':oo,
        'max_residual':max_residual,
        'tolerance':float(tol),
        'multiplicity_recomputed':False,
        'shell_recomputed':False,
        'orbit_list_recomputed':False,
        'structural_closure_claimed':False,
        'actual_q_coverage_certified':False,
        'master_rank_certified':False,
        'claim_boundary':(
            'Numerical coordinate basis inside the already-frozen irrep identity only. '
            'No multiplicity, shell, orbit, structural-support, actual-q, master-rank or closure claim.'
        ),
    }
    return Q,cert
