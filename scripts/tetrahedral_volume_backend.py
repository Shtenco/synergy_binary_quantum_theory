#!/usr/bin/env python3
"""Tetrahedrally covariant four-leg volume completion on charged sectors.

The historical magnetic/covariant implementation used q_123 on the first three
local legs.  On a Gauss J=0 four-valent intertwiner this is sufficient (up to an
overall orientation/sign normalization), but a Thiemann word inserts V after
holonomy hits where the source is charged and total J need not vanish.  There a
fixed q_123 extension selects a preferred fourth leg.

Define instead the normalized oriented tetrahedral grasping

  Q_tet = (1/4) sum_{r=0}^3 (-1)^r q_{all legs except r}

in the canonical local-slot orientation.  The vertex orientation multiplies all
four terms by one common sign and therefore drops out of V=sqrt(abs(Q_tet)).
The factor 1/4 preserves the already frozen Gauss-sector absolute-volume
normalization: on J=0, Q_tet is +/- q_123 and hence V_tet=V_123.

The implementation works blockwise in exact total-(J,M) recoupling sectors and
uses the same zero-aware spectral functional calculus as the v1.2 production
volume.  It never constructs the exponentially larger global spin-network
Hilbert space.
"""
from __future__ import annotations
import functools,itertools,math
from contextlib import contextmanager
import numpy as np
import k5_peter_weyl_safe_hda_column as PW
import charged_intertwiner_recoupling_gate as CH
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_lorentzian_K_block_gate as KG
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

NORM=0.25
HERMITICITY_TOL=1e-10
M_INDEPENDENCE_TOL=1e-10
RECONSTRUCTION_TOL=1e-10
TETRA_Q_DIAGNOSTICS={}
TETRA_VOLUME_DIAGNOSTICS={}

def apply_q_tetra(T,spins_local):
    spins_local=tuple(spins_local);mats=[PW.spin_mats_cached(s) for s in spins_local]
    Y=np.zeros_like(T,dtype=complex)
    for omitted in range(4):
        tri=tuple(r for r in range(4) if r!=omitted);pref=NORM*((-1)**omitted)
        for a,b,c in itertools.product(range(3),repeat=3):
            e=PW.EPS3[a,b,c]
            if not e:continue
            X=T
            # Axis applications commute across distinct legs. Reverse order only
            # avoids any accidental interpretation as same-leg operator order.
            for leg,axis in reversed(tuple(zip(tri,(a,b,c)))):
                X=PW.apply_axis_np(X,leg,mats[leg][axis])
            Y += pref*e*X
    return Y

@functools.lru_cache(maxsize=None)
def tetra_q_block(spins_local,J2,M2):
    spins_local=tuple(spins_local);labels=CH.allowed_charged_labels(spins_local,J2)
    basis=[CH.charged_tensor(spins_local,a,b,J2,M2) for a,b in labels]
    Q=np.zeros((len(labels),len(labels)),complex);max_abs=0.0;max_rel=0.0
    for j,B in enumerate(basis):
        QB=apply_q_tetra(B,spins_local)
        coeff=np.asarray([np.vdot(A,QB) for A in basis],complex)
        recon=sum((c*A for c,A in zip(coeff,basis)),np.zeros_like(QB))
        Q[:,j]=coeff;err=float(np.linalg.norm(QB-recon));rel=err/max(float(np.linalg.norm(QB)),1e-30)
        max_abs=max(max_abs,err);max_rel=max(max_rel,rel)
    herm=float(np.linalg.norm(Q-Q.conj().T,ord="fro"))
    scale=max(float(np.linalg.norm(Q,ord="fro")),1.0)
    if not math.isfinite(herm) or herm > HERMITICITY_TOL*scale:
        raise RuntimeError(f"raw tetrahedral Q block is not Hermitian before symmetrization: {herm}")
    if (not math.isfinite(max_rel) or not math.isfinite(max_abs)
            or max_rel > RECONSTRUCTION_TOL):
        raise RuntimeError(
            f"raw tetrahedral Q leaks out of its total-(J,M) sector: "
            f"absolute={max_abs}, relative={max_rel}"
        )
    TETRA_Q_DIAGNOSTICS[(spins_local,int(J2),int(M2))]={
        "raw_hermiticity_defect":herm,"raw_frobenius_norm":scale,
        "absolute_sector_leakage":max_abs,"relative_sector_leakage":max_rel}
    # Return the raw reduced block: the caller must compare M sectors before averaging/symmetrizing.
    return Q,max_abs,max_rel

@functools.lru_cache(maxsize=None)
def tetra_volume_block_general(spins_local,J2):
    spins_local=tuple(spins_local);m_values=tuple(PW.m2vals_t(J2));blocks=[]
    for M2 in m_values:blocks.append(tetra_q_block(spins_local,J2,M2)[0])
    Q0=blocks[0]
    scale=max(float(np.linalg.norm(Q0,ord="fro")),1.0)
    m_defects=[float(np.linalg.norm(Q-Q0,ord="fro")) for Q in blocks]
    max_m_defect=max(m_defects,default=0.0)
    if not math.isfinite(max_m_defect) or max_m_defect > M_INDEPENDENCE_TOL*scale:
        raise RuntimeError(f"raw tetrahedral Q block depends on M before averaging: {max_m_defect}")
    Q=sum(blocks)/len(blocks)
    # Symmetrization is now only a roundoff cleanup after raw Hermiticity passed.
    Q=0.5*(Q+Q.conj().T)
    eigenvalues=np.linalg.eigvalsh(Q)
    V=ZVM.zeroaware_sqrt_abs(Q)[0]
    Vherm=float(np.linalg.norm(V-V.conj().T,ord="fro"))
    Vmin=float(np.min(np.linalg.eigvalsh(V))) if V.size else 0.0
    V2=V@V; V4=V2@V2; Q2=Q@Q
    functional_error=float(np.linalg.norm(V4-Q2,ord="fro")/max(float(np.linalg.norm(Q2,ord="fro")),1.0))
    if (not np.isfinite(eigenvalues).all() or not np.isfinite(V).all()
        or Vherm>HERMITICITY_TOL*max(float(np.linalg.norm(V,ord="fro")),1.0)
        or Vmin < -HERMITICITY_TOL
        or functional_error > 1e-9):
        raise RuntimeError("tetrahedral volume spectral functional-calculus validation failed")
    TETRA_VOLUME_DIAGNOSTICS[(spins_local,int(J2))]={
        "M_values_doubled":list(m_values),
        "raw_Q_M_block_defects":m_defects,
        "max_raw_Q_M_block_defect":max_m_defect,
        "max_Q_sector_leakage_absolute":max(TETRA_Q_DIAGNOSTICS[(spins_local,int(J2),int(M))]["absolute_sector_leakage"] for M in m_values),
        "max_Q_sector_leakage_relative":max(TETRA_Q_DIAGNOSTICS[(spins_local,int(J2),int(M))]["relative_sector_leakage"] for M in m_values),
        "raw_Q_hermiticity_defects":[TETRA_Q_DIAGNOSTICS[(spins_local,int(J2),int(M))]["raw_hermiticity_defect"] for M in m_values],
        "Q_eigenvalues":[float(x) for x in eigenvalues],
        "Q_min_eigenvalue":float(np.min(eigenvalues)) if eigenvalues.size else 0.0,
        "Q_max_eigenvalue":float(np.max(eigenvalues)) if eigenvalues.size else 0.0,
        "V_eigenvalues":[float(x) for x in np.linalg.eigvalsh(V)],
        "V_min_eigenvalue":Vmin,"V_hermiticity_defect":Vherm,
        "V4_minus_Q2_relative_error":functional_error,
        "zero_aware_functional_calculus":True}
    return V

def tetra_volume_diagnostics(spins_local,J2):
    key=(tuple(spins_local),int(J2))
    tetra_volume_block_general(*key)
    return dict(TETRA_VOLUME_DIAGNOSTICS[key])

def tetra_charged_volume(spins_local):
    return tetra_volume_block_general(tuple(spins_local),1)

def apply_tetra_volume_tensor(T,spins_local):
    """Apply V_tet to an arbitrary local magnetic tensor via complete J blocks."""
    spins_local=tuple(spins_local)
    T=np.asarray(T,dtype=complex)
    expected_shape=tuple(s+1 for s in spins_local)
    if T.shape!=expected_shape:
        raise ValueError(f"volume tensor shape {T.shape} does not match {expected_shape}")
    if not np.isfinite(T.real).all() or not np.isfinite(T.imag).all():
        raise ValueError("volume input tensor contains non-finite values")
    Y=np.zeros_like(T,dtype=complex);recon=np.zeros_like(T,dtype=complex)
    for J2 in CV.all_total_J2(spins_local):
        labels=CH.allowed_charged_labels(spins_local,J2)
        if not labels:continue
        V=tetra_volume_block_general(spins_local,J2)
        for M2 in PW.m2vals_t(J2):
            basis=[CH.charged_tensor(spins_local,a,b,J2,M2) for a,b in labels]
            coeff=np.asarray([np.vdot(B,T) for B in basis],complex)
            for c,B in zip(coeff,basis):recon+=c*B
            out=V@coeff
            for c,B in zip(out,basis):Y+=c*B
    residual=float(np.linalg.norm(T-recon))
    relative=residual/max(float(np.linalg.norm(T)),1e-30)
    if not math.isfinite(relative) or relative>RECONSTRUCTION_TOL:
        raise RuntimeError(f"tetrahedral volume magnetic-basis reconstruction incomplete: {relative}")
    if not np.isfinite(Y.real).all() or not np.isfinite(Y.imag).all():
        raise RuntimeError("tetrahedral volume output contains non-finite values")
    return Y

def clear_tetra_caches():
    tetra_q_block.cache_clear();tetra_volume_block_general.cache_clear()
    TETRA_Q_DIAGNOSTICS.clear();TETRA_VOLUME_DIAGNOSTICS.clear()

@contextmanager
def install_tetrahedral_volume_backend():
    """Install only the volume extension; restore every historical function on exit."""
    old_apply=PW.apply_volume_tensor
    old_charged=CV.canonical_charged_volume
    old_general=CK.canonical_volume_block_general
    # Ensure the common spectral zero convention is active before replacing the
    # q-construction itself.
    ZVM.patch_and_clear();clear_tetra_caches()
    PW.apply_volume_tensor=apply_tetra_volume_tensor
    CV.canonical_charged_volume=tetra_charged_volume
    CK.canonical_volume_block_general=tetra_volume_block_general
    for obj in (getattr(KG,'local_volume_column',None),getattr(PW,'apply_H_cached',None),
                getattr(PW,'apply_H_cached_state',None),getattr(CK,'HE_complete_cached',None)):
        if hasattr(obj,'cache_clear'):obj.cache_clear()
    try:
        yield
    finally:
        PW.apply_volume_tensor=old_apply
        CV.canonical_charged_volume=old_charged
        CK.canonical_volume_block_general=old_general
        for obj in (getattr(KG,'local_volume_column',None),getattr(PW,'apply_H_cached',None),
                    getattr(PW,'apply_H_cached_state',None),getattr(CK,'HE_complete_cached',None)):
            if hasattr(obj,'cache_clear'):obj.cache_clear()
        clear_tetra_caches()
