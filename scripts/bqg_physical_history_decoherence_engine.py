#!/usr/bin/env python3
"""Finite reference engine for BQG physical overlaps, histories and causal fronts.

All physical inputs (rigging-map inner product, physical states, clock, distance,
and retarded kernels) must be supplied by a future theory-specific construction.
Passing these reference tests does not close any physical BQG gate.
"""
from __future__ import annotations
from typing import Hashable, Mapping, Sequence
import json
import math
import numpy as np


def _hermitian_psd(name: str, matrix: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    if not math.isfinite(float(tol)) or tol < 0:
        raise ValueError(f"{name} tolerance must be finite and nonnegative")
    a = np.asarray(matrix, dtype=complex)
    if not np.all(np.isfinite(a.real)) or not np.all(np.isfinite(a.imag)):
        raise ValueError(f"{name} must contain only finite values")
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError(f"{name} must be square")
    if not np.allclose(a, a.conj().T, atol=tol, rtol=0):
        raise ValueError(f"{name} must be Hermitian")
    if np.linalg.eigvalsh((a + a.conj().T) / 2).min() < -tol:
        raise ValueError(f"{name} must be positive semidefinite")
    return a


def physical_overlap_matrix(states: Sequence[np.ndarray],
                            physical_inner_product: np.ndarray,
                            tol: float = 1e-10) -> dict:
    """Gram matrix G_ab=<E_a|E_b>_phys for supplied physical vectors.

    Positivity follows exactly from c†Gc=||sum_a c_a E_a||²_phys >= 0,
    provided the supplied physical form is positive semidefinite. Null vectors
    must be quotiented before interpreting states as distinct physical records.
    """
    eta = _hermitian_psd("physical_inner_product", physical_inner_product, tol)
    vectors = [np.asarray(v, dtype=complex).reshape(-1) for v in states]
    if any(not np.all(np.isfinite(v.real)) or not np.all(np.isfinite(v.imag)) for v in vectors):
        raise ValueError("states must contain only finite values")
    if not vectors or any(v.size != eta.shape[0] for v in vectors):
        raise ValueError("states must be nonempty and match physical metric dimension")
    G = np.array([[np.vdot(a, eta @ b) for b in vectors] for a in vectors])
    G = (G + G.conj().T) / 2
    eig = np.linalg.eigvalsh(G)
    cs_violation = max((abs(G[i,j])**2 - max(0.0, G[i,i].real*G[j,j].real)
                        for i in range(len(vectors)) for j in range(len(vectors))), default=0.0)
    return {
        "gram": G,
        "minimum_eigenvalue": float(eig.min()),
        "rank": int(np.linalg.matrix_rank(G, tol=tol)),
        "null_directions": int(np.count_nonzero(np.abs(eig) <= tol)),
        "cauchy_schwarz_violation": float(max(0.0, cs_violation)),
        "positive_semidefinite": bool(eig.min() >= -tol),
        "scope": "Positivity conditional on supplied physical inner product; no BQG rigging map or actual BQG physical states are constructed here."
    }


def validate_density(rho: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    r = _hermitian_psd("density matrix", rho, tol)
    if abs(np.trace(r) - 1) > tol:
        raise ValueError("density matrix must have unit trace")
    return r


def decoherence_functional(class_ops: Mapping[Hashable, np.ndarray],
                           rho: np.ndarray, tol: float = 1e-10):
    """D_ab=Tr(C_a rho C_b†); positivity is a Gram-matrix identity."""
    if not class_ops:
        raise ValueError("at least one class operator is required")
    r = validate_density(rho, tol)
    mats = [np.asarray(x, dtype=complex) for x in class_ops.values()]
    n = r.shape[0]
    if any(x.shape != (n,n) for x in mats):
        raise ValueError("class operators must match rho")
    if any(not np.all(np.isfinite(x.real)) or not np.all(np.isfinite(x.imag)) for x in mats):
        raise ValueError("class operators must contain only finite values")
    w, u = np.linalg.eigh(r)
    root = (u * np.sqrt(np.clip(w, 0, None))) @ u.conj().T
    xs = [x @ root for x in mats]
    D = np.array([[np.trace(xs[a] @ xs[b].conj().T)
                   for b in range(len(xs))] for a in range(len(xs))])
    D = (D + D.conj().T) / 2
    if np.linalg.eigvalsh(D).min() < -tol:
        raise ArithmeticError("decoherence functional lost strong positivity")
    return list(class_ops), D


def history_sector_report(paths: Sequence[tuple[int, Hashable, complex]],
                          relative_zero_tol: float = 0.0) -> dict:
    if not math.isfinite(float(relative_zero_tol)) or relative_zero_tol < 0:
        raise ValueError("relative_zero_tol must be finite and nonnegative")
    groups = {}
    for depth, sector, amplitude in paths:
        if isinstance(depth, (bool, np.bool_)) or not isinstance(depth, (int, np.integer)):
            raise ValueError("depth must be an integer, not a boolean or coerced numeric value")
        if depth < 0:
            raise ValueError("depth must be nonnegative")
        z = complex(amplitude)
        if not math.isfinite(z.real) or not math.isfinite(z.imag):
            raise ValueError("path amplitudes must contain only finite values")
        groups.setdefault((int(depth), sector), []).append(z)
    rows = []
    for (depth, sector), amps in sorted(groups.items(), key=lambda z:(z[0][0],repr(z[0][1]))):
        total, l1 = sum(amps, 0j), sum(abs(a) for a in amps)
        rows.append({"depth":depth, "sector":str(sector), "path_count":len(amps),
                     "amplitude":[total.real,total.imag], "l1":l1,
                     "cancellation_ratio":0.0 if l1 == 0 else 1-abs(total)/l1,
                     "coherently_nonzero":abs(total)>relative_zero_tol*l1})
    support_depth = min((d for d,_,_ in paths), default=None)
    coherent_depth = min((r["depth"] for r in rows if r["coherently_nonzero"]), default=None)
    return {"rows":rows,"first_support_depth":support_depth,
            "first_coherent_depth":coherent_depth,
            "warning":"Finite-precision coherent support depends on the declared tolerance."}


def channel_frontier_report(samples: Mapping[str, Sequence[tuple[float,float,complex]]],
                            response_floor: float) -> dict:
    if not math.isfinite(float(response_floor)) or response_floor < 0:
        raise ValueError("response_floor must be finite and nonnegative")
    out = {}
    for channel, triples in samples.items():
        by_t = {}
        for t,x,g in triples:
            t, x, z = float(t), float(x), complex(g)
            if not math.isfinite(t) or not math.isfinite(x):
                raise ValueError("time and distance must be finite")
            if not math.isfinite(z.real) or not math.isfinite(z.imag):
                raise ValueError("retarded-response samples must be finite")
            if t < 0 or x < 0:
                raise ValueError("time and distance must be nonnegative")
            by_t.setdefault(t,[]).append((x,abs(z)))
        front = [{"time":t,"max_detected_distance":max((x for x,a in vals if a>response_floor),default=None)}
                 for t,vals in sorted(by_t.items())]
        pts = [(z["time"],z["max_detected_distance"]) for z in front
               if z["max_detected_distance"] is not None and z["time"]>0]
        speeds = [(x1-x0)/(t1-t0) for (t0,x0),(t1,x1) in zip(pts,pts[1:]) if t1>t0]
        out[channel] = {"frontier":front,"max_observed_secant_speed":max(speeds) if speeds else None}
    v = [x["max_observed_secant_speed"] for x in out.values() if x["max_observed_secant_speed"] is not None]
    return {"channels":out,"cross_channel_spread":max(v)-min(v) if len(v)>1 else None,
            "status":"FINITE_THRESHOLD_DIAGNOSTIC_ONLY",
            "not_proven":["exact microcausality","physical distance and clock","continuum/refinement limit"]}


def refinement_channel_frontier_report(levels: Mapping[str, Mapping[str, Sequence[tuple[float,float,complex]]]],
                                       response_floor: float) -> dict:
    reports = []
    names = sorted({c for channels in levels.values() for c in channels})
    sequences = {c:[] for c in names}
    for level, channels in levels.items():
        r = channel_frontier_report(channels, response_floor)
        reports.append({"refinement":level,**r})
        for c in names:
            v = r["channels"].get(c,{}).get("max_observed_secant_speed")
            if v is not None:
                sequences[c].append(v)
    return {"levels":reports,"speeds_by_channel":sequences,
            "last_step_change":{c:(abs(v[-1]-v[-2]) if len(v)>1 else None) for c,v in sequences.items()},
            "status":"FINITE_REFINEMENT_DIAGNOSTIC_ONLY"}


def run_self_test() -> dict:
    # Physical Gram matrix control.
    eta = np.array([[2,.2j,0],[-.2j,1.5,.1],[0,.1,1]],complex)
    states = [np.array([1,0,1j]),np.array([.5j,1,1]),np.array([1,1j,.5])]
    gram = physical_overlap_matrix(states,eta)
    assert gram["positive_semidefinite"] and gram["cauchy_schwarz_violation"] < 1e-10

    # Both histories reach the same system output. Environment records are
    # orthogonal, and both branch weights are nonzero.
    k0=np.array([[1],[0]],complex); k1=np.array([[0],[1]],complex)
    F0=k0@k0.conj().T; F1=k0@k1.conj().T
    e0=np.array([[1,0],[0,0]],complex); e1=np.array([[0,0],[1,0]],complex)
    rho_s=np.full((2,2),.5,complex); rho_e=e0
    _, Denv=decoherence_functional({"a":np.kron(F0,e0),"b":np.kron(F1,e1)},np.kron(rho_s,rho_e))
    _, Dcoh=decoherence_functional({"a":F0,"b":F1},rho_s)
    assert abs(Denv[0,0].real-.5)<1e-10 and abs(Denv[1,1].real-.5)<1e-10
    assert abs(Denv[0,1])<1e-10 and abs(Dcoh[0,1].real-.5)<1e-10

    paths=[(1,"sector",1+1j),(1,"sector",-1-1j),(2,"sector",.2+.1j)]
    history=history_sector_report(paths)
    assert history["first_support_depth"]==1 and history["first_coherent_depth"]==2

    refinements={}
    for label,h in (("coarse",1.),("medium",.5),("fine",.25)):
        v=2+.2*h
        refinements[label]={c:[(.5,.5*v,1),(.5,.75*v,0),(1.,v,1),(1.,1.25*v,0),(1.5,1.5*v,1)]
                            for c in ("photon","matter","tensor")}
    cone=refinement_channel_frontier_report(refinements,1e-12)
    assert all(abs(x["cross_channel_spread"] or 0)<1e-10 for x in cone["levels"])
    return {
        "status":"PASS_REFERENCE_ONLY",
        "physical_overlap":{"minimum_eigenvalue":gram["minimum_eigenvalue"],"rank":gram["rank"],
                            "cauchy_schwarz_violation":gram["cauchy_schwarz_violation"]},
        "environment_decoherence_matrix":[[[z.real,z.imag] for z in row] for row in Denv],
        "coherent_decoherence_matrix":[[[z.real,z.imag] for z in row] for row in Dcoh],
        "history_support_depth":history["first_support_depth"],
        "history_coherent_depth":history["first_coherent_depth"],
        "refinement_fronts":cone,
        "physical_gates":"OPEN_PHYSICAL: no BQG rigging-map metric, physical clock, actual physical retarded kernels, or continuum cone supplied."
    }


if __name__ == "__main__":
    print(json.dumps(run_self_test(),indent=2,sort_keys=True))
