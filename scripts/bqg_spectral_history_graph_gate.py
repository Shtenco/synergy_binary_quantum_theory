#!/usr/bin/env python3
"""Fail-closed spectral-history graph gate for a finite positive master constraint.

Given certified matrix moments mu_n = V^† M^n V, construct the finite Krylov
spectral quotient associated with the seed block V. When the cyclic Krylov
subspace is certified invariant, the quotient is exact for V^† f(M) V,
including the master heat kernel and the zero spectral projector.

This gate is deliberately LIMITED to the finite spectral statement. It never
self-certifies a theory-specific physical BQG history from inline booleans.
Physical promotion requires an independent evidence-linked certificate that
binds this spectral result to the complete master/HDA/domain/refinement/source
artifacts on the same habitat.

This gate never interprets master-constraint eigenvalues as particle masses or
physical frequencies. Its heat parameter sigma is projector/constraint flow,
not relational proper time.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping

import numpy as np

SCHEMA = "BQG_SPECTRAL_HISTORY_MOMENTS_V1"
RESULT_SCHEMA = "BQG_SPECTRAL_HISTORY_GRAPH_V1"


def _c(v):
    if isinstance(v, (int, float)):
        return complex(float(v), 0.0)
    if isinstance(v, (list, tuple)) and len(v) == 2 and all(isinstance(x, (int, float)) for x in v):
        return complex(float(v[0]), float(v[1]))
    if isinstance(v, Mapping):
        return complex(float(v.get("re", 0.0)), float(v.get("im", 0.0)))
    raise ValueError(f"bad complex value {v!r}")


def decode_array(obj) -> np.ndarray:
    a = np.asarray([[_c(x) for x in row] for row in obj], dtype=complex)
    if a.ndim != 2:
        raise ValueError("matrix data must be rectangular 2D")
    return a


def decode_matrix(obj) -> np.ndarray:
    a = decode_array(obj)
    if a.shape[0] != a.shape[1]:
        raise ValueError("moment matrices must be square")
    return a


def encode_matrix(a: np.ndarray, imag_tol: float = 1e-13):
    out = []
    for row in np.asarray(a):
        rr = []
        for z in row:
            z = complex(z)
            if abs(z.imag) <= imag_tol:
                rr.append(float(z.real))
            else:
                rr.append([float(z.real), float(z.imag)])
        out.append(rr)
    return out


def _herm(a):
    return 0.5 * (a + a.conj().T)


def _opnorm(a):
    if a.size == 0:
        return 0.0
    return float(np.linalg.norm(a, 2))


def _block_matrix(mus, r: int, shift: int = 0):
    return np.block([[mus[i + j + shift] for j in range(r + 1)] for i in range(r + 1)])


def _rank_psd(a, rtol):
    h = _herm(a)
    ev = np.linalg.eigvalsh(h)
    scale = max(1.0, float(np.max(np.abs(ev))) if ev.size else 1.0)
    cut = rtol * scale
    rank = int(np.sum(ev > cut))
    return rank, ev, cut


def _spectral_weights(T: np.ndarray, C: np.ndarray):
    ev, U = np.linalg.eigh(_herm(T))
    weights = []
    for j in range(len(ev)):
        u = U[:, j:j + 1]
        W = C.conj().T @ (u @ u.conj().T) @ C
        weights.append(_herm(W))
    return ev, U, weights


def construct(packet: Mapping):
    if packet.get("schema") != SCHEMA:
        raise ValueError(f"expected schema {SCHEMA}")
    depth = int(packet.get("depth", -1))
    if depth < 0:
        raise ValueError("depth must be >=0")

    raw = packet.get("moments", {})
    mus = {int(k): decode_matrix(v) for k, v in raw.items()}
    need = 2 * depth + 1
    missing = [k for k in range(need + 1) if k not in mus]
    if missing:
        raise ValueError(f"missing moments through mu_{need}: {missing}")
    d = mus[0].shape[0]
    if any(a.shape != (d, d) for a in mus.values()):
        raise ValueError("all moments must have the same boundary dimension")

    tol = packet.get("tolerances", {})
    herm_tol = float(tol.get("hermiticity", 1e-10))
    mu0_tol = float(tol.get("mu0_identity", 1e-10))
    rank_rtol = float(tol.get("hankel_rank_rtol", 1e-11))
    psd_tol = float(tol.get("psd", 1e-10))
    moment_tol = float(tol.get("moment_reproduction", 5e-9))
    zero_tol = float(tol.get("zero_eigenvalue", 1e-10))
    residual_tol = float(tol.get("termination_residual", 1e-10))

    herm_errors = {str(k): _opnorm(a - a.conj().T) for k, a in mus.items()}
    herm_ok = max(herm_errors.values(), default=0.0) <= herm_tol
    mu0_error = _opnorm(mus[0] - np.eye(d))
    mu0_ok = mu0_error <= mu0_tol

    # Scale powers before the moment-Hankel construction to reduce conditioning.
    moment_scale = max(1.0, _opnorm(mus[1])) if need >= 1 else 1.0
    nu = {k: mus[k] / (moment_scale ** k) for k in mus}
    G = _herm(_block_matrix(nu, depth, 0))
    H = _herm(_block_matrix(nu, depth, 1))

    gev, GU = np.linalg.eigh(G)
    gscale = max(1.0, float(np.max(np.abs(gev))) if gev.size else 1.0)
    gcut = rank_rtol * gscale
    support = gev > gcut
    rank = int(np.sum(support))
    hankel_psd = bool(float(np.min(gev)) >= -psd_tol * gscale)
    if rank == 0:
        raise ValueError("zero-rank moment Hankel matrix")

    Z = GU[:, support] / np.sqrt(gev[support])
    T_scaled = _herm(Z.conj().T @ H @ Z)
    T = _herm(moment_scale * T_scaled)
    g0 = np.vstack([nu[k] for k in range(depth + 1)])
    C = Z.conj().T @ g0

    # Block-Gauss/Krylov moment reproduction on the supplied depth.
    reproduction = []
    max_rep = 0.0
    for n in range(2 * depth + 2):
        approx = C.conj().T @ np.linalg.matrix_power(T, n) @ C
        ref = mus[n]
        rel = _opnorm(approx - ref) / max(1.0, _opnorm(ref))
        max_rep = max(max_rep, rel)
        reproduction.append({"n": n, "relative_operator_error": rel})
    reproduction_ok = max_rep <= moment_tol

    tev, TU, weights = _spectral_weights(T, C)
    tscale = max(1.0, float(np.max(np.abs(tev))) if tev.size else 1.0)
    quotient_psd = bool(float(np.min(tev)) >= -psd_tol * tscale)

    sigmas = [float(x) for x in packet.get("heat_sigma", [0.0, 0.01, 0.1, 1.0, 10.0])]
    heat_rows = []
    for sigma in sigmas:
        if sigma < 0:
            raise ValueError("heat sigma must be >=0")
        heat = np.zeros((d, d), complex)
        p = 0.0
        lp = 0.0
        for lam, W in zip(tev, weights):
            fac = math.exp(-sigma * max(float(lam), 0.0))
            heat += fac * W
            trw = float(np.trace(W).real)
            p += fac * trw / d
            lp += max(float(lam), 0.0) * fac * trw / d
        ds = (2.0 * sigma * lp / p) if p > 0.0 else None
        heat_rows.append({
            "sigma": sigma,
            "seed_return_probability": p,
            "constraint_spectral_dimension": ds,
            "projected_heat": encode_matrix(_herm(heat)),
        })

    zero_weight = np.zeros((d, d), complex)
    spectral_rows = []
    for lam, W in zip(tev, weights):
        if abs(float(lam)) <= zero_tol * tscale:
            zero_weight += W
        spectral_rows.append({
            "eigenvalue": float(lam),
            "boundary_weight_trace": float(np.trace(W).real),
        })

    zeta_rows = []
    for s in [float(x) for x in packet.get("zeta_s", [0.5, 1.0, 2.0])]:
        val = 0.0
        for lam, W in zip(tev, weights):
            if float(lam) > zero_tol * tscale:
                val += float(np.trace(W).real) * float(lam) ** (-s)
        zeta_rows.append({"s": s, "seed_weighted_zeta": val})

    # Closure is fail-closed. Rank stabilization is only used when an upstream
    # producer explicitly certifies the moments strongly enough for that claim.
    rank_next = None
    rank_stable = False
    next_hankel_psd = None
    if 2 * (depth + 1) in mus:
        Gnext = _herm(_block_matrix(nu, depth + 1, 0))
        rank_next, evn, _ = _rank_psd(Gnext, rank_rtol)
        nscale = max(1.0, float(np.max(np.abs(evn))) if evn.size else 1.0)
        next_hankel_psd = bool(float(np.min(evn)) >= -psd_tol * nscale)
        rank_stable = rank_next == rank

    term = packet.get("termination_certificate", {})
    mode = str(term.get("mode", "none"))
    termination_ok = False
    termination_detail = {"mode": mode}
    if mode == "direct_block_residual":
        residual = float(term.get("residual_norm", math.inf))
        upstream = bool(term.get("certified", False))
        termination_ok = upstream and residual <= residual_tol
        termination_detail.update({"residual_norm": residual, "upstream_certified": upstream})
    elif mode == "certified_hankel_rank_stabilization":
        upstream = bool(term.get("certified", False))
        termination_ok = upstream and rank_stable and bool(next_hankel_psd)
        termination_detail.update({
            "upstream_certified": upstream,
            "rank_current": rank,
            "rank_next": rank_next,
            "rank_stable": rank_stable,
            "next_hankel_psd": next_hankel_psd,
        })

    finite_spectral_history_closed = bool(
        herm_ok and mu0_ok and hankel_psd and quotient_psd and reproduction_ok and termination_ok
    )

    # These are retained only as DECLARATIONS for diagnostics/backward
    # compatibility. They are never sufficient to certify physical history.
    pre = packet.get("physical_preconditions", {})
    declared_physical_preconditions = {
        "domain_complete": bool(pre.get("domain_complete", False)),
        "master_constraint_certified": bool(pre.get("master_constraint_certified", False)),
        "quantum_hda_or_explicit_dtarget_certified": bool(pre.get("quantum_hda_or_explicit_dtarget_certified", False)),
        "source_seed_complete_for_claim": bool(pre.get("source_seed_complete_for_claim", False)),
    }
    physical_history_candidate_ready = bool(
        finite_spectral_history_closed and all(declared_physical_preconditions.values())
    )

    # Critical fail-closed rule: this moment/quotient gate cannot verify those
    # declarations against independent artifacts, so it NEVER emits P_phys.
    physical_history_closed = False

    if physical_history_candidate_ready:
        status = "FINITE_SPECTRAL_HISTORY_CLOSED_PHYSICAL_CERTIFICATE_REQUIRED"
    elif finite_spectral_history_closed:
        status = "FINITE_SPECTRAL_HISTORY_CLOSED_BUT_PHYSICAL_PRECONDITIONS_OPEN"
    else:
        status = "SPECTRAL_HISTORY_GRAPH_OPEN"

    return {
        "schema": RESULT_SCHEMA,
        "status": status,
        "seed_label": packet.get("seed_label", "unspecified"),
        "depth": depth,
        "boundary_dimension": d,
        "moment_scale": moment_scale,
        "provenance": dict(packet.get("provenance", {})),
        "hankel": {
            "dimension": int(G.shape[0]),
            "rank": rank,
            "rank_next": rank_next,
            "rank_stable": rank_stable,
            "min_eigenvalue": float(np.min(gev)),
            "rank_cut": gcut,
            "positive_semidefinite": hankel_psd,
        },
        "quotient": {
            "dimension": int(T.shape[0]),
            "min_eigenvalue": float(np.min(tev)),
            "max_eigenvalue": float(np.max(tev)),
            "positive_semidefinite": quotient_psd,
            "operator_matrix": encode_matrix(T),
            "seed_embedding": encode_matrix(C),
            "spectral_nodes": spectral_rows,
        },
        "moment_checks": {
            "hermiticity_errors": herm_errors,
            "mu0_identity_error": mu0_error,
            "max_reproduction_relative_operator_error": max_rep,
            "reproduction": reproduction,
        },
        "termination": termination_detail,
        "finite_spectral_history_closed": finite_spectral_history_closed,
        "physical_preconditions": declared_physical_preconditions,
        "declared_physical_preconditions": declared_physical_preconditions,
        "physical_preconditions_verified": False,
        "physical_history_candidate_ready": physical_history_candidate_ready,
        "physical_history_closed": physical_history_closed,
        "physical_projector_emitted": False,
        "physical_projector_boundary_gram": None,
        "candidate_zero_spectral_weight": encode_matrix(_herm(zero_weight)),
        "heat_history": heat_rows,
        "spectral_diagnostics": {
            "weighted_zeta_positive_spectrum": zeta_rows,
            "note": "constraint/master spectral diagnostics only; not spacetime dimension, particle mass, physical omega, or relational time",
        },
        "claim_boundary": (
            "This gate certifies only the finite master spectral quotient/history. "
            "Inline physical_preconditions are unverified declarations and can never emit P_phys here. "
            "Physical BQG promotion requires an independent evidence-linked certificate binding this result "
            "to the complete master/HDA/domain/refinement/source artifacts on the same habitat."
        ),
    }


def _moments(M, V, nmax):
    out = {}
    X = V.copy()
    for n in range(nmax + 1):
        out[n] = V.conj().T @ X
        X = M @ X
    return out


def self_test():
    rng = np.random.default_rng(5092026)
    n, d = 6, 2
    Q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    lam = np.array([0.0, 0.2, 0.55, 1.1, 1.7, 2.4])
    M = _herm(Q @ np.diag(lam) @ Q.T)
    V, _ = np.linalg.qr(rng.normal(size=(n, d)))
    mus = _moments(M, V, 6)
    packet = {
        "schema": SCHEMA,
        "depth": 2,
        "seed_label": "synthetic_full_cyclic_control",
        "moments": {str(k): encode_matrix(v) for k, v in mus.items()},
        "termination_certificate": {"mode": "certified_hankel_rank_stabilization", "certified": True},
        "physical_preconditions": {
            "domain_complete": True,
            "master_constraint_certified": True,
            "quantum_hda_or_explicit_dtarget_certified": True,
            "source_seed_complete_for_claim": True,
        },
        "tolerances": {"moment_reproduction": 1e-8, "hankel_rank_rtol": 1e-10},
        "heat_sigma": [0.0, 0.1, 1.0, 10.0],
        "provenance": {"control": "synthetic_only"},
    }
    got = construct(packet)
    if not got["finite_spectral_history_closed"]:
        raise AssertionError(got)
    if not got["physical_history_candidate_ready"]:
        raise AssertionError("declared-precondition candidate readiness regression failed")
    if got["physical_history_closed"] or got["physical_projector_emitted"]:
        raise AssertionError("inline physical booleans incorrectly promoted finite spectral closure")

    T = decode_matrix(got["quotient"]["operator_matrix"])
    C = decode_array(got["quotient"]["seed_embedding"])
    te, tu = np.linalg.eigh(_herm(T))
    direct_e, direct_u = np.linalg.eigh(M)
    for sigma in (0.1, 1.0, 10.0):
        direct = V.T @ ((direct_u * np.exp(-sigma * direct_e)) @ direct_u.T) @ V
        approx = C.conj().T @ ((tu * np.exp(-sigma * te)) @ tu.conj().T) @ C
        if _opnorm(direct - approx) > 2e-9:
            raise AssertionError("heat mismatch")

    # Fail-closed negative: identical moments but no upstream termination cert.
    bad = dict(packet)
    bad["termination_certificate"] = {"mode": "none"}
    badgot = construct(bad)
    if badgot["finite_spectral_history_closed"] or badgot["physical_projector_emitted"]:
        raise AssertionError("missing termination certificate did not fail closed")

    # Fail-closed negative: finite spectral closure but one declared upstream
    # physical precondition is open. The finite quotient remains valid, but even
    # candidate readiness must be false.
    hda = json.loads(json.dumps(packet))
    hda["physical_preconditions"]["quantum_hda_or_explicit_dtarget_certified"] = False
    hgot = construct(hda)
    if not hgot["finite_spectral_history_closed"] or hgot["physical_history_candidate_ready"]:
        raise AssertionError("HDA declared-precondition fail-closed regression failed")
    if hgot["physical_history_closed"] or hgot["physical_projector_emitted"]:
        raise AssertionError("physical promotion bypass regression failed")

    return {
        "passed": True,
        "closed_control_quotient_dimension": got["quotient"]["dimension"],
        "closed_control_hankel_rank": got["hankel"]["rank"],
        "candidate_zero_weight_trace": float(np.trace(decode_matrix(got["candidate_zero_spectral_weight"])).real),
        "inline_physical_preconditions_do_not_emit_projector": True,
        "negative_missing_termination_status": badgot["status"],
        "negative_hda_status": hgot["status"],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--packet", type=Path)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        out = self_test()
    else:
        if args.packet is None:
            ap.error("--packet is required unless --self-test")
        packet = json.loads(args.packet.read_text(encoding="utf-8"))
        out = construct(packet)
    txt = json.dumps(out, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(txt, encoding="utf-8")
    print(txt, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
