#!/usr/bin/env python3
"""Exact Fubini--Study symplectic gate on the q=2 C8 relational carrier.

This is a positive-control test, not the physical BQG gravity projector.
It constructs the rank-two invariant subspace of G=S_clock tensor J_geom,
embeds CP^1 into the 16-dimensional history Hilbert space, and verifies:
  * V^dagger V=I and P=VV^dagger is a rank-two orthogonal projector;
  * the pullback Fubini--Study form is (hbar/2) sin(theta) dtheta wedge dphi;
  * Bloch expectation functions obey the induced Poisson algebra;
  * the bracket agrees with the Pauli commutator using the explicit sign
    convention stated below.

Convention:
  omega_FS = 2 hbar Im <d_i psi|(1-P)|d_j psi> d lambda^i wedge d lambda^j/2
  {f,g}_FS = - omega_FS^{-1}(df,dg), so {q,p}=+1 for omega=dq wedge dp.
With this convention {n_a,n_b}=+(2/hbar) eps_abc n_c and
{f_A,f_B}=(1/(i hbar)) <[A,B]>.

The q=2 system step J is a declared positive-control step, not derived
physical gravitational evolution. The script does not close PHYSICAL_PROJECTOR_HISTORY.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import sympy as sp

def zero(M: sp.Matrix) -> bool:
    return all(sp.simplify(sp.expand_complex(v)) == 0 for v in M)

def run() -> dict:
    h, th, ph = sp.symbols("hbar theta phi", positive=True, real=True)
    S = sp.zeros(8)
    for t in range(8):
        S[(t + 1) % 8, t] = 1
    J = sp.Matrix([[0, -1], [1, 0]])
    V = sp.zeros(16, 2)
    for t in range(8):
        V[2*t:2*t+2, :] = J**t / sp.sqrt(8)
    P = sp.simplify(V * V.conjugate().T)
    I2 = sp.eye(2)
    checks = {
        "J8_equals_identity": zero(J**8-I2),
        "V_is_isometry": zero(V.conjugate().T*V-I2),
        "P_is_hermitian": zero(P-P.conjugate().T),
        "P_is_idempotent": zero(P*P-P),
        "P_has_rank_two": P.rank() == 2,
        "P_commutes_with_combined_step": zero(sp.kronecker_product(S,J)*P-P*sp.kronecker_product(S,J)),
    }

    # Generic normalized qubit and its embedding into relational histories.
    psi = sp.Matrix([
        sp.cos(th/2),
        sp.exp(sp.I*ph)*sp.sin(th/2)
    ])
    Psi = V*psi
    dth, dph = Psi.diff(th), Psi.diff(ph)
    Q = sp.eye(16)-Psi*Psi.conjugate().T
    omega_thph = sp.simplify(2*h*sp.im((dth.conjugate().T*Q*dph)[0]))
    omega_expected = h*sp.sin(th)/2
    checks["normalized_state"] = sp.simplify((Psi.conjugate().T*Psi)[0]-1) == 0
    checks["FS_pullback_exact"] = sp.simplify(omega_thph-omega_expected) == 0

    # Pauli expectations and the inverse-form Poisson convention.
    sx = sp.Matrix([[0,1],[1,0]])
    sy = sp.Matrix([[0,-sp.I],[sp.I,0]])
    sz = sp.Matrix([[1,0],[0,-1]])
    pauli = [sx, sy, sz]
    nvec = [sp.simplify((psi.conjugate().T*A*psi)[0]) for A in pauli]
    # omega_{theta phi}=a; our convention {f,g}=(f_theta*g_phi-f_phi*g_theta)/a.
    a = omega_expected
    brackets = {}
    commutator_residuals = {}
    eps = {(0,1,2):1,(1,2,0):1,(2,0,1):1,(1,0,2):-1,(2,1,0):-1,(0,2,1):-1}
    all_brackets = True
    for i in range(3):
        for j in range(i+1,3):
            b = sp.simplify((sp.diff(nvec[i],th)*sp.diff(nvec[j],ph)
                             - sp.diff(nvec[i],ph)*sp.diff(nvec[j],th))/a)
            k = 3-i-j
            target = sp.simplify(2*eps[(i,j,k)]*nvec[k]/h)
            brackets[f"{i}{j}"] = str(sp.trigsimp(b))
            all_brackets &= sp.trigsimp(b-target) == 0
            comm = pauli[i]*pauli[j]-pauli[j]*pauli[i]
            comm_expect = sp.simplify((psi.conjugate().T*comm*psi)[0]/(sp.I*h))
            commutator_residuals[f"{i}{j}"] = str(sp.trigsimp(b-comm_expect))
            all_brackets &= sp.trigsimp(b-comm_expect) == 0
    checks["Bloch_Poisson_matches_su2"] = bool(all_brackets)

    # Embedding preserves every expectation and commutator on the code subspace.
    embedded_expectations = []
    for A in pauli:
        Aemb = V*A*V.conjugate().T
        embedded_expectations.append(sp.simplify((Psi.conjugate().T*Aemb*Psi)[0] - (psi.conjugate().T*A*psi)[0]))
    checks["embedding_preserves_pauli_expectations"] = all(sp.simplify(x)==0 for x in embedded_expectations)

    return {
        "status": "exact q2 relational Fubini-Study positive-control gate",
        "passed": bool(all(checks.values())),
        "hilbert_dimension": 16,
        "invariant_subspace_dimension": 2,
        "projector_rank": P.rank(),
        "omega_theta_phi": str(omega_thph),
        "omega_expected": str(omega_expected),
        "bloch_expectations": [str(x) for x in nvec],
        "bloch_poisson_brackets": brackets,
        "commutator_residuals": commutator_residuals,
        "checks": checks,
        "claim_boundary": (
            "Exact symplectic and commutator correspondence on the finite q=2 "
            "relational positive-control subspace only. J is not the physical "
            "graph-changing gravity evolution operator; this does not establish "
            "the full BQG physical projector, continuum limit, TT gravitons, or Einstein IR."
        )
    }

def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=Path("verification_results/Q2_RELATIONAL_FS_SYMPLECTIC.json"))
    args=ap.parse_args()
    out=run()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))
    return 0 if out["passed"] else 1

if __name__=="__main__":
    raise SystemExit(main())
