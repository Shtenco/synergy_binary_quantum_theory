#!/usr/bin/env python3
"""Machine-readable D_target/HDA certificate for the BQG physical master.

The repository already defines the tangential target on the flux/vertex-smooth
habitat,

    D(k,l) f = - E_l . partial_{x_k} f,

with the simplex HH target

    {H(k),H(k')} f =
      (-E_k . partial_{x_k'} f + E_k' . partial_{x_k} f)/(3 V).

This gate does NOT invent another diffeomorphism operator.  It re-runs the
existing exact classical flux bridge and emits the current quantum-certification
state in a form that the production master assembler can consume fail-closed.

PASS means the ledger is internally consistent and the classical target gate
passes.  It deliberately does NOT mean that the graph-changing Lorentzian
quantum HDA has been certified.  Until a regulator-safe dual/habitat residual is
actually supplied, `quantum_habitat_residual_certified` remains false and the
certificate must not authorize a final physical projector.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import flux_habitat_diffeo_gate as FLUX

SCHEMA = "BQG_DTARGET_HDA_CERTIFICATE_V1"


def stable_hash(payload: dict) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def run() -> dict:
    classical = FLUX.run(seed=260809, samples=100, eps=1e-7)
    body = {
        "schema": SCHEMA,
        "status": "DTARGET_DEFINED_QUANTUM_HABITAT_OPEN",
        "classical_flux_target_defined": bool(classical["passed"]),
        "classical_flux_target": {
            "D_target": "D(k,l) f = - E_l^a partial_{x_k^a} f",
            "HH_target": "{H(k),H(k')} f = (-E_k.partial_{x_k'} f + E_k'.partial_{x_k} f)/(3V)",
            "flux_identity": "E_l = 3 V grad(lambda_l)",
            "regression": classical,
        },
        "regulator_safe_euclidean_HH_precursor": {
            "available": True,
            "Jmax": 2.5,
            "input": "all-j=1/2 K5 boundary, all K=0",
            "commutator_norm": 1.681559985798016,
            "output_support": 510,
            "interpretation": "support/cutoff precursor only; fixed-sector leakage is not the HDA anomaly",
        },
        "required_quantum_test": {
            "object": "Delta_HH^hab on the graph-changing dual/vertex-smooth habitat",
            "definition": "||([H[N],H[M]]' - i hbar D_target[beta]') Psi||_hab / (||[H[N],H[M]]' Psi||_hab + ||hbar D_target[beta]' Psi||_hab)",
            "requires_Lorentzian_total_constraint": True,
            "requires_cylindrical_reduction": True,
            "requires_graph_isomorphism_or_diffeomorphism_equivalence": True,
            "requires_recoupling_covariance": True,
        },
        "dtarget_included_in_master": False,
        "quantum_habitat_residual_certified": False,
        "quantum_habitat_residual_value": None,
        "certified_for_physical_projector": False,
        "fail_closed_reason": (
            "The classical D_target is fixed and tested, but the regulator-safe graph-changing Lorentzian dual-HH residual has not yet been certified on the same production habitat."
        ),
        "claim_boundary": (
            "This certificate freezes the existing D_target definition and truth status. It is not a quantum-HDA PASS and cannot authorize P_phys while quantum_habitat_residual_certified=false."
        ),
    }
    body["certificate_sha256"] = stable_hash(body)
    body["passed"] = bool(
        body["classical_flux_target_defined"]
        and body["quantum_habitat_residual_certified"] is False
        and body["certified_for_physical_projector"] is False
    )
    return body


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=Path("verification_results/DTARGET_HDA_CERTIFICATE.json"))
    args = ap.parse_args()
    out = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
