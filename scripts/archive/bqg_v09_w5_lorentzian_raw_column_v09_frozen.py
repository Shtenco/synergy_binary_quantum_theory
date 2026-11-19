#!/usr/bin/env python3
"""Actual W5 j=2 raw Lorentzian K-K-V column.

Unlike bqg_v08_lorentzian_24term_raw_column.py, this driver uses the exact
five W5 basis states used by the Euclidean certificate:
  spins=(4,)*10, K_labels=(K,0,0,0,0), K in {0,2,4,6,8}.
It returns both the sine-Hermitian Euclidean H_E column and the raw coherent
24-term Lorentzian K-K-V column in the same Gauss output-key convention.

IMPORTANT: raw K-K-V is not yet the Hermitian-completed full H_L and omits the
fixed beta/kappa/hbar prefactor. This driver therefore produces structural
cross-Gram evidence only, never a physical H_E+H_L kernel verdict.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
from pathlib import Path
import sys


def permutation_sign(base, permutation):
    ix = [base.index(q) for q in permutation]
    inversions = sum(ix[i] > ix[j] for i in range(3) for j in range(i + 1, 3))
    return -1 if inversions % 2 else 1


def ordered_epsilon_triples(neighbors):
    out = []
    for r, omitted in enumerate(neighbors):
        tri = tuple(x for x in neighbors if x != omitted)
        for p in itertools.permutations(tri):
            out.append(((-1) ** r * permutation_sign(tri, p), p))
    if len(out) != 24:
        raise AssertionError("K5 tetrahedral epsilon expansion must contain 24 ordered triples")
    return out


def w5_input_key(K):
    if K not in (0, 2, 4, 6, 8):
        raise ValueError("W5 K must be one of 0,2,4,6,8")
    return ((4,) * 10, (K, 0, 0, 0, 0))


def covariant_scalar_to_gauss(key, source_v):
    spins, Kother, J2, M2, K12, K34 = key
    if J2 != 0 or M2 != 0 or K12 != K34:
        raise ValueError(
            "raw epsilon-traced KKV output is not a Gauss scalar: "
            f"J2={J2}, M2={M2}, K12={K12}, K34={K34}"
        )
    labels = list(Kother)
    if len(labels) != 5 or labels[source_v] != -1:
        raise ValueError("malformed covariant sentinel or K5 vertex count")
    labels[source_v] = K12
    if any(x < 0 for x in labels):
        raise ValueError("unclosed non-source Gauss label in KKV output")
    return (tuple(spins), tuple(labels))


def norm_sparse(state):
    return math.sqrt(sum(abs(z) ** 2 for z in state.values()))


def add_sparse(dst, src, scale=1.0):
    for key, value in src.items():
        dst[key] = dst.get(key, 0j) + scale * value


def serialize_sparse(state):
    return [
        [repr(key), float(value.real), float(value.imag)]
        for key, value in sorted(state.items(), key=lambda item: repr(item[0]))
        if abs(value) > 1e-13
    ]


def run_column(source_v, input_K, Jmax2=10, start=0, count=24):
    if Jmax2 < 10:
        raise ValueError(
            "W5 j=2 KKV requires Jmax2>=10 as a conservative support wall; "
            "the old Jmax2=7 wall was preregistered for all-j=1/2 inputs"
        )
    if not (0 <= start < 24 and 1 <= count <= 24 - start):
        raise ValueError("invalid 24-term shard range")

    import k5_peter_weyl_safe_hda_column as PW
    import peter_weyl_euclidean_sine_ordering_gate as ESINE
    import peter_weyl_lorentzian_sine_ordered_triple_gate as LSINE
    import peter_weyl_lorentzian_ordered_triple_gate as RAW
    import peter_weyl_covariant_volume_leg_gate as CV
    import peter_weyl_zeroaware_volume_migration_experiment as ZVM

    ZVM.patch_and_clear()
    ini = w5_input_key(input_K)
    psi_gauss = {ini: 1 + 0j}
    psi_cov = CV.gauss_to_covariant(psi_gauss, source_v)

    # Euclidean column: same physical sine-Hermitian ordering as the W5 source.
    A = ESINE.safe_H_sine(psi_gauss, source_v, Jmax2)

    terms = ordered_epsilon_triples(PW.NEIG[source_v])
    B_cov = {}
    term_rows = []
    first_leg_cache = {}

    def executor():
        for term_index in range(start, min(start + count, 24)):
            epsilon, (a, b, c) = terms[term_index]
            term_out = {}
            max_leak = 0.0
            stage = {
                "after_CV_support_sum": 0,
                "after_CK1_support_sum": 0,
                "after_scalar_restriction_support_sum": 0,
                "after_CK2_support_sum": 0,
                "nonzero_CV_paths": 0,
                "nonzero_CK1_paths": 0,
                "nonzero_CK2_paths": 0,
            }
            for i, j, ell in itertools.product(range(2), repeat=3):
                cache_key = (c, ell, i)
                if cache_key not in first_leg_cache:
                    first_leg_cache[cache_key] = RAW.COMP.C_volume_component(
                        psi_cov, source_v, c, ell, i, Jmax2
                    )
                q, leak = first_leg_cache[cache_key]
                q = dict(q)
                max_leak = max(max_leak, float(leak))
                if q:
                    stage["after_CV_support_sum"] += len(q)
                    stage["nonzero_CV_paths"] += 1
                if not q:
                    continue

                q, d1 = RAW.KCOMP.C_K_component(q, source_v, b, j, ell, Jmax2)
                max_leak = max(
                    max_leak,
                    float(d1["outer_complete_basis_leakage"]),
                    float(d1["internal_volume_sector_leakage"]),
                    float(d1.get("complete_charge_basis_leakage", 0.0)),
                )
                if q:
                    stage["after_CK1_support_sum"] += len(q)
                    stage["nonzero_CK1_paths"] += 1

                # Exact scalar-channel selection rule for the final rank-(0+1) leg.
                q = {key: z for key, z in q.items() if key[2] in (0, 2)}
                stage["after_scalar_restriction_support_sum"] += len(q)
                if not q:
                    continue

                q, d2 = RAW.KCOMP.C_K_component(q, source_v, a, i, j, Jmax2)
                max_leak = max(
                    max_leak,
                    float(d2["outer_complete_basis_leakage"]),
                    float(d2["internal_volume_sector_leakage"]),
                    float(d2.get("complete_charge_basis_leakage", 0.0)),
                )
                if q:
                    stage["after_CK2_support_sum"] += len(q)
                    stage["nonzero_CK2_paths"] += 1
                add_sparse(term_out, q)

            add_sparse(B_cov, term_out, epsilon)
            row = {
                "term_index": term_index,
                "epsilon_sign": epsilon,
                "ordered_edges": [a, b, c],
                "output_support": len(term_out),
                "output_norm": norm_sparse(term_out),
                "max_basis_leakage": max_leak,
                "stage_counts": stage,
            }
            term_rows.append(row)
            print("W5_L_KKV_TERM_FINISHED " + json.dumps(row), flush=True)

    # Installs the exact same sine-H_E -> K_sine=[V,H_E^sine] engine used by
    # the all-jhalf raw KKV code; only the incoming W5 state and cutoff differ.
    LSINE.run(source_v, executor=executor)

    B_gauss = {}
    scalar_rejection = []
    for cov_key, amp in B_cov.items():
        try:
            gauss_key = covariant_scalar_to_gauss(cov_key, source_v)
        except ValueError:
            scalar_rejection.append(repr(cov_key))
            continue
        B_gauss[gauss_key] = B_gauss.get(gauss_key, 0j) + amp

    non_scalar_norm = math.sqrt(sum(abs(z) ** 2 for key, z in B_cov.items() if key[2] != 0))
    result = {
        "status": "COMPUTED_RAW_W5_KKV_COLUMN" if len(term_rows) == count else "PARTIAL_RAW_W5_KKV_COLUMN",
        "input_basis": "j=2 K5 all links; W5 basis K=(0,2,4,6,8) at vertex 0",
        "input_K": input_K,
        "input_key": repr(ini),
        "source_vertex": source_v,
        "Jmax2": Jmax2,
        "Jmax": Jmax2 / 2,
        "epsilon_terms_requested": count,
        "epsilon_terms_completed": len(term_rows),
        "term_diagnostics": term_rows,
        "A_Euclidean_sine_Gauss": serialize_sparse(A),
        "B_Lorentzian_raw_KKV_Gauss": serialize_sparse(B_gauss),
        "B_Lorentzian_raw_KKV_covariant": serialize_sparse(B_cov),
        "A_support": len(A),
        "B_support": len(B_gauss),
        "A_norm": norm_sparse(A),
        "B_raw_norm": norm_sparse(B_gauss),
        "B_non_scalar_covariant_norm": non_scalar_norm,
        "B_scalar_projection_rejections": scalar_rejection[:20],
        "output_basis_convention": "Gauss key (spins, K_labels); same node-local output Hilbert space for A and B",
        "normalization": {
            "A": "H_E^sine structural normalization from the existing sine-ordering gate",
            "B": "raw epsilon-summed K_sine-K_sine-V structural amplitude; no beta/kappa/hbar prefactor",
        },
        "physical_status": {
            "hermitian_completion_of_HL": False,
            "full_physical_HL": False,
            "full_H_E_plus_H_L": False,
            "physical_joint_kernel_verdict": False,
        },
        "scope_note": "This is a W5-input raw structural column. It is not yet the Hermitian Lorentzian Hamiltonian or a physical kernel test."
    }
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--node", type=int, default=0)
    ap.add_argument("--input-K", type=int, choices=(0,2,4,6,8), default=0)
    ap.add_argument("--Jmax2", type=int, default=10)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--count", type=int, default=24)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    result = run_column(args.node, args.input_K, args.Jmax2, args.start, args.count)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k not in (
        "A_Euclidean_sine_Gauss","B_Lorentzian_raw_KKV_Gauss",
        "B_Lorentzian_raw_KKV_covariant")}, indent=2))
    return 0 if result["epsilon_terms_completed"] == args.count else 1


if __name__ == "__main__":
    raise SystemExit(main())
