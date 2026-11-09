#!/usr/bin/env python3
"""Aggregate complete W5 raw Lorentzian KKV shards into structural Gram matrices.

Requires every (node,K) pair and all 24 unique epsilon-term indices. Multiple
shard JSONs per (node,K) are summed coherently. This computes raw structural
G_EE, G_EL and G_LL; it does NOT certify Hermitian completion or physical H_L.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np

W5_K = (0, 2, 4, 6, 8)
NODES = (0, 1, 2, 3, 4)


def decode_sparse(rows):
    out = {}
    for key, re, im in rows:
        z = complex(float(re), float(im))
        out[key] = out.get(key, 0j) + z
    return out


def sparse_close(a, b, tol=1e-10):
    keys = set(a) | set(b)
    err = sum(abs(a.get(k, 0j) - b.get(k, 0j)) ** 2 for k in keys) ** 0.5
    scale = max(sum(abs(z) ** 2 for z in b.values()) ** 0.5, 1.0)
    return err <= tol * scale


def matrix_sparse(columns, row_keys):
    M = np.zeros((len(row_keys), 5), dtype=np.complex128)
    index = {key: i for i, key in enumerate(row_keys)}
    for col, state in enumerate(columns):
        for key, amp in state.items():
            M[index[key], col] = amp
    return M


def encode_complex(M):
    return [[[float(z.real), float(z.imag)] for z in row] for row in M]


def assemble_w5_gram(records):
    grouped = {(v, K): [] for v in NODES for K in W5_K}
    for rec in records:
        v, K = int(rec["source_vertex"]), int(rec["input_K"])
        if (v, K) not in grouped:
            raise ValueError(f"unexpected (node,K)=({v},{K})")
        grouped[(v, K)].append(rec)

    missing = [key for key, shards in grouped.items() if not shards]
    if missing:
        raise ValueError(f"missing W5 columns: {missing}")

    GEE = np.zeros((5, 5), dtype=np.complex128)
    GEL = np.zeros((5, 5), dtype=np.complex128)
    GLL = np.zeros((5, 5), dtype=np.complex128)
    node_meta = []
    common_cutoff = None

    for v in NODES:
        Acols, Bcols = [], []
        for K in W5_K:
            shards = grouped[(v, K)]
            term_ids = []
            Bsum = {}
            Aref = None
            for shard in shards:
                cutoff = int(shard["Jmax2"])
                if common_cutoff is None:
                    common_cutoff = cutoff
                if cutoff != common_cutoff:
                    raise ValueError("mixed Jmax2 cutoffs in one W5 Gram")
                if shard["physical_status"]["hermitian_completion_of_HL"]:
                    raise ValueError("unexpected physical-HL claim in raw KKV shard")
                Acur = decode_sparse(shard["A_Euclidean_sine_Gauss"])
                if Aref is None:
                    Aref = Acur
                elif not sparse_close(Acur, Aref):
                    raise ValueError(f"Euclidean A column differs across shards at node={v}, K={K}")
                Bpart = decode_sparse(shard["B_Lorentzian_raw_KKV_Gauss"])
                for key, amp in Bpart.items():
                    Bsum[key] = Bsum.get(key, 0j) + amp
                for row in shard["term_diagnostics"]:
                    term_ids.append(int(row["term_index"]))
                if shard["B_scalar_projection_rejections"]:
                    raise ValueError(f"non-scalar KKV outputs at node={v}, K={K}")
            if sorted(term_ids) != list(range(24)):
                raise ValueError(
                    f"epsilon terms incomplete/duplicated at node={v}, K={K}: "
                    f"{sorted(term_ids)}"
                )
            Acols.append(Aref or {})
            Bcols.append(Bsum)

        rows = sorted(set().union(*(set(x) for x in Acols), *(set(x) for x in Bcols)))
        Amat = matrix_sparse(Acols, rows)
        Bmat = matrix_sparse(Bcols, rows)
        GEE += Amat.conj().T @ Amat
        GEL += Amat.conj().T @ Bmat
        GLL += Bmat.conj().T @ Bmat
        node_meta.append({
            "node": v,
            "common_output_row_count": len(rows),
            "A_column_norms": np.linalg.norm(Amat, axis=0).tolist(),
            "B_raw_column_norms": np.linalg.norm(Bmat, axis=0).tolist(),
            "A_shape": list(Amat.shape),
            "B_shape": list(Bmat.shape),
        })

    GLE = GEL.conj().T
    herm = {
        "GEE": float(np.linalg.norm(GEE - GEE.conj().T)),
        "GLL": float(np.linalg.norm(GLL - GLL.conj().T)),
        "GLE_vs_GEL_dagger": float(np.linalg.norm(GLE - GEL.conj().T)),
    }
    return {
        "status": "COMPLETE_RAW_W5_STRUCTURAL_GRAM",
        "input_basis_K": list(W5_K),
        "nodes": list(NODES),
        "Jmax2": common_cutoff,
        "node_metadata": node_meta,
        "gram_matrices": {
            "G_EE": encode_complex(GEE),
            "G_EL_raw": encode_complex(GEL),
            "G_LE_raw": encode_complex(GLE),
            "G_LL_raw": encode_complex(GLL),
        },
        "spectra": {
            "G_EE_eigenvalues": np.linalg.eigvalsh((GEE + GEE.conj().T) / 2).tolist(),
            "G_LL_raw_eigenvalues": np.linalg.eigvalsh((GLL + GLL.conj().T) / 2).tolist(),
        },
        "hermiticity_diagnostics": herm,
        "claim_scope": "raw K_sine-K_sine-V structural maps only; no Hermitian H_L completion, no physical prefactor, no physical joint-kernel verdict",
        "next_required": [
            "prove/build Hermitian completion of the full Lorentzian operator",
            "fix beta/kappa/hbar and regulator normalization from the convention ledger",
            "repeat cutoff convergence and precision/error analysis",
            "only then form the physical full-constraint Gram and test its kernel",
        ],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", action="append", required=True, help="column/shard JSON; repeat for all files")
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    records = [json.loads(Path(p).read_text(encoding="utf-8")) for p in args.input]
    result = assemble_w5_gram(records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "Jmax2": result["Jmax2"],
        "G_EE_eigenvalues": result["spectra"]["G_EE_eigenvalues"],
        "G_LL_raw_eigenvalues": result["spectra"]["G_LL_raw_eigenvalues"],
        "hermiticity_diagnostics": result["hermiticity_diagnostics"],
        "output": str(args.output),
        "claim_scope": result["claim_scope"],
    }, indent=2))


if __name__ == "__main__":
    main()
