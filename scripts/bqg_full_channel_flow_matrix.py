#!/usr/bin/env python3
"""Full BQG multiplicity-channel flow matrix across j -> j+1/2.

For every master eigenchannel a in the source [2,2] isotypic sector:
  1. tensor it with the full fresh binary-ancilla [2,2] singlet carrier;
  2. apply the exact microscopic symmetric edge-blocking tensor;
  3. extract the generated target [2,2] copy;
  4. compute its projector overlap with every isolated target master channel b.

The result is a row-stochastic overlap matrix F_ab (up to numerical error):
    F_ab = 1/2 Tr(P_block[a] P_master_target[b]).
This is NOT a transition probability in physical time.  It is a representation-
RG compatibility matrix between microscopic ancestry and target local master
eigenchannels.
"""
from __future__ import annotations
import argparse,itertools,json,math,sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0,str(HERE))

import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_euclidean_sine_ordering_gate as SINE
import peter_weyl_j1_s4_block_gate as J1
import peter_weyl_zeroaware_volume_migration_experiment as ZVM
from bqg_ancilla_peter_weyl_refinement_gate import inclusion

PERMS=tuple(itertools.permutations(range(4)))
CHI22={(1,1,1,1):2,(2,1,1):0,(2,2):2,(3,1):-1,(4,):0}

def inner(a,b):
    if len(a)>len(b): return np.conj(inner(b,a))
    return sum(np.conj(v)*b.get(k,0j) for k,v in a.items())

def gram(imgs):
    n=len(imgs); G=np.zeros((n,n),complex)
    for i in range(n):
        for j in range(i,n):
            z=inner(imgs[i],imgs[j]); G[i,j]=z; G[j,i]=np.conj(z)
    return (G+G.conj().T)/2

def validate_master_spectrum(eigenvalues, pair_spreads, context):
    """Require isolated, numerically degenerate two-dimensional master channels."""
    values = np.asarray(eigenvalues, dtype=float)
    spreads = np.asarray(pair_spreads, dtype=float)
    if values.ndim != 1 or values.size == 0 or spreads.shape != values.shape:
        raise RuntimeError(f"{context}: malformed master-channel spectrum")
    if not np.isfinite(values).all() or not np.isfinite(spreads).all():
        raise RuntimeError(f"{context}: nonfinite master-channel spectrum")
    if np.any(spreads < 0.0) or float(np.max(spreads)) >= 1e-8:
        raise RuntimeError(
            f"{context}: [2,2] eigenvalue-pair splitting exceeds 1e-8: {spreads.tolist()}"
        )
    gaps = np.diff(values)
    scale = max(float(np.max(np.abs(values))), 1.0)
    gap_tol = 1e-10 * scale
    if gaps.size and np.any(gaps <= gap_tol):
        raise RuntimeError(
            f"{context}: master channels are not spectrally isolated; "
            f"gaps={gaps.tolist()}, tolerance={gap_tol}"
        )
    return float(np.min(gaps)) if gaps.size else None


def validate_generated_copy_spectrum(singular_values, context):
    """Require the microscopic image to be exactly one numerical 2D copy."""
    s = np.asarray(singular_values, dtype=float)
    if s.ndim != 1 or s.size < 2 or not np.isfinite(s).all():
        raise RuntimeError(f"{context}: malformed generated-copy singular spectrum")
    rank = int(np.sum(s > 1e-10))
    if rank != 2:
        raise RuntimeError(
            f"{context}: generated target [2,2] copy has numerical rank {rank}, expected 2; "
            f"singular_values={s.tolist()}"
        )
    pair_spread = float(abs(s[0] - s[1]))
    if pair_spread >= 1e-8:
        raise RuntimeError(
            f"{context}: generated [2,2] pair splitting exceeds 1e-8: {pair_spread}"
        )
    return rank, pair_spread


def p22(s2):
    d=s2+1; P=np.zeros((d,d),complex)
    for p in PERMS:
        P += CHI22[J1.cycle_type(p)]*J1.permutation_matrix(s2,p)
    P*=2/24
    return (P+P.conj().T)/2

def master_channels(s2):
    spins=(s2,)*len(PW.EDGES)
    K2=tuple(range(0,2*s2+1,2))
    imgs=[]
    for K in K2:
        imgs.append(SINE.safe_H_sine({(spins,(K,0,0,0,0)):1+0j},0,s2+3))
    M=gram(imgs)
    Ms=np.zeros_like(M)
    for p in PERMS:
        U=J1.permutation_matrix(s2,p)
        Ms += U.conj().T@M@U
    Ms/=24; Ms=(Ms+Ms.conj().T)/2

    P=p22(s2)
    wp,Bfull=np.linalg.eigh(P)
    B=Bfull[:,wp>0.5]
    M22=B.conj().T@Ms@B
    M22=(M22+M22.conj().T)/2
    w,U=np.linalg.eigh(M22)
    order=np.argsort(w); w=w[order]; U=U[:,order]

    if len(w) % 2:
        raise RuntimeError(
            f"j={s2 / 2}: [2,2] isotypic master multiplicity is odd ({len(w)})"
        )
    channels=[]
    eigvals=[]
    spreads=[]
    for a in range(len(w)//2):
        pair=w[2*a:2*a+2]
        S=B@U[:,2*a:2*a+2]
        channels.append(S)
        eigvals.append(float(np.mean(pair)))
        spreads.append(float(np.ptp(pair)))
    min_channel_gap = validate_master_spectrum(
        eigvals, spreads, f"j={s2 / 2} master"
    )
    return {
        "channels":channels,
        "eigenvalues":eigvals,
        "pair_spreads":spreads,
        "minimum_interchannel_gap":min_channel_gap,
        "P22_basis":B,
        "P22":P,
    }

def block_tensor(old,anc,Bedge):
    return np.einsum(
        "Aae,Bbf,Ccg,Ddh,abcd,efgh->ABCD",
        Bedge,Bedge,Bedge,Bedge,old,anc,optimize=True
    )

def refinement_tensor(s2):
    _,old=J1.local_basis(s2)
    _,target=J1.local_basis(s2+1)
    _,anc=J1.local_basis(1)
    Bedge=inclusion(s2).conj().T.reshape(s2+2,s2+1,2)
    T=np.zeros((len(target),len(old),len(anc)),complex)
    leak=0.0
    for b,ob in enumerate(old):
        for a,aa in enumerate(anc):
            out=block_tensor(ob,aa,Bedge)
            coeff=np.array([np.vdot(t,out) for t in target],complex)
            T[:,b,a]=coeff
            proj=sum(coeff[r]*target[r] for r in range(len(target)))
            leak += float(np.vdot(out-proj,out-proj).real)
    return T,leak

def refinement_tensor_diagnostics(T, context="refinement tensor"):
    """Report the induced tensor's actual Gram/singular spectrum without assuming it is an isometry."""
    A = np.asarray(T, dtype=complex).reshape(T.shape[0], -1)
    if A.ndim != 2 or not np.isfinite(A.real).all() or not np.isfinite(A.imag).all():
        raise RuntimeError(f"{context}: tensor is non-finite or malformed")
    singular = np.linalg.svd(A, compute_uv=False)
    scale = max(float(singular[0]) if singular.size else 0.0, 1.0)
    rank_tol = 1e-12 * scale
    rank = int(np.sum(singular > rank_tol))
    row_gram = A @ A.conj().T
    col_gram = A.conj().T @ A
    return {
        "shape": list(A.shape),
        "numerical_rank": rank,
        "rank_tolerance": rank_tol,
        "singular_values": [float(x) for x in singular],
        "frobenius_norm": float(np.linalg.norm(A, ord="fro")),
        "row_gram_eigenvalues": [float(x) for x in np.linalg.eigvalsh((row_gram + row_gram.conj().T) / 2)],
        "column_gram_nonzero_eigenvalues": [float(x) for x in np.linalg.eigvalsh((col_gram + col_gram.conj().T) / 2) if x > rank_tol**2],
        "row_gram_identity_defect": float(np.linalg.norm(row_gram - np.eye(A.shape[0]), ord="fro")),
        "column_gram_idempotence_defect": float(np.linalg.norm(col_gram @ col_gram - col_gram, ord="fro")),
        "normalization_claim": "diagnostic only; no unit-isometry assumption is imposed on the induced singlet-to-singlet map",
    }


def generated_copy(T,source_basis,target_P22_basis):
    # source_basis: d_source x 2, ancilla logical dimension 2
    G=np.einsum("rba,bu->rua",T,source_basis).reshape(T.shape[0],4)
    B=target_P22_basis
    G22=B@(B.conj().T@G)
    U,s,Vh=np.linalg.svd(G22,full_matrices=False)
    rank, _ = validate_generated_copy_spectrum(s, "microscopic blocking")
    Q=U[:,:2]
    P=Q@Q.conj().T
    return P,s,rank

def run(s2):
    ZVM.patch_and_clear()
    src=master_channels(s2)
    dst=master_channels(s2+1)
    T,leak=refinement_tensor(s2)
    refinement_diagnostics = refinement_tensor_diagnostics(T, f"s2={s2}")
    if not math.isfinite(float(leak)) or float(leak) < 0.0:
        raise RuntimeError(f"s2={s2}: invalid refinement leakage {leak}")

    rows=[]
    singulars=[]
    row_sums=[]
    max_pair_spread=0.0
    generated_copy_ranks=[]

    dst_projectors=[S@S.conj().T for S in dst["channels"]]

    for a,Ssrc in enumerate(src["channels"]):
        Pblock,svals,generated_rank=generated_copy(T,Ssrc,dst["P22_basis"])
        singulars.append([float(x) for x in svals])
        generated_copy_ranks.append(generated_rank)
        if len(svals)>=2:
            max_pair_spread=max(max_pair_spread,float(abs(svals[0]-svals[1])))
        row=[]
        for Ptar in dst_projectors:
            row.append(float(np.trace(Pblock@Ptar).real/2.0))
        rows.append(row)
        row_sums.append(float(sum(row)))

    F=np.array(rows,float)
    # best target for each source channel
    best=np.argmax(F,axis=1)
    best_overlaps=np.max(F,axis=1)

    # overlap completeness: target master channels resolve target [22] isotypic sector
    stochastic_err=float(np.max(np.abs(np.array(row_sums)-1.0)))
    passed=(
        leak<1e-8
        and max_pair_spread<1e-8
        and all(rank == 2 for rank in generated_copy_ranks)
        and max(src["pair_spreads"], default=0.0)<1e-8
        and max(dst["pair_spreads"], default=0.0)<1e-8
        and stochastic_err<1e-8
        and np.min(F)>-1e-9
        and np.max(F)<1+1e-9
    )

    out={
      "status":"full microscopic/master multiplicity channel-flow matrix",
      "passed":bool(passed),
      "source_s2":s2,"source_j":s2/2,
      "target_s2":s2+1,"target_j":(s2+1)/2,
      "source_master_eigenvalues":src["eigenvalues"],
      "target_master_eigenvalues":dst["eigenvalues"],
      "source_pair_spreads":src["pair_spreads"],
      "target_pair_spreads":dst["pair_spreads"],
      "source_minimum_interchannel_gap":src["minimum_interchannel_gap"],
      "target_minimum_interchannel_gap":dst["minimum_interchannel_gap"],
      "channel_flow_matrix":F.tolist(),
      "row_sums":row_sums,
      "row_stochasticity_error":stochastic_err,
      "best_target_index_by_source":[int(x) for x in best],
      "best_overlap_by_source":[float(x) for x in best_overlaps],
      "generated_singular_values":singulars,
      "generated_copy_ranks":generated_copy_ranks,
      "generated_pair_spread_max":max_pair_spread,
      "Gauss_projection_leakage":leak,
      "refinement_tensor_diagnostics": refinement_diagnostics,
      "claim_boundary":"F_ab is a representation-RG projector overlap, not a physical-time probability. It diagnoses which target master eigenchannel is naturally inherited from each source channel under the exact binary-ancilla blocking tensor."
    }
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--s2",type=int,required=True)
    ap.add_argument("--output",type=Path)
    a=ap.parse_args()
    out=run(a.s2)
    txt=json.dumps(out,indent=2)
    print(txt)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(txt+"\n",encoding="utf-8")
    return 0 if out["passed"] else 2

if __name__=="__main__":
    raise SystemExit(main())
