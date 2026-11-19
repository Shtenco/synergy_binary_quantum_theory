#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path
import numpy as np

TOL = 1e-10

def herm(M):
    return (M + M.conj().T) / 2

def dense_block_lanczos(M, Q0, max_depth=16, tol=TOL):
    Qblocks=[Q0]
    Ablocks=[]
    Bblocks=[]
    steps=[]
    closed=False
    for n in range(max_depth+1):
        Qn=Qblocks[n]
        MQ=M@Qn
        A=Qn.conj().T@MQ
        h_err=float(np.linalg.norm(A-A.conj().T))
        A=herm(A)
        R=MQ-Qn@A
        if n>0:
            R=R-Qblocks[n-1]@Bblocks[n-1].conj().T
        reorth=0.0
        for Qk in Qblocks:
            C=Qk.conj().T@R
            reorth=max(reorth,float(np.linalg.norm(C)))
            R=R-Qk@C
        G=herm(R.conj().T@R)
        ew,U=np.linalg.eigh(G)
        scale=max(float(np.max(np.abs(ew))) if ew.size else 0.0, float(np.linalg.norm(MQ))**2, 1.0)
        rank=int(np.sum(ew>tol*scale))
        residual_norm=float(np.sqrt(max(float(np.trace(G).real),0.0)))
        rel=residual_norm/max(float(np.linalg.norm(MQ)),1e-30)
        Ablocks.append(A)
        steps.append({"n":n,"block_rank":Qn.shape[1],"A_hermiticity_error":h_err,
                      "reorth_max_correction":reorth,"next_residual_norm":residual_norm,
                      "next_residual_relative":rel,"next_rank":rank,
                      "next_gram_eigenvalues":[float(x) for x in ew]})
        if rank==0:
            closed=True
            break
        keep=np.where(ew>tol*scale)[0]
        lam=ew[keep]
        Ur=U[:,keep]
        Qnext=R@Ur@np.diag(1/np.sqrt(lam))
        Bnext=np.diag(np.sqrt(lam))@Ur.conj().T
        Qblocks.append(Qnext)
        Bblocks.append(Bnext)
    return Qblocks,Ablocks,Bblocks,steps,closed

def assemble_J(Ablocks,Bblocks):
    dims=[A.shape[0] for A in Ablocks]
    off=np.cumsum([0]+dims)
    J=np.zeros((off[-1],off[-1]),complex)
    for n,A in enumerate(Ablocks):
        J[off[n]:off[n+1],off[n]:off[n+1]]=A
        if n < len(Bblocks) and n+1 < len(Ablocks):
            B=Bblocks[n]
            J[off[n+1]:off[n+2],off[n]:off[n+1]]=B
            J[off[n]:off[n+1],off[n+1]:off[n+2]]=B.conj().T
    return herm(J),off

def spectral_matrix_function(H, fn):
    w,U=np.linalg.eigh(herm(H))
    return (U*fn(w))@U.conj().T

def run_regression():
    rng=np.random.default_rng(260906)
    nreach=7
    r0=3
    X=rng.normal(size=(nreach,nreach))+1j*rng.normal(size=(nreach,nreach))
    U,_=np.linalg.qr(X)
    evals=np.array([0.0,0.17,0.42,0.9,1.8,3.1,5.0])
    Mreach=(U*evals)@U.conj().T
    M=np.zeros((10,10),complex)
    M[:nreach,:nreach]=Mreach
    M[nreach:,nreach:]=np.diag([6.0,7.0,8.0])
    Q0=np.eye(10,dtype=complex)[:,:r0]
    Q,A,B,steps,closed=dense_block_lanczos(M,Q0,max_depth=8,tol=1e-11)
    J,_=assemble_J(A,B)
    Qall=np.concatenate(Q[:len(A)],axis=1)
    compression_err=float(np.linalg.norm(J-Qall.conj().T@M@Qall))
    orth_err=float(np.linalg.norm(Qall.conj().T@Qall-np.eye(Qall.shape[1])))
    E0=np.zeros((J.shape[0],r0),complex); E0[:r0,:]=np.eye(r0)
    heat_errs={}
    for sigma in (0.01,0.1,1.0,5.0):
        direct=Q0.conj().T@spectral_matrix_function(M,lambda w:np.exp(-sigma*w))@Q0
        red=E0.conj().T@spectral_matrix_function(J,lambda w:np.exp(-sigma*w))@E0
        heat_errs[str(sigma)]=float(np.linalg.norm(direct-red))
    pd=Q0.conj().T@spectral_matrix_function(M,lambda w:(np.abs(w)<1e-9).astype(float))@Q0
    pr=E0.conj().T@spectral_matrix_function(J,lambda w:(np.abs(w)<1e-9).astype(float))@E0
    proj_err=float(np.linalg.norm(pd-pr))
    wj,Uj=np.linalg.eigh(J)
    seed=Uj[:r0,:]
    trace_weights=np.sum(np.abs(seed)**2,axis=0).real
    sum_rule=float(abs(np.sum(trace_weights)-r0))
    passed=bool(closed and Qall.shape[1]==nreach and orth_err<1e-9 and compression_err<1e-9
                and max(heat_errs.values())<1e-9 and proj_err<1e-9 and sum_rule<1e-9
                and np.min(wj)>-1e-9)
    return {"mode":"regression","passed":passed,"closed":closed,"reachable_dimension":nreach,
            "reduced_dimension":int(J.shape[0]),"orthogonality_error":orth_err,
            "compression_error":compression_err,"heat_kernel_errors":heat_errs,
            "zero_projector_error":proj_err,"spectral_weight_sum_rule_error":sum_rule,
            "J_min_eigenvalue":float(np.min(wj)),"steps":steps}

def _sparse_add(dst, src, scale=1.0, tol=1e-11):
    for k, a in src.items():
        z = dst.get(k, 0j) + scale * a
        if abs(z) > tol:
            dst[k] = z
        elif k in dst:
            del dst[k]

def _sparse_inner(a, b):
    if len(a) > len(b):
        a, b = b, a
    return sum(np.conj(v) * b.get(k, 0j) for k, v in a.items())

def _sparse_norm(a):
    return math.sqrt(float(sum(abs(v) ** 2 for v in a.values())))

def _sparse_gram(cols_a, cols_b=None):
    if cols_b is None:
        cols_b = cols_a
    G = np.zeros((len(cols_a), len(cols_b)), dtype=complex)
    for i, a in enumerate(cols_a):
        for j, b in enumerate(cols_b):
            G[i, j] = _sparse_inner(a, b)
    return G

def _combine_sparse(columns, coeffs, tol=1e-11):
    out = {}
    for col, c in zip(columns, coeffs):
        if abs(c) > tol:
            _sparse_add(out, col, c, tol=tol)
    return out

def run_peter_weyl_euclidean_first_edge(jmax2=5, tol=1e-10):
    """First actual block-Lanczos edge of the finite K5 Euclidean BQG master."""
    HERE = Path(__file__).resolve().parent
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    try:
        import k5_peter_weyl_safe_hda_column as PW
        import peter_weyl_logical_anisotropy_gate as AN
        import peter_weyl_euclidean_sine_ordering_gate as SINE
        import peter_weyl_zeroaware_volume_migration_experiment as ZVM
    except Exception as exc:
        return {"mode":"peter-weyl-euclidean-first-edge","passed":False,
                "science_status":"IMPORT_FAILURE","error":repr(exc)}

    ZVM.patch_and_clear()
    labels=[]
    q0=[]
    for env in AN.ENV_STATES:
        for pair in AN.PAIR_STATES:
            key=AN.logical_key(pair[0],pair[1],env)
            q0.append({key:1.0+0j})
            labels.append({"environment_K234":list(env),"pair_K01":list(pair)})
    q0_gram=_sparse_gram(q0)
    q0_orth_error=float(np.linalg.norm(q0_gram-np.eye(len(q0))))

    def apply_master(state):
        out={}
        for v in PW.VERT:
            hv=SINE.safe_H_sine(state,v,jmax2)
            h2v=SINE.safe_H_sine(hv,v,jmax2)
            _sparse_add(out,h2v)
        return out

    mq0=[apply_master(q) for q in q0]
    supports=[len(x) for x in mq0]
    norms=[_sparse_norm(x) for x in mq0]
    A0_raw=_sparse_gram(q0,mq0)
    A0_herm_error=float(np.linalg.norm(A0_raw-A0_raw.conj().T))
    A0=herm(A0_raw)
    a0_eigs=np.linalg.eigvalsh(A0).real
    mu2=herm(_sparse_gram(mq0))
    residual=[]
    for j,mq in enumerate(mq0):
        r=dict(mq)
        for i,qi in enumerate(q0):
            _sparse_add(r,qi,-A0[i,j])
        residual.append(r)
    R1_direct=herm(_sparse_gram(residual))
    R1_moment=herm(mu2-A0.conj().T@A0)
    moment_identity_error=float(np.linalg.norm(R1_direct-R1_moment))
    ew,U=np.linalg.eigh(R1_direct)
    scale=max(float(np.max(np.abs(ew))) if ew.size else 0.0,float(np.linalg.norm(mu2)),1.0)
    rank_tol=tol*scale
    keep=np.where(ew>rank_tol)[0]
    rank=int(len(keep))
    negative_min=float(np.min(ew)) if ew.size else 0.0
    residual_norm=math.sqrt(max(float(np.trace(R1_direct).real),0.0))
    mq_norm=math.sqrt(sum(x*x for x in norms))
    residual_relative=residual_norm/max(mq_norm,1e-30)
    recurrence_error=0.0
    q1_orth_error=0.0
    q0q1_error=0.0
    B1_shape=[0,len(q0)]
    if rank:
        lam=ew[keep]
        Ur=U[:,keep]
        q1=[]
        for a in range(rank):
            coeff=Ur[:,a]/math.sqrt(float(lam[a]))
            q1.append(_combine_sparse(residual,coeff))
        B1=np.diag(np.sqrt(lam))@Ur.conj().T
        B1_shape=list(B1.shape)
        q1_gram=_sparse_gram(q1)
        q0q1=_sparse_gram(q0,q1)
        q1_orth_error=float(np.linalg.norm(q1_gram-np.eye(rank)))
        q0q1_error=float(np.linalg.norm(q0q1))
        errs=[]
        for j,mq in enumerate(mq0):
            recon={}
            for i,qi in enumerate(q0):
                _sparse_add(recon,qi,A0[i,j])
            for a,qa in enumerate(q1):
                _sparse_add(recon,qa,B1[a,j])
            diff=dict(mq)
            _sparse_add(diff,recon,-1.0)
            errs.append(_sparse_norm(diff)/max(_sparse_norm(mq),1e-30))
        recurrence_error=max(errs,default=0.0)
    closure=(rank==0)
    algebra_pass=(q0_orth_error<1e-12 and A0_herm_error<2e-8
                  and float(np.min(a0_eigs))>-2e-8*max(float(np.max(np.abs(a0_eigs))),1.0)
                  and negative_min>-2e-8*scale
                  and moment_identity_error<2e-7*max(float(np.linalg.norm(R1_direct)),1.0)
                  and q1_orth_error<2e-7 and q0q1_error<2e-7 and recurrence_error<2e-7)
    return {
        "mode":"peter-weyl-euclidean-first-edge","passed":bool(algebra_pass),
        "science_status":"FINITE_EUCLIDEAN_MASTER_FIRST_SPECTRAL_EDGE",
        "operator":"M_E=sum_v (H_E,v^sine)^dagger H_E,v^sine = sum_v (H_E,v^sine)^2",
        "constraint_vertices":len(PW.VERT),"seed":"complete 32-state all-j=1/2 logical K5 sector",
        "seed_dimension":len(q0),"labels":labels,"jmax2":int(jmax2),"rank_tolerance":float(rank_tol),
        "q0_orthogonality_error":q0_orth_error,"master_image_support_min":min(supports,default=0),
        "master_image_support_max":max(supports,default=0),"master_image_norm_min":min(norms,default=0.0),
        "master_image_norm_max":max(norms,default=0.0),"A0_hermiticity_error":A0_herm_error,
        "A0_eigenvalue_min":float(np.min(a0_eigs)),"A0_eigenvalue_max":float(np.max(a0_eigs)),
        "mu2_minus_A0dagA0_identity_error":moment_identity_error,"R1_gram_min_eigenvalue":negative_min,
        "R1_gram_max_eigenvalue":float(np.max(ew)) if ew.size else 0.0,"B1_rank":rank,"B1_shape":B1_shape,
        "q1_orthogonality_error":q1_orth_error,"q0_q1_overlap_error":q0q1_error,
        "first_recurrence_relative_error":recurrence_error,"next_residual_relative":residual_relative,
        "closure_certified":bool(closure),
        "closure_statement":("B1=0: the 32-state seed space is invariant and the finite Euclidean seed history closes at depth 0"
                             if closure else "B1!=0: the master operator generates new directions; exact finite Euclidean history is NOT yet closed and block Lanczos must continue"),
        "claim_boundary":("Even if finite Euclidean Krylov closure is reached, this certifies only the chosen finite-regulator constraint spectral history. "
                          "It does not by itself close the Lorentzian/global refinement-compatible physical history, physical clock, connected W[J], or graviton propagator.")}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=("regression","peter-weyl-euclidean-first-edge"),default="regression")
    ap.add_argument("--jmax2",type=int,default=5)
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    out=run_regression() if args.mode=="regression" else run_peter_weyl_euclidean_first_edge(args.jmax2)
    text=json.dumps(out,indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text+"\n")
    return 0 if out["passed"] else 1

if __name__=="__main__":
    raise SystemExit(main())
