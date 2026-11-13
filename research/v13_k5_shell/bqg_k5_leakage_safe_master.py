#!/usr/bin/env python3
"""Leakage-safe restricted K5 master built from all outgoing Peter-Weyl columns.

No physical projector is claimed without global constraint completeness,
shell/refinement stability, and adjoint/domain certification.
"""
from __future__ import annotations
import argparse,json,sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,vstack

HERE=Path(__file__).resolve().parent
SCRIPTS=HERE.parent.parent/"scripts"
for d in (HERE,SCRIPTS):
    if str(d) not in sys.path:sys.path.insert(0,str(d))

def assemble_columns(input_basis,vertices,image,signatures,dense_limit=256,
                     threshold=1e-9,reverse_image=None,reverse_limit=0):
    """Keep all output rows, including those outside the input shell."""
    inp=list(input_basis)
    if not inp or len(set(inp))!=len(inp):raise ValueError("empty or duplicate input shell")
    n=len(inp);pos={k:i for i,k in enumerate(inp)}
    cache={};all_outputs=set(inp)
    for v in vertices:
        for k in inp:
            out=image(k,v)
            if not isinstance(out,dict):raise TypeError("image must return dictionary")
            if not all(np.isfinite(a.real) and np.isfinite(a.imag) for a in out.values()):
                raise ValueError("non-finite amplitude")
            cache[(v,k)]=out
            all_outputs.update(ko for ko,a in out.items() if a!=0)
    targets=sorted(all_outputs,key=repr);tpos={k:i for i,k in enumerate(targets)}
    maps=[];stats=[];reverse_pairs=[]
    for v in vertices:
        row=[];col=[];dat=[];inside2=outside2=0.
        col_stats=[];reverse_candidates=[]
        for j,k in enumerate(inp):
            inside=outside=0.
            for ko,a in cache[(v,k)].items():
                if a==0:continue
                row.append(tpos[ko]);col.append(j);dat.append(complex(a))
                w=abs(a)**2
                if ko in pos:inside+=w;inside2+=w
                else:
                    outside+=w;outside2+=w
                    reverse_candidates.append((w,k,ko,complex(a)))
            col_stats.append({"input":j,"inside_norm2":inside,"outside_norm2":outside})
        C=coo_matrix((np.asarray(dat,dtype=complex),(row,col)),
                     shape=(len(targets),n)).tocsr()
        maps.append(C)
        stats.append({"vertex":v,"inside_norm2":inside2,"outside_norm2":outside2,
                      "outside_fraction":outside2/(inside2+outside2) if inside2+outside2 else 0.,
                      "columns_with_leakage":sum(x["outside_norm2"]>1e-24 for x in col_stats),
                      "column_statistics":col_stats})
        if reverse_image is not None and reverse_limit>0:
            for w,ki,ko,a in sorted(reverse_candidates,key=lambda x:-x[0])[:reverse_limit]:
                rev=reverse_image(ko,v).get(ki,0j)
                reverse_pairs.append({"vertex":v,"direct_norm2":w,
                   "direct_reversed_error":abs(rev-a.conjugate()),
                   "direct_amplitude":[a.real,a.imag],
                   "reverse_amplitude":[rev.real,rev.imag],
                   "input_repr":repr(ki),"outside_repr":repr(ko)})
    Cstack=vstack(maps,format="csr")
    master=(Cstack.conj().T@Cstack).tocsr()
    within_rows=[tpos[k] for k in inp]
    Cin=vstack([C[within_rows,:] for C in maps],format="csr")
    comp=(Cin.conj().T@Cin).tocsr()
    omitted=master-comp
    tf=float(master.diagonal().real.sum())
    tc=float(comp.diagonal().real.sum())
    result={"input_dimension":n,"output_union_dimension":len(targets),
       "external_output_dimension":len(targets)-n,"constraints":len(maps),
       "total_column_norm2":tf,"discarded_norm2_in_old_compression":tf-tc,
       "fraction_discarded":(tf-tc)/tf if tf else 0.,
       "per_vertex":stats,"reverse_pairs":reverse_pairs,
       "reverse_pairs_checked":len(reverse_pairs),
       "adjoint_definition":"Exact finite rectangular adjoint C†; global H_v adjoint not certified",
       "master_identity":"M_full=C_all†C_all; M_comp=C_inside†C_inside; M_full-M_comp=C_outside†C_outside >= 0",
       "spectrum":{"relative_threshold":threshold,"dense_limit":dense_limit},
       "status":"finite restricted master; not a physical BQG projector"}
    if n<=dense_limit:
        F=master.toarray();B=comp.toarray();O=omitted.toarray()
        ef,U=np.linalg.eigh(F);eb=np.linalg.eigvalsh(B)
        omin=float(np.linalg.eigvalsh(O)[0])
        tol=threshold*max(float(ef[-1]),1.)
        zero=ef<=tol;zb=eb<=tol
        P=U[:,zero]@U[:,zero].conj().T if zero.any() else np.zeros_like(F)
        groups=defaultdict(list)
        for j,k in enumerate(inp):groups[signatures(k)].append(j)
        blocks=[]
        for ga,ia in groups.items():
            for gb,ib in groups.items():
                if ga==gb:continue
                norm=float(np.linalg.norm(P[np.ix_(ib,ia)]))
                if norm>threshold:
                    blocks.append({"source":repr(ga),"target":repr(gb),"frobenius":norm})
        result["spectrum"].update({"computed":True,"full_min":float(ef[0]),
          "compressed_min":float(eb[0]),"full_nullity_at_threshold":int(zero.sum()),
          "compressed_nullity_at_threshold":int(zb.sum()),
          "false_null_modes_eliminated":int(zb.sum()-zero.sum()),
          "omitted_psd_min_eigenvalue":omin,
          "projector_residual_norm":float(np.linalg.norm(Cstack@P)),
          "full_nullity_threshold_scan":{str(r):int((ef<=r*max(float(ef[-1]),1.)).sum())
             for r in (1e-8,1e-10,1e-12)},
          "finite_shell_cross_graph_blocks":blocks,
          "projector_interpretation":"finite restricted kernel only, not physical projector"})
        if omin < -1e-8:raise AssertionError("omitted term not PSD")
        if zb.sum()<zero.sum():raise AssertionError("kernel monotonicity violated")
    else:
        result["spectrum"].update({"computed":False,
            "reason":"input shell exceeds dense_limit; no global projector asserted"})
    result["passed_structural_checks"]=bool(
        tf>=tc-1e-8 and len(targets)>=n and
        result["spectrum"].get("omitted_psd_min_eigenvalue",0.)>=-1e-8)
    return result

def run_k5(depth=0,jmax2=5,dense_limit=256,reverse_limit=0):
    import bqg_k5_graph_sector_master_shell as SHELL
    import k5_peter_weyl_safe_hda_column as K5
    initial,levels,basis=SHELL.generate_shell(depth,jmax2)
    out=assemble_columns(basis,K5.VERT,
       lambda k,v:SHELL.apply_H_to_key(k,v,jmax2),SHELL.graph_signature,
       dense_limit=dense_limit,
       reverse_image=(lambda ko,v:SHELL.apply_H_to_key(ko,v,jmax2)) if reverse_limit else None,
       reverse_limit=reverse_limit)
    out.update({"backend":"actual K5 Peter-Weyl H_v via existing BQG backend",
      "depth":depth,"jmax2":jmax2,"initial_state_repr":repr(initial),
      "level_sizes":[len(x) for x in levels],
      "regularization_note":"Existing K5 image routines prune small amplitudes; no exact symbolic certificate"})
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--depth",type=int,default=0)
    p.add_argument("--jmax2",type=int,default=5)
    p.add_argument("--dense-limit",type=int,default=256)
    p.add_argument("--reverse-limit",type=int,default=0)
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    result=run_k5(a.depth,a.jmax2,a.dense_limit,a.reverse_limit)
    text=json.dumps(result,indent=2,ensure_ascii=False)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(text+"\n",encoding="utf-8")
    print(text)
    return 0 if result["passed_structural_checks"] else 1
if __name__=="__main__":raise SystemExit(main())
