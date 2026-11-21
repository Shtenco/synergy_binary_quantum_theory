#!/usr/bin/env python3
"""Serialize and certify the complete 5x32 Euclidean q=2 boundary packet.

This is the reusable microscopic E-sector input for the BQG master assembler.
All five node constraints H_E,v are applied to all 32 frozen all-j=1/2 Gauss
boundary states with the same regulator-safe Peter-Weyl implementation used by
the preregistered five-node boundary master. Every sparse outgoing column is
written once, hashed for provenance, and the resulting M_EE Gram is audited
without re-running the microscopic action.

Acceptance is deliberately based on physics/numerics invariants of the packet,
not on a byte-level SHA of one JSON serialization.  In addition, the exact norm
of every vector discarded by sparse pruning is measured from the same unpruned
microscopic action.  These norms give a rigorous Frobenius/operator upper bound
on the difference between the retained M_EE and the unpruned M_EE, which is the
relevant numerical floor for later near-zero spectral claims.

The 32D boundary is deliberately marked domain_complete=false: a full-rank
compressed M_EE is a boundary diagnostic, not the full physical projector.
"""
from __future__ import annotations

import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_euclidean_constraint_column as COL
import bqg_constraint_master_assembler_gate as MASTER

LEGACY_NODE0_INPUT0_SHA="a5b5461cdaeedd1baf49dcfac881eda96e3d04cea182b8b6b639f5a6a585edbf"
TOL=2.0e-9


def canonical_packet_hash(rows):
    h=hashlib.sha256()
    for r in rows:
        h.update(f"{r['family']}:{r['node']}:{r['input_index']}:{r['sha256']}\n".encode())
    return h.hexdigest()


def sparse_inner(a,b):
    if len(a)>len(b): return np.conj(sparse_inner(b,a))
    return sum(np.conj(z)*b.get(k,0.0j) for k,z in a.items())


def sparse_norm(a): return float(np.sqrt(max(float(sparse_inner(a,a).real),0.0)))


def is_boundary_key(key):
    spins,Ks=key
    return all(int(s)==1 for s in spins) and all(int(k) in (0,2) for k in Ks)


def boundary_projection_norm(state): return sparse_norm({k:z for k,z in state.items() if is_boundary_key(k)})

def finite_state(state): return all(np.isfinite(complex(z).real) and np.isfinite(complex(z).imag) for z in state.values())


def run(outdir:Path,jmax2=5,prune=1e-8):
    basis=PW.basis_full_jhalf()
    if len(basis)!=32: raise RuntimeError(f"expected 32 boundary states, got {len(basis)}")
    coldir=outdir/'columns'; coldir.mkdir(parents=True,exist_ok=True)
    packet_rows=[];images={v:[] for v in range(5)};support=[];norms=[]
    retained_norm_by_node={v:[] for v in range(5)};discarded_norm_by_node={v:[] for v in range(5)}
    direct_return_max=0.0;all_finite=True;all_column_passed=True;parameter_consistent=True
    prune_pyth_max=0.0;max_discarded_amp=0.0;max_relative_discard=0.0;discarded_support_total=0

    for v in range(5):
        for i in range(32):
            payload=COL.run(v,i,jmax2,prune)
            p=coldir/f'E_node{v}_input{i:02d}.json';p.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
            rows=payload['complete_gauss_outgoing_column']['state'];st=MASTER.decode_state_rows(rows);pa=payload['pruning_audit']
            images[v].append(st);support.append(payload['complete_gauss_outgoing_column']['support']);norms.append(payload['complete_gauss_outgoing_column']['norm'])
            retained_norm_by_node[v].append(float(pa['retained_norm']));discarded_norm_by_node[v].append(float(pa['discarded_norm']))
            prune_pyth_max=max(prune_pyth_max,float(pa['pythagorean_relative_error']));max_discarded_amp=max(max_discarded_amp,float(pa['max_discarded_amplitude']))
            max_relative_discard=max(max_relative_discard,float(pa['relative_discarded_norm']));discarded_support_total+=int(pa['discarded_support'])
            direct_return_max=max(direct_return_max,boundary_projection_norm(st));all_finite=all_finite and finite_state(st);all_column_passed=all_column_passed and bool(payload.get('passed',False))
            parameter_consistent=parameter_consistent and abs(float(payload['Jmax'])-jmax2/2)<1e-15 and abs(float(payload['prune_threshold'])-prune)<1e-20
            packet_rows.append({'family':'E','node':v,'input_index':i,'path':str(p.relative_to(outdir)),'sha256':payload['column_sha256'],
                                'support':payload['complete_gauss_outgoing_column']['support'],'norm':payload['complete_gauss_outgoing_column']['norm'],
                                'discarded_norm':float(pa['discarded_norm']),'relative_discarded_norm':float(pa['relative_discarded_norm'])})

    identities=[(r['family'],r['node'],r['input_index']) for r in packet_rows]
    exact_domain_coverage=(len(packet_rows)==160 and len(set(identities))==160 and set(identities)=={('E',v,i) for v in range(5) for i in range(32)})

    MEE_raw=np.zeros((32,32),complex);node_rows=[];node_herm_max=0.0;prune_bound_total=0.0
    for v in range(5):
        G=MASTER.gram(images[v]);MEE_raw+=G;gh=float(np.linalg.norm(G-G.conj().T));node_herm_max=max(node_herm_max,gh)
        BF=float(np.sqrt(sum(x*x for x in retained_norm_by_node[v])));DF=float(np.sqrt(sum(x*x for x in discarded_norm_by_node[v])))
        bound=2.0*BF*DF+DF*DF;prune_bound_total+=bound
        node_rows.append({'node':v,'trace':float(np.trace(G).real),'frobenius_norm':float(np.linalg.norm(G)),'hermiticity_error':gh,
            'support_min':int(min(packet_rows[v*32+i]['support'] for i in range(32))),'support_max':int(max(packet_rows[v*32+i]['support'] for i in range(32))),
            'support_mean':float(np.mean([packet_rows[v*32+i]['support'] for i in range(32)])),
            'retained_column_map_frobenius_norm':BF,'discarded_column_map_frobenius_norm':DF,
            'unpruned_minus_retained_Gram_operator_norm_upper_bound':bound})

    raw_herm=float(np.linalg.norm(MEE_raw-MEE_raw.conj().T));MEE=.5*(MEE_raw+MEE_raw.conj().T);audit=MASTER.spectral_audit(MEE)
    ev=np.asarray(audit['eigenvalues'],float);scale=max(float(np.max(np.abs(ev))),1.0);traces=np.asarray([r['trace'] for r in node_rows],float)
    trace_rel_spread=float((np.max(traces)-np.min(traces))/max(abs(float(np.mean(traces))),1.0))
    eig_min=float(np.min(ev));eig_max=float(np.max(ev));unpruned_eig_min_lower=eig_min-prune_bound_total

    rng=np.random.default_rng(5092026);quadratic_checks=[];quad_max_rel=0.0
    for _ in range(5):
        c=rng.normal(size=32)+1j*rng.normal(size=32);c/=np.linalg.norm(c);lhs=float(np.vdot(c,MEE@c).real);rhs=0.0
        for v in range(5):
            acc={}
            for ci,img in zip(c,images[v]):
                for k,z in img.items(): acc[k]=acc.get(k,0.0j)+ci*z
            rhs+=float(sparse_inner(acc,acc).real)
        rel=abs(lhs-rhs)/max(abs(rhs),1e-30);quad_max_rel=max(quad_max_rel,rel);quadratic_checks.append({'matrix_form':lhs,'sum_image_norms':rhs,'relative_error':rel})

    checks={
        'exact_5x32_domain_coverage':bool(exact_domain_coverage),'all_column_gates_passed':bool(all_column_passed),'all_sparse_amplitudes_finite':bool(all_finite),
        'column_parameters_consistent':bool(parameter_consistent),'node_grams_hermitian':bool(node_herm_max<TOL),'master_hermitian_before_symmetrization':bool(raw_herm<TOL),
        'master_positive_semidefinite':bool(eig_min>-TOL*scale),'direct_boundary_return_zero_by_parity':bool(direct_return_max<1e-10),
        'node_trace_permutation_covariance':bool(trace_rel_spread<2e-9),'quadratic_form_identity':bool(quad_max_rel<2e-10),
        'pruning_pythagorean_identity':bool(prune_pyth_max<2e-12),'discarded_amplitudes_respect_threshold':bool(max_discarded_amp<=prune*(1+1e-12)),
        'full_rank_robust_under_pruning_bound':bool(unpruned_eig_min_lower>0.0),
    }

    manifest={
        'schema':'BQG_MICROSCOPIC_CONSTRAINT_PACKET_V2','status':'complete Euclidean q=2 boundary outgoing-column packet with invariant and pruning-error certification',
        'family':'E','domain_label':'q2_all_jhalf_K5_boundary','domain_dimension':32,'domain_complete':False,'nodes':[0,1,2,3,4],'Jmax':jmax2/2,'prune_threshold':prune,
        'columns':packet_rows,'packet_sha256':canonical_packet_hash(packet_rows),
        'M_EE':{'definition':'M_EE[i,j] = sum_v <H_E,v b_i | H_E,v b_j>','rank':audit['rank'],'nullity':audit['nullity'],'rank_tolerance':audit['rank_tolerance'],
            'hermiticity_error_before_symmetrization':raw_herm,'eigenvalue_min':eig_min,'eigenvalue_max':eig_max,'smallest_positive':audit['smallest_positive'],
            'condition_number_on_support':audit['condition_number_on_support'],'trace':float(np.trace(MEE).real),'frobenius_norm':float(np.linalg.norm(MEE)),
            'hash':MASTER.hash_arrays(MEE),'per_node':node_rows},
        'pruning_error_certificate':{
            'derivation':'For each node A=B+D, so ||A^dag A-B^dag B||_2 <= 2||B||_2||D||_2+||D||_2^2 <= 2||B||_F||D||_F+||D||_F^2; node bounds are summed.',
            'unpruned_minus_retained_M_EE_operator_norm_upper_bound':prune_bound_total,
            'relative_to_retained_max_eigenvalue':prune_bound_total/max(abs(eig_max),1e-300),
            'retained_eigenvalue_min':eig_min,'certified_unpruned_eigenvalue_min_lower_bound':unpruned_eig_min_lower,
            'full_rank_robust_under_bound':bool(unpruned_eig_min_lower>0.0),'total_discarded_support':discarded_support_total,
            'max_column_relative_discarded_norm':max_relative_discard,'max_discarded_amplitude':max_discarded_amp,'max_pythagorean_relative_error':prune_pyth_max,
            'near_zero_rule':'Any future low-eigenvalue/refinement claim approaching this certified operator-norm floor must rerun at a lower prune threshold before physical interpretation.'
        },
        'support_summary':{'min':int(min(support)),'max':int(max(support)),'mean':float(np.mean(support))},
        'norm_summary':{'min':float(min(norms)),'max':float(max(norms)),'mean':float(np.mean(norms))},'direct_boundary_projection_max_norm':direct_return_max,
        'node_trace_relative_spread':trace_rel_spread,'quadratic_form_regression':quadratic_checks,'quadratic_form_max_relative_error':quad_max_rel,
        'serialization_provenance':{'legacy_node0_input0_sha':LEGACY_NODE0_INPUT0_SHA,'observed_node0_input0_sha':packet_rows[0]['sha256'],
            'matches_legacy_serialization':packet_rows[0]['sha256']==LEGACY_NODE0_INPUT0_SHA,'hard_acceptance':False,
            'note':'SHA records serialized payload provenance only; physics acceptance is governed by invariant_checks and the explicit pruning-error certificate.'},
        'invariant_checks':checks,'passed':bool(all(checks.values())),
        'claim_boundary':'Complete reusable E-sector boundary packet with a rigorous pruning perturbation bound. domain_complete=false: no physical-projector claim follows from the compressed 32D spectrum.'
    }
    (outdir/'euclidean_packet_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8');return manifest


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,default=Path('euclidean_boundary_packet'));ap.add_argument('--jmax2',type=int,default=5);ap.add_argument('--prune',type=float,default=1e-8)
    a=ap.parse_args();out=run(a.output_dir,a.jmax2,a.prune);print(json.dumps({k:v for k,v in out.items() if k!='columns'},indent=2));return 0 if out['passed'] else 1
if __name__=='__main__': raise SystemExit(main())
