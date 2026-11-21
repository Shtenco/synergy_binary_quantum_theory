#!/usr/bin/env python3
"""Build the first actual projected-heat/block-Krylov moment from BQG E columns.

For the orthonormal q=2 boundary injection V0 and the Euclidean constraint map
C_E=(H_E,v)_v, the positive master contribution is

    M_E = C_E^dagger C_E.

The complete 5 x 32 outgoing-column packet determines exactly

    mu_0 = V0^dagger V0 = I_32,
    mu_1 = V0^dagger M_E V0
         = (C_E V0)^dagger (C_E V0)
         = M_EE.

Therefore the actual projected heat kernel has the theory-specific short-tau
expansion

    K_E(tau) = V0^dagger exp(-tau M_E) V0
             = I_32 - tau mu_1 + O(tau^2).

This gate deliberately DOES NOT replace K_E(tau) by exp(-tau M_EE).  Higher
moments require further action of the enlarged graph-changing master and are
not present in the one-hit boundary packet.  Nor does positive mu_1 imply a
zero physical boundary projector: a boundary vector may have both zero-mode
and positive-spectrum components.
"""
from __future__ import annotations

import argparse, hashlib, json, math
from pathlib import Path
import numpy as np

SCHEMA_IN="BQG_MICROSCOPIC_CONSTRAINT_PACKET_V2"
SCHEMA_OUT="BQG_PROJECTED_HEAT_FIRST_MOMENT_V1"


def decode(rows):
    out={}
    for r in rows:
        key=(tuple(int(x) for x in r['spins']), tuple(int(x) for x in r['K_labels']))
        z=complex(float(r['amp'][0]),float(r['amp'][1]))
        out[key]=out.get(key,0j)+z
    return out


def inner(a,b):
    if len(a)>len(b):
        return np.conj(inner(b,a))
    return sum(np.conj(z)*b.get(k,0j) for k,z in a.items())


def canonical_matrix_hash(M):
    payload={
        'shape':list(M.shape),
        'real':np.asarray(M.real,dtype=np.float64).tolist(),
        'imag':np.asarray(M.imag,dtype=np.float64).tolist(),
    }
    raw=json.dumps(payload,sort_keys=True,separators=(',',':'),allow_nan=False).encode('utf-8')
    return hashlib.sha256(raw).hexdigest()


def run(packet_dir:Path):
    packet_dir=Path(packet_dir)
    manifest_path=packet_dir/'euclidean_packet_manifest.json'
    if not manifest_path.exists():
        raise RuntimeError(f'missing manifest {manifest_path}')
    manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
    if manifest.get('schema')!=SCHEMA_IN or manifest.get('family')!='E' or not manifest.get('passed',False):
        raise RuntimeError('invalid Euclidean microscopic packet manifest')
    if int(manifest.get('domain_dimension',-1))!=32 or manifest.get('nodes')!=[0,1,2,3,4]:
        raise RuntimeError('expected frozen q=2 K5 5x32 Euclidean packet')

    cols={};column_gate_pass=True;parameter_consistency=True
    for v in range(5):
        for i in range(32):
            p=packet_dir/'columns'/f'E_node{v}_input{i:02d}.json'
            if not p.exists(): raise RuntimeError(f'missing column {p}')
            d=json.loads(p.read_text(encoding='utf-8'))
            column_gate_pass &= bool(d.get('passed',False))
            parameter_consistency &= d.get('family')=='E' and int(d.get('node',-1))==v and int(d.get('input_index',-1))==i
            cols[v,i]=decode(d['complete_gauss_outgoing_column']['state'])

    node_grams=[];mu1=np.zeros((32,32),complex)
    for v in range(5):
        G=np.zeros((32,32),complex)
        for i in range(32):
            for j in range(32): G[i,j]=inner(cols[v,i],cols[v,j])
        node_grams.append(G);mu1+=G

    herm_err=float(np.linalg.norm(mu1-mu1.conj().T))
    H=.5*(mu1+mu1.conj().T)
    ev=np.linalg.eigvalsh(H);scale=max(float(np.max(np.abs(ev))),1.0);rank_tol=scale*1e-10
    rank=int(np.sum(ev>rank_tol));nullity=32-rank
    trace=float(np.trace(H).real);frob=float(np.linalg.norm(H));emin=float(ev[0]);emax=float(ev[-1])
    node_traces=[float(np.trace(G).real) for G in node_grams]
    node_rel_spread=(max(node_traces)-min(node_traces))/max(abs(sum(node_traces)/5),1e-300)

    ref=manifest['M_EE']
    regression={
        'rank_matches_manifest':rank==int(ref['rank']),
        'nullity_matches_manifest':nullity==int(ref['nullity']),
        'eigenvalue_min_matches_manifest':abs(emin-float(ref['eigenvalue_min']))<2e-10,
        'eigenvalue_max_matches_manifest':abs(emax-float(ref['eigenvalue_max']))<2e-10,
        'trace_matches_manifest':abs(trace-float(ref['trace']))<2e-10,
        'frobenius_matches_manifest':abs(frob-float(ref['frobenius_norm']))<2e-10,
    }
    hard={
        'manifest_domain_complete_is_false':manifest.get('domain_complete') is False,
        'exact_5x32_column_coverage':len(cols)==160,
        'all_column_gates_passed':bool(column_gate_pass),
        'column_parameters_consistent':bool(parameter_consistency),
        'mu1_hermitian':herm_err<2e-10,
        'mu1_positive_semidefinite':emin>-2e-10*scale,
        'reconstructed_mu1_matches_manifest_invariants':all(regression.values()),
    }
    mu1_hash=canonical_matrix_hash(H)
    return {
        'schema':SCHEMA_OUT,
        'passed':bool(all(hard.values())),
        'status':'actual E-sector first projected-heat / block-Krylov moment',
        'source_packet_schema':manifest['schema'],
        'source_packet_sha256':manifest.get('packet_sha256'),
        'family':'E',
        'domain_label':manifest.get('domain_label'),
        'boundary_dimension':32,
        'boundary_injection':'V0 = orthonormal q=2 all-j=1/2 K5 logical boundary basis',
        'mu0':{'definition':'V0^dagger V0','dimension':32,'identity':True},
        'mu1':{
            'definition':'V0^dagger M_E V0 = sum_v (H_E,v V0)^dagger(H_E,v V0) = M_EE',
            'rank':rank,'nullity':nullity,'rank_tolerance':rank_tol,
            'eigenvalue_min':emin,'eigenvalue_max':emax,'trace':trace,'frobenius_norm':frob,
            'hermiticity_error':herm_err,'matrix_hash':mu1_hash,
            'real':H.real.tolist(),'imag':H.imag.tolist(),
        },
        'projected_heat_short_tau':{
            'object':'K_E(tau)=V0^dagger exp(-tau M_E) V0',
            'K_at_tau0':'I_32',
            'first_derivative_at_tau0':'-mu1',
            'trace_at_tau0':32.0,
            'trace_first_derivative_at_tau0':-trace,
            'directional_first_derivative_range_for_unit_boundary_vectors':[-emax,-emin],
            'expansion':'K_E(tau)=I_32-tau*mu1+O(tau^2)',
            'compressed_exponential_used':False,
            'invalid_substitution_exp_minus_tau_mu1_used':False,
        },
        'block_krylov_status':{
            'actual_moments_available':[0,1],
            'higher_moments_available':False,
            'mu2_requires_enlarged_master_action_beyond_one_hit_packet':True,
            'physical_projector_limit_available':False,
        },
        'production_frontier':{
            'euclidean_first_moment_measured':True,
            'lorentzian_cross_and_quadratic_master_blocks_complete':False,
            'quantum_HH_Dtarget_habitat_closure_complete':False,
            'full_constraint_family_complete':False,
            'physical_projector_emitted':False,
            'source_dressed_W_BQG_emitted':False,
        },
        'node_trace_relative_spread':node_rel_spread,
        'manifest_regression_checks':regression,
        'hard_integrity_checks':hard,
        'claim_boundary':(
            'This is the first actual theory-specific projected-history datum from the real BQG Euclidean constraint packet: '
            'mu1 and dK_E/dtau at tau=0. It is not exp(-tau M_EE), not the tau->infinity projector, not P_phys, and not a physical-time propagator. '
            'The full BQG history still requires Lorentzian master data, quantum HH<->Dtarget habitat certification, and higher enlarged-space Krylov moments.'
        ),
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--packet-dir',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    out=run(a.packet_dir);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='mu1'},indent=2));print(json.dumps({'mu1':{k:v for k,v in out['mu1'].items() if k not in ('real','imag')}},indent=2));return 0 if out['passed'] else 1
if __name__=='__main__': raise SystemExit(main())
