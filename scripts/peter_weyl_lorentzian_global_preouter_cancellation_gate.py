#!/usr/bin/env python3
"""Measure the globally signed Lorentzian middle sum before the final C(K).

For a certified all-middle packet define, for each auxiliary (i,j),

    S_ij = sum_a sum_{b,c,k} epsilon_{abc} Xi_bc^{ijk},
    Xi_bc^{ijk}=C_b(K)_{jk} C_c(V)_{ki}|psi>.

No outer C(K) is evaluated here.  The gate measures whether the signed middle
sum cancels before the last operator.  This matters computationally because the
``delta_ij K`` part of C_a(K)_ij is independent of the outer edge a; if S_ii
cancels, the complete direct-K contribution cancels globally before any
edge-dependent h_a K h_a^-1 term is considered.

Cancellation is a measured finite-packet result, not assumed for other logical
inputs or refinements. PASS checks packet integrity; the cancellation
certificate is reported separately.
"""
from __future__ import annotations

import argparse, json, math
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_lorentzian_middle_prefix_gate as MID
import peter_weyl_lorentzian_pruned_prefix_worker as PLAN

ABS_CERT_TOL=1.0e-13
REL_CERT_TOL=1.0e-12


def decode(rows):
    out={}
    for r in rows:
        key=(tuple(int(x) for x in r['spins']),tuple(int(x) for x in r['Kother']),
             int(r['J2']),int(r['M2']),int(r['K12']),int(r['K34']))
        out[key]=out.get(key,0j)+complex(float(r['amp'][0]),float(r['amp'][1]))
    return out


def add_exact(dst,src,scale=1):
    for k,z in src.items(): dst[k]=dst.get(k,0j)+complex(scale)*z
    for k in [k for k,z in dst.items() if z==0j]: del dst[k]


def norm(s): return math.sqrt(sum(abs(z)**2 for z in s.values()))
def scalar(s): return {k:z for k,z in s.items() if int(k[2]) in (0,2)}


def run(prefix_dir:Path):
    summary=json.loads((prefix_dir/'middle_prefix_summary.json').read_text(encoding='utf-8'))
    if summary.get('schema')!='BQG_LORENTZIAN_ALL_MIDDLE_PREFIXES_V1' or summary.get('passed') is not True:
        raise RuntimeError('invalid all-middle summary')
    source=int(summary['source_node']);input_index=int(summary['input_logical_basis_index'])
    neighbors=tuple(FULL.RAW.PW.NEIG[source]);pairs=MID.ordered_pairs(source)
    if len(neighbors)!=4 or len(pairs)!=12: raise RuntimeError('expected four-valent source / 12 ordered prefixes')
    prefixes=[]
    for idx in range(12):
        p=json.loads((prefix_dir/f'prefix_{idx}.json').read_text(encoding='utf-8'))
        if p.get('passed') is not True or int(p['ordered_pair']['pair_index'])!=idx: raise RuntimeError(f'invalid prefix {idx}')
        if int(p['source_node'])!=source or int(p['input_logical_basis_index'])!=input_index: raise RuntimeError(f'prefix {idx} provenance mismatch')
        if (int(p['ordered_pair']['b']),int(p['ordered_pair']['c']))!=pairs[idx]: raise RuntimeError(f'prefix {idx} pair ordering mismatch')
        prefixes.append(p)

    groups={};outer_rows=[];all_triples=[]
    for a in neighbors:
        for i in range(2):
            for j in range(2):
                g={};contributors=[]
                for idx,p in enumerate(prefixes):
                    b=int(p['ordered_pair']['b']);c=int(p['ordered_pair']['c'])
                    rows=[r for r in PLAN.epsilon_outer_terms(source,b,c) if int(r['a'])==a]
                    if not rows: continue
                    if len(rows)!=1: raise RuntimeError('duplicate fixed-a epsilon row')
                    row=rows[0];sign=int(row['sign']);all_triples.append(tuple(int(x) for x in row['ordered_edges'])) if (i,j)==(0,0) else None
                    pm={tuple(int(x) for x in q['indices']):q for q in p['paths']}
                    before=norm(g)
                    for k in range(2): add_exact(g,scalar(decode(pm[(i,j,k)]['middle_state'])),scale=sign)
                    contributors.append({'prefix_index':idx,'sign':sign,'group_norm_before':before,'group_norm_after':norm(g)})
                groups[(a,i,j)]=g
                outer_rows.append({'outer_a':a,'indices_ij':[i,j],'support':len(g),'norm':norm(g),'contributors':contributors})

    global_rows=[];all_cert=True
    for i in range(2):
        for j in range(2):
            s={}
            quadrature=0.0
            for a in neighbors:
                add_exact(s,groups[(a,i,j)])
                quadrature+=norm(groups[(a,i,j)])**2
            n=norm(s);den=math.sqrt(quadrature);mx=max((abs(z) for z in s.values()),default=0.0);rel=n/max(den,1e-300)
            cert=mx<ABS_CERT_TOL and rel<REL_CERT_TOL
            all_cert &= cert
            global_rows.append({'indices_ij':[i,j],'residual_support':len(s),'residual_norm':n,'max_abs_amplitude':mx,
                                'quadrature_reference_norm':den,'relative_residual':rel,'cancellation_certified':bool(cert)})

    hard={
        'summary_passed':summary.get('passed') is True,
        'all_12_prefixes_loaded':len(prefixes)==12,
        'all_24_epsilon_triples_accounted_once':len(all_triples)==24 and len(set(all_triples))==24,
        'all_16_fixed_aij_groups_constructed':len(groups)==16,
        'no_pre_outer_tolerance_pruning':True,
    }
    diagonal_cert=all(r['cancellation_certified'] for r in global_rows if r['indices_ij'][0]==r['indices_ij'][1])
    return {
        'schema':'BQG_LORENTZIAN_GLOBAL_PREOUTER_CANCELLATION_V1',
        'passed':bool(all(hard.values())),
        'science_status':'GLOBAL_PREOUTER_SUM_CANCELS_WITHIN_CERT_TOL' if all_cert else 'GLOBAL_PREOUTER_SUM_NONZERO',
        'source_node':source,'input_logical_basis_index':input_index,'Jmax':float(summary['Jmax']),
        'absolute_certificate_tolerance':ABS_CERT_TOL,'relative_certificate_tolerance':REL_CERT_TOL,
        'fixed_outer_groups':outer_rows,'global_signed_middle_sums':global_rows,
        'all_aux_channels_cancellation_certified':bool(all_cert),
        'direct_K_global_cancellation_certified':bool(diagonal_cert),
        'hard_integrity_checks':hard,
        'claim_boundary':'Measured pre-outer cancellation on one declared finite source/input packet only. If direct_K_global_cancellation_certified=true, only the outer-edge-independent delta_ij K contribution is certified to cancel for this packet; edge-dependent h_a K h_a^-1 actions, the full H_L column, HDA, P_phys and cosmology remain separate computations.'
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--prefix-dir',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    out=run(a.prefix_dir);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in out.items() if k!='fixed_outer_groups'},indent=2));return 0 if out['passed'] else 1
if __name__=='__main__': raise SystemExit(main())
