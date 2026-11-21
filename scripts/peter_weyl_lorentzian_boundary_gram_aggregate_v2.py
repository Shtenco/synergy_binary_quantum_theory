#!/usr/bin/env python3
"""Order-invariant production wrapper for the Lorentzian boundary Gram aggregate.

The V1 aggregate correctly verifies that the 32 logical input indices form the
set {0,...,31}, but its final integrity flag also compared the *incoming file
order* to range(32).  A glob is not required to be numerically ordered, so that
check can create a false fail-closed result for a complete packet.

This wrapper canonicalizes input artifacts by their declared
`input_logical_basis_index` before calling the frozen V1 physics aggregate.  It
changes no amplitudes, operators, thresholds, Gram construction, parity test or
S4 test.  Duplicate/missing indices remain hard errors in V1.
"""
from __future__ import annotations

import argparse, glob, json
from pathlib import Path

import peter_weyl_lorentzian_boundary_gram_aggregate as V1


def canonicalize(paths):
    keyed=[]
    for p in paths:
        p=Path(p)
        d=json.loads(p.read_text(encoding='utf-8'))
        if 'input_logical_basis_index' not in d:
            raise RuntimeError(f'{p}: missing input_logical_basis_index')
        keyed.append((int(d['input_logical_basis_index']),p))
    idx=[i for i,_ in keyed]
    if len(idx)!=32 or set(idx)!=set(range(32)) or len(set(idx))!=32:
        raise RuntimeError(f'expected exactly one artifact for each input 0..31, got {sorted(idx)}')
    keyed.sort(key=lambda x:x[0])
    return [p for _,p in keyed]


def run(paths,e_packet_dir=None):
    ordered=canonicalize(paths)
    out=V1.run(ordered,e_packet_dir)
    out['schema']='BQG_LORENTZIAN_BOUNDARY_GRAM_V2'
    out['input_order_canonicalization']={
        'applied':True,
        'canonical_order':list(range(32)),
        'physics_changed':False,
        'reason':'prevent filesystem/glob ordering from affecting exact_32_input_coverage',
    }
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--column',type=Path,action='append')
    ap.add_argument('--column-glob')
    ap.add_argument('--euclidean-packet-dir',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    paths=list(a.column or [])
    if a.column_glob:
        paths.extend(Path(x) for x in glob.glob(a.column_glob))
    out=run(paths,a.euclidean_packet_dir)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,indent=2))
    return 0 if out['passed'] else 1

if __name__=='__main__':
    raise SystemExit(main())
