#!/usr/bin/env python3
"""Pure V21 Lorentzian workflow dependency/census checks; no SU(2) physics."""
from itertools import product


def operations(targets,mode,aux):
    a,b,c=targets;i,j,k=aux
    if mode=='forward':return [('CV',c,k,i),('CK',b,j,k),('CK',a,i,j)]
    if mode=='adjoint':return [('CK',a,j,i),('CK',b,k,j),('CV',c,i,k)]
    raise ValueError('INVALID_MODE')


def all_32_specs(abc):
    a,b,c=abc
    for mode,order in product(('forward','adjoint'),('p','q')):
        t=(a,b,c) if order=='p' else (b,a,c)
        for aux in product((0,1),repeat=3):
            yield {'mode':mode,'order':order,'aux':list(aux),'targets':t,
                   'operations':operations(t,mode,aux)}


def first_leg_census(abc):
    specs=list(all_32_specs(abc))
    first_ops={tuple(row['operations'][0]) for row in specs}
    pairs={(row['mode'],row['order'],tuple(row['aux'])) for row in specs}
    if len(pairs)!=32 or len(first_ops)!=12:
        raise ValueError('INVALID_32_PATH_CENSUS')
    return len(first_ops),len(pairs)


def predecessor_ok(prev,cur,stage,first_operation,second_operation=None):
    if prev.get('stage')!=stage-1 or prev.get('passed') is not True:
        raise ValueError('STAGE_NOT_PASSED')
    if stage not in (2,3):raise ValueError('INVALID_STAGE')
    for key in ('source_commit','source_node','omit','cycle','mode','Jmax2',
                'volume_backend','volume_backend_sha256','input_state_sha256'):
        if prev.get(key)!=cur.get(key):raise ValueError('MISMATCH_'+key)
    if stage==2:
        # First operator can be shared only when operator signatures match;
        # p/q and unused auxiliary indices need not match.
        if tuple(prev.get('operator_signature',()))!=tuple(first_operation):
            raise ValueError('MISMATCH_FIRST_OPERATOR_SIGNATURE')
    else:
        for key in ('order','aux'):
            if prev.get(key)!=cur.get(key):raise ValueError('MISMATCH_'+key)
        if not prev.get('post_second_scalar_filter'):
            raise ValueError('STAGE2_FILTER_NOT_APPLIED')
        if tuple(prev.get('operator_signature',()))!=tuple(second_operation):
            raise ValueError('MISMATCH_SECOND_OPERATOR_SIGNATURE')
    return True
