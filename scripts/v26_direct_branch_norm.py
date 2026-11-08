#!/usr/bin/env python3
"""V26 independent postprocessing of original tensor-product SU(2) branches.
This improves accumulation, NOT the accuracy of original binary64 SU(2) values.
Failure is explicit when a product basis is too large. No physical HDA claims.
"""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction
from math import isfinite, prod
import numpy as np
import mpmath as mp

class CertificationUnavailable(RuntimeError): pass

def _parse(branches, vertices):
    if not vertices or len(set(vertices))!=len(vertices):
        raise ValueError("INVALID_VERTICES")
    groups=defaultdict(list)
    for spin,tensors,c0 in branches:
        c=complex(c0)
        if not isfinite(c.real) or not isfinite(c.imag):
            raise ValueError("NONFINITE_COEFFICIENT")
        arrays=[np.asarray(tensors[v],dtype=np.complex128).ravel() for v in vertices]
        if any(len(a)==0 or not np.isfinite(a).all() for a in arrays):
            raise ValueError("INVALID_LOCAL_TENSOR")
        groups[tuple(spin)].append((arrays,c))
    return groups

def _mp_float(x):
    n,d=float(x).as_integer_ratio()
    return mp.mpf(n)/d

def _mp(z):
    z=complex(z)
    return mp.mpc(_mp_float(z.real),_mp_float(z.imag))

def squared_norm(branches,vertices,*,digits=80,max_group_dim=4096,max_work=2000000):
    if digits<35: raise ValueError("PRECISION_TOO_LOW")
    grouped=_parse(branches,vertices)
    work=0
    for group in grouped.values():
        dims=tuple(len(x) for x in group[0][0])
        if any(tuple(map(len,arrs))!=dims for arrs,c in group):
            raise ValueError("INCONSISTENT_LOCAL_PRODUCT_BASIS")
        dimension=prod(dims)
        if dimension>max_group_dim:
            raise CertificationUnavailable("BASIS_DIMENSION_LIMIT")
        work+=dimension*len(group)
    if work>max_work:
        raise CertificationUnavailable("WORK_BUDGET_EXCEEDED")
    with mp.workdps(digits):
        contributions=[]
        sector_info=[]
        for spin,group in sorted(grouped.items(),key=lambda x:repr(x[0])):
            dim=prod(map(len,group[0][0]))
            columns=[[] for _ in range(dim)]
            for arrays,c in group:
                values=[_mp(c)]
                for a in arrays:
                    values=[v*_mp(z) for v in values for z in a]
                for j,z in enumerate(values): columns[j].append(z)
            final=[mp.mpc(mp.fsum(z.real for z in column),
                          mp.fsum(z.imag for z in column)) for column in columns]
            sn=mp.fsum(abs(z)**2 for z in final)
            contributions.append(sn)
            sector_info.append({'spin':repr(spin),'branches':len(group),'basis_dim':dim,
                                'squared_norm':mp.nstr(sn,digits-5)})
        squared=mp.fsum(contributions)
        return {'status':'BINARY64_INPUT_POSTPROCESSING_ONLY',
                'squared_norm':mp.nstr(squared,digits-5),
                'norm':mp.nstr(mp.sqrt(squared),digits-5),
                'number_of_branches':sum(map(len,grouped.values())),
                'number_of_spin_sectors':len(grouped),
                'digits':digits,'sectors':sector_info,
                'original_SU2_rounding_certified':False,
                'coherent_charge_projection_proved':False,
                'full_HDA_proved':False}

def exact_small_squared_norm(branches,vertices,*,max_terms=1000):
    """Exact rational result for previously rounded binary64 input only."""
    groups=_parse(branches,vertices)
    def cfrac(z):
        z=complex(z)
        return (Fraction(*z.real.as_integer_ratio()),Fraction(*z.imag.as_integer_ratio()))
    def mul(x,y):return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
    def add(x,y):return x[0]+y[0],x[1]+y[1]
    work=0;total=Fraction(0)
    for spin,group in groups.items():
        dims=tuple(map(len,group[0][0]))
        if any(tuple(map(len,arrs))!=dims for arrs,c in group):
            raise ValueError("INCONSISTENT_LOCAL_PRODUCT_BASIS")
        dim=prod(dims);work+=dim*len(group)
        if work>max_terms:raise CertificationUnavailable("EXACT_RATIONAL_WORK_LIMIT")
        coords=[(Fraction(0),Fraction(0)) for _ in range(dim)]
        for arrays,c in group:
            values=[cfrac(c)]
            for a in arrays: values=[mul(v,cfrac(z)) for v in values for z in a]
            coords=[add(a,b) for a,b in zip(coords,values)]
        total+=sum((a*a+b*b for a,b in coords),Fraction(0))
    return total
