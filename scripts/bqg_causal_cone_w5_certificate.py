#!/usr/bin/env python3
"""Exact rational certificate for the published finite W5 Gram payload.

The matrix entries below are copied verbatim from
research/node1_fullstate_gate_2026_10_31/V04_CROSS_BRANCH_REAL_EUCLIDEAN_W5_LOG_EVIDENCE.json
on branch research/bqg-v26-multiprecision-charge-20261109.

This certifies the Hermitian decimal matrix in that artifact. It does not
retroactively certify upstream SU(2) floating-point amplitude generation.
"""
from fractions import Fraction as F
from itertools import permutations
from math import sqrt, asinh

RE = [
 ["34.040571789989606","0.39404535517790373","1.3150915174897833","8.436719205195953e-16","1.345420211432997"],
 ["0.39404535517790373","33.26065884643677","0.3083892701385429","0.5387850728686964","1.400247105765029e-16"],
 ["1.3150915174897833","0.3083892701385429","33.66542393826306","0.23072933451788813","0.8051819063978833"],
 ["8.436719205195953e-16","0.5387850728686964","0.23072933451788813","33.49945906751966","0.12898161909632236"],
 ["1.345420211432997","1.400247105765029e-16","0.8051819063978833","0.12898161909632236","35.41445901968747"],
]
IM = [
 ["0","8.892800195574525e-16","-2.5858964809436897e-17","1.0087846118837086e-16","4.1601438882811575e-16"],
 ["-8.892800195574525e-16","0","-1.2826380347628653e-16","-9.02521422332931e-16","7.416687809231958e-17"],
 ["2.5858964809436897e-17","1.2826380347628653e-16","0","1.265563501359149e-15","-1.1956371511799464e-15"],
 ["-1.0087846118837086e-16","9.02521422332931e-16","-1.265563501359149e-15","0","1.089957273706888e-15"],
 ["-4.1601438882811575e-16","-7.416687809231958e-17","1.1956371511799464e-15","-1.089957273706888e-15","0"],
]

def zadd(a,b): return (a[0]+b[0], a[1]+b[1])
def zmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def zconj(a): return (a[0],-a[1])

def permutation_sign(p):
    inv=sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
    return -1 if inv%2 else 1

def determinant_exact(A):
    n=len(A)
    total=(F(0),F(0))
    for p in permutations(range(n)):
        term=(F(permutation_sign(p)),F(0))
        for i,j in enumerate(p):
            term=zmul(term,A[i][j])
        total=zadd(total,term)
    return total

def main():
    n=5
    G=[[(F(RE[i][j]),F(IM[i][j])) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            assert G[i][j]==zconj(G[j][i]), (i,j,"matrix is not exactly Hermitian")
    q=F("32.41356")
    shifted=[[G[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        shifted[i][i]=(shifted[i][i][0]-q,shifted[i][i][1])
    minors=[]
    for k in range(1,n+1):
        d=determinant_exact([row[:k] for row in shifted[:k]])
        assert d[1]==0 and d[0]>0, (k,d)
        minors.append(d[0])
    print("W5 exact Hermitian check: PASS")
    print("threshold q =",q,"=",float(q))
    print("Sylvester leading principal minors of G-qI:")
    for k,d in enumerate(minors,1):
        print(f"  D{k} = {d} = {float(d):.16g} > 0")
    print("CERTIFIED: G_EE - 32.41356 I is positive definite")
    print("therefore lambda_min(G_EE) > 32.41356")
    sig=sqrt(float(q))
    print(f"sigma_min(A) > sqrt(32.41356) = {sig:.15g}")
    print(f"for H=-(2/3)A+Q: ||Q||_op > {2*sig/3:.15g} is necessary for a zero")
    print(f"for H=A+lambda B: |lambda| ||B||_op > {sig:.15g} is necessary for a zero")
    print("NOTE: B is not present in the W5 evidence; this is a necessary threshold, not a measured Lorentzian cancellation.")
    print("\nReduced free-fermion propagation cone:")
    print("m       Delta_sp      mu_c             xi               v_front")
    for m in (0.0,0.2,0.5,1.0,2.0,5.0):
        gap=2*abs(m)
        if m==0: muc=0.0; xi=float("inf")
        else: muc=asinh(gap/3.0); xi=1.0/muc
        v=sqrt(gap*gap+9.0)-gap
        xis="inf" if xi==float("inf") else f"{xi:.12f}"
        print(f"{m:<7g}{gap:<14.8g}{muc:<17.12f}{xis:<17}{v:.12f}")
    print("Exact identities: mu_c=asinh(Delta_sp/3), xi=1/mu_c,")
    print("v_front=sqrt(Delta_sp^2+9)-Delta_sp=3 exp(-1/xi).")
    print("This cone is proved for the reduced free-fermion model, not full physical BQG.")

if __name__=="__main__":
    main()
