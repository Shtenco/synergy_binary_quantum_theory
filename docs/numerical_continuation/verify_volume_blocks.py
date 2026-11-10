"""Regenerate exact local SU(2) volume blocks; no physical-limit claim."""
import json,itertools
from pathlib import Path
import sympy as sp
from sympy.physics.wigner import clebsch_gordan


def require(condition, message):
    """Keep exact-volume gate invariants active under python -O."""
    if not condition:
        raise AssertionError(message)


root=Path(__file__).resolve().parent

def labels_for(spins,J):
    return [(a,b) for a in range(abs(spins[0]-spins[1]),spins[0]+spins[1]+1,2)
            for b in range(abs(spins[2]-spins[3]),spins[2]+spins[3]+1,2)
            if abs(a-b)<=J<=a+b and (a+b+J)%2==0]

def spin(j):
    m=[j-i for i in range(int(2*j)+1)];up=sp.zeros(len(m))
    for col,x in enumerate(m):
        if x+1 in m:up[m.index(x+1),col]=sp.sqrt(j*(j+1)-x*(x+1))
    return m,[(up+up.T)/2,(up-up.T)/(2*sp.I),sp.diag(*m)]

rows=[]
for sp2,J2 in [((2,2,2,2),0),((1,2,2,2),1),((3,2,2,2),1)]:
    js=[sp.Rational(x,2) for x in sp2];data=[spin(j) for j in js];mag=[d[0] for d in data];labels=labels_for(sp2,J2);cols=[]
    for aa,bb in labels:
        a,b=sp.Rational(aa,2),sp.Rational(bb,2);col=[]
        for m1,m2,m3,m4 in itertools.product(*mag):
            x=m1+m2;y=m3+m4
            if abs(x)>a or abs(y)>b:col.append(0)
            else:col.append(clebsch_gordan(js[0],js[1],a,m1,m2,x)*clebsch_gordan(js[2],js[3],b,m3,m4,y)*clebsch_gordan(a,b,sp.Rational(J2,2),x,y,sp.Rational(J2,2)))
        cols.append(sp.SparseMatrix(col))
    B=sp.SparseMatrix.hstack(*cols);n=B.rows;Q=sp.zeros(n)
    for omit in range(4):
        legs=[v for v in range(4) if v!=omit]
        for xyz in itertools.permutations(range(3)):
            sign=(-1)**sum(xyz[i]>xyz[j] for i in range(3) for j in range(i+1,3));f=[sp.eye(len(m)) for m in mag]
            for leg,a in zip(legs,xyz):f[leg]=data[leg][1][a]
            Q+=sp.Rational((-1)**(omit+1)*sign,4)*sp.kronecker_product(*f)
    q=sp.simplify(B.conjugate().T*Q*B);require(sp.simplify(B.conjugate().T*B)==sp.eye(len(labels)), 'recoupling basis is not isometric');require(sp.simplify(Q*B-B*q)==sp.zeros(n,len(labels)), 'volume operator does not preserve the selected intertwiner carrier')
    S=sp.simplify(sp.trace(q*q)/2);P=sp.simplify((S*S-sp.trace(q**4)/2)/2);lm=(S-sp.sqrt(S*S-4*P))/2;lp=(S+sp.sqrt(S*S-4*P))/2
    if len(labels)==3:
        require(sp.simplify(q**3-S*q)==sp.zeros(3), 'three-dimensional characteristic identity failed')
        V=q*q/S**sp.Rational(3,4)
    elif len(labels)==4:
        a=(lp**sp.Rational(1,4)-lm**sp.Rational(1,4))/(lp-lm);b=(lp*lm**sp.Rational(1,4)-lm*lp**sp.Rational(1,4))/(lp-lm);V=a*q*q+b*sp.eye(4)
        require(sp.simplify(q**4-S*q*q+P*sp.eye(4))==sp.zeros(4), 'four-dimensional characteristic identity failed')
    elif len(labels)==5:
        a=(lp**sp.Rational(-3,4)-lm**sp.Rational(-3,4))/(lp-lm);b=(lp*lm**sp.Rational(-3,4)-lm*lp**sp.Rational(-3,4))/(lp-lm);V=a*q**4+b*q*q
        require(sp.simplify(q**5-S*q**3+P*q)==sp.zeros(5), 'five-dimensional characteristic identity failed')
    else:raise AssertionError('Unexpected charged multiplicity')
    saved=json.loads((root/('exact_source_volume.json' if J2==0 else 'exact_initial_charged_volumes.json')).read_text())
    if J2: saved=next(row for row in saved if tuple(row['spins_doubled'])==sp2)
    expected=sp.Matrix([[sp.sympify(x) for x in row] for row in saved['Q_block']])
    require(sp.simplify(q-expected)==sp.zeros(len(labels)), 'recomputed local volume block differs from the frozen exact source')
    row={'spins_doubled':sp2,'J2':J2,'dimension':len(labels),'exact_Q_and_charge_invariance_checked':True,'polynomial_checked':True,'physicalization_complete':False}
    rows.append(row)
print(json.dumps({'exact_local_blocks':rows,'scope':'regular-volume candidate, selected three blocks only'},indent=2))
