#!/usr/bin/env python3
"""K5 SU2 initial-source vertex swap + archived source0 Walsh conditional node1 predictor.
Use original PW V5 tensor; NOT an independent node1 Lorentzian contraction.
"""
import itertools,json,math,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import k5_peter_weyl_safe_hda_column as PW
def kr(*ops):
    z=np.array([[1]],complex)
    for a in ops:z=np.kron(z,a)
    return z
I=np.eye(2);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1,-1]);X=np.array([[0,1],[1,0]])
spins=(1,)*4
inter=[PW.intertwiner_tensor_cached(spins,k) for k in (0,2)]
local=np.array([[np.vdot(inter[i],inter[j].transpose(1,0,2,3)) for j in range(2)] for i in range(2)])
assert np.linalg.norm(local-np.diag([-1,1]))<1e-12
bits=list(itertools.product(range(2),repeat=5))
U=np.zeros((32,32),complex)
for j,ks in enumerate(bits):
    target=(ks[1],ks[0],ks[2],ks[3],ks[4])
    U[bits.index(target),j]=np.prod([-1 if ks[v]==0 else 1 for v in (2,3,4)])
v=PW.v5_tensor()
assert np.linalg.norm(U@v+v)/np.linalg.norm(v)<1e-12
assert np.linalg.norm(U@U.conj().T-np.eye(32))<1e-12
# Frozen genuine archived raw Pauli-imag coefficients from the full 24-term source0
# verification_results/PETER_WEYL_LORENTZIAN_ENVTRACE_WALSH_NODE012.json
evidence=Path("verification_results/PETER_WEYL_LORENTZIAN_ENVTRACE_WALSH_NODE012.json")
if not evidence.is_file():raise SystemExit("BLOCKED: missing source0 archived Walsh evidence")
raw=json.loads(evidence.read_text(encoding="utf-8"))
assert raw.get("passed") and raw.get("source_node")==0
keys=["YI1I2","YZ1I2","YI1Z2","YZ1Z2"]
a,b,c,d=[-32/9*raw["raw_pauli"][key][1] for key in keys]
H0=a*kr(Y,I,I)+b*kr(Y,Z,I)+c*kr(Y,I,Z)+d*kr(Y,Z,Z)
H1=a*kr(I,Y,I)+b*kr(Z,Y,I)+c*kr(I,Y,Z)+d*kr(Z,Y,Z)
V=np.zeros((8,8),complex);ib=list(itertools.product(range(2),repeat=3))
for j,ks in enumerate(ib):
    target=(ks[1],ks[0],ks[2])
    V[ib.index(target),j]=-1 if ks[2]==0 else 1
assert np.linalg.norm(V@H0@V.conj().T-H1)<1e-12
C=-1j*(H0@H1-H1@H0)
u2=2*(a*b+c*d);u3=2*(a*d+b*c)
Cdecomp=u2*(kr(X,Y,I)-kr(Y,X,I))+u3*(kr(X,Y,Z)-kr(Y,X,Z))
assert np.linalg.norm(C-Cdecomp)<1e-12
print(json.dumps({"kinematic_swap_control":"PASS","node0_raw_evidence":"ARCHIVED_MEASURED",
  "node1":"SYMMETRY_PREDICTION_ONLY_NOT_DIRECT_MICROSCOPIC",
  "V5_swap_eigen_minus1_residual":float(np.linalg.norm(U@v+v)/np.linalg.norm(v)),
  "source0_signed_pauli_coeff":[a,b,c,d],
  "conditional_two_node_C_coeff":u2,"conditional_three_node_C_coeff":u3,
  "conditional_C_spectral_norm":float(np.linalg.norm(C,ord=2)),
  "conditional_C_frobenius_norm":float(np.linalg.norm(C,"fro")),
  "conditional_C_rank":int(np.linalg.matrix_rank(C)),
  "physical_full_HDA":"BLOCKED: no off-code Q returns, no independent D"},indent=2))
