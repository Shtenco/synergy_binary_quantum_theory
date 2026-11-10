"""Replay saved finite Euclidean evidence; does not regenerate SU(2) amplitudes.

Only load the trusted pickle files committed alongside this script.
"""
import json,sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'proof_tools'))
from safe_pickle import safe_load_bytes
import numpy as np
from scipy.sparse import load_npz
p=Path(__file__).resolve().parent
A=load_npz(p/'fine_H0_full_images.npz').toarray();rows=[(tuple(a),tuple(b)) for a,b in json.loads((p/'fine_H0_full_keys.json').read_text())]
h0={k:complex(A[i,0]) for i,k in enumerate(rows) if abs(A[i,0])>1e-12}
h1={(tuple(a),tuple(b)):complex(re,im) for a,b,re,im in json.loads((p/'node1_00.json').read_text())}
outs=[]
for name,inp in [('H0H1',h1),('H1H0',h0)]:
 out={}
 for j,(k,c) in enumerate(sorted(inp.items())):
  im=safe_load_bytes((p/f'{name}_{j:03d}.pkl').read_bytes())
  for key,z in im.items():out[key]=out.get(key,0j)+c*z
 out={k:z for k,z in out.items() if abs(z)>1e-10}
 saved=safe_load_bytes((p/f'{name}.pkl').read_bytes())
 assert max(abs(out.get(k,0j)-saved.get(k,0j)) for k in set(out)|set(saved))<1e-12
 outs.append(out)
c={k:outs[0].get(k,0j)-outs[1].get(k,0j) for k in set(outs[0])|set(outs[1])}
norm=float(np.sqrt(sum(abs(z)**2 for z in c.values())))
assert abs(norm-json.loads((p/'results.json').read_text())['commutator_norm'])<1e-12
print(json.dumps({'saved_evidence_replay_passed':True,'commutator_norm':norm,'fresh_operator_generation':False}))
