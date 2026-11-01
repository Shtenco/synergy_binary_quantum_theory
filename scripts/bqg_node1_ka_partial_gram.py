#!/usr/bin/env python3
"""Real K5 node1 KB→KA partial-wave coherent Gram audit.

Consumes ONLY original saved GitHub Actions .npz/.json checkpoints; verifies
each SHA256, same frozen operator code, source spin columns and KB parent.
This is NOT the complete Lorentzian operator or a physical HDA/TT result.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np

STUDIED=(0,1,2,3,4,5,12,13,16,17,32,33,58)

def checksum(path):
 h=hashlib.sha256()
 with open(path,"rb") as f:
  for blob in iter(lambda:f.read(1<<20),b""):h.update(blob)
 return h.hexdigest()

def vector(root,index):
 stem=root/f"ka_000_part{index:03d}_of120"
 m=json.loads(stem.with_suffix(".json").read_text(encoding="utf-8"))
 assert m["channel"]==[0,0,0] and m["part_index"]==index and m["parts"]==120
 assert m["rows"]==1 and m["sha256"]==checksum(stem.with_suffix(".npz"))
 with np.load(stem.with_suffix(".npz"),allow_pickle=False) as d:
  assert tuple(d.files)==("spins","Kother","charges","amp")
  assert d["spins"].shape==(m["support"],10)
  assert d["Kother"].shape==(m["support"],5)
  assert d["charges"].shape==(m["support"],4)
  keys=list(zip(map(tuple,d["spins"]),map(tuple,d["Kother"]),map(tuple,d["charges"])))
  assert len(set(keys))==len(keys)
  vals=np.asarray(d["amp"],dtype=np.complex128)
  assert np.isfinite(vals.real).all() and np.isfinite(vals.imag).all()
  assert abs(float(np.linalg.norm(vals))-m["norm"])<1e-12
  return dict(zip(keys,vals)),m

def calc(root,study=STUDIED):
 kb=json.loads((root/"kb_000.json").read_text(encoding="utf-8"))
 assert kb["stage"]=="kb" and kb["sha256"]==checksum(root/"kb_000.npz")
 with np.load(root/"kb_000.npz",allow_pickle=False) as b:
  all_spins=b["spins"];all_amplitudes=b["amp"];assert len(all_spins)==len(all_amplitudes)==120
 basis=[];meta=[]
 for idx in study:
  vec,mm=vector(root,idx)
  assert mm["freeze"]==kb["freeze"],"mixed source/operator checkpoint"
  assert mm["kb_sha256"]==kb["sha256"],"mixed KB parent"
  basis.append(vec);meta.append(mm)
 union=set().union(*(x.keys() for x in basis))
 gram=np.zeros((len(study),len(study)),complex)
 for i,one in enumerate(basis):
  for j in range(i,len(basis)):
   other=basis[j]
   val=sum(np.conjugate(one[k])*other[k] for k in one.keys()&other.keys())
   gram[i,j]=val;gram[j,i]=np.conjugate(val)
 assert np.linalg.norm(gram-gram.conj().T)<1e-12
 eig=np.linalg.eigvalsh(gram)
 errtol=1e-11
 assert eig.min()>-errtol,("negative Gram",eig.min())
 total={}
 for item in basis:
  for key,amp in item.items():total[key]=total.get(key,0)+amp
 coherent_norm=float(np.sqrt(sum(abs(v)**2 for v in total.values())))
 trace=float(np.trace(gram).real)
 coherence_from_gram=float(np.sqrt(np.ones(len(study))@gram@np.ones(len(study))).real)
 assert abs(coherent_norm-coherence_from_gram)<1e-10
 overlaps=[]
 for i in range(len(study)):
  for j in range(i+1,len(study)):
   v=gram[i,j]/np.sqrt(gram[i,i]*gram[j,j])
   overlaps.append((float(abs(v)),study[i],study[j]))
 frac=float(sum(abs(all_amplitudes[i])**2 for i in study)/sum(abs(all_amplitudes)**2))
 result={
  "status":"PARTIAL_TRUE_PETER_WEYL_RESULTS_ONLY",
  "physical_HDA":"OPEN","physical_TT_pole":"OPEN","full_24_term_HL":"OPEN",
  "source":"node1, triple (0,2,3), auxiliary 000, fixed source column0",
  "measured_KB_input_keys":len(all_spins),"completed_KA_input_keys":len(study),
  "completed_key_indices":list(study),
  "distinct_complete_covariant_outgoing_keys":len(union),
  "coherent_output_norm":coherent_norm,
  "quadrature_output_norm":float(np.sqrt(trace)),
  "squared_norm_interference":float(coherent_norm**2-trace),
  "Gram_spectrum":eig.tolist(),
  "Gram_rank_numpy":int(np.linalg.matrix_rank(gram)),
  "largest_absolute_normalized_overlap":max(overlaps) if overlaps else None,
  "input_KB_squared_weight_covered":frac,
  "input_KB_squared_weight_unprocessed":1-frac,
  "CAUTION":"Unprocessed input norm squared fraction is NOT an outgoing error bound without ||C_a(K)|| operator bound.",
  "primitive_basis_leakage_status":"NOT CLEARED — frozen max primitive CK leakage=1.0",
  "source_kb_sha256":kb["sha256"]
 }
 return result

if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--output",type=Path)
 a=p.parse_args();result=calc(a.root)
 message=json.dumps(result,indent=2,ensure_ascii=False)
 print(message,flush=True)
 if a.output:a.output.write_text(message+"\n",encoding="utf-8")
