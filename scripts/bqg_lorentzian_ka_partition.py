#!/usr/bin/env python3
"""Partition JUST the expensive second C(K) action over KB covariant keys.
Exact source/intertwiner operator, with explicit control of numerical summation.
Not a physical Lorentzian HDA or TT proof.

Each job reads the *same SHA-checked KB* and computes a disjoint source-key
partition. Merge hashes, combine all pieces, then scalar Gauss closure.
Because the upstream C_K_component and LP.add use fixed tolerance pruning,
the partitioned sum has different floating-point addition order; do not claim
bitwise equality with the original monolithic contraction.
"""
import argparse,hashlib,itertools,json,math,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import bqg_lorentzian_checkpointed_triple as Base
import peter_weyl_lorentzian_logical_projection_gate as LP
import peter_weyl_lorentzian_gauss_action_gate as LGA
import peter_weyl_lorentzian_ordered_triple_gate as RAW

def cli():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--mode",choices=("part","merge"),required=True)
 p.add_argument("--checkpoint-dir",type=Path,required=True)
 p.add_argument("--parts",type=int,default=8)
 p.add_argument("--part-index",type=int,default=0)
 p.add_argument("--channel",choices=["".join(map(str,a)) for a in itertools.product(range(2),repeat=3)],default="000")
 p.add_argument("--source",type=int,default=1)
 p.add_argument("--column-index",type=int,default=0)
 p.add_argument("--a",type=int,default=0);p.add_argument("--b",type=int,default=2);p.add_argument("--c",type=int,default=3)
 p.add_argument("--coefficient",type=int,default=-1)
 p.add_argument("--json-output",type=Path)
 p.add_argument("--state-output",type=Path)
 a=p.parse_args()
 if a.parts<1 or not 0<=a.part_index<a.parts:p.error("invalid partition index")
 return a
def parse_key(s):return tuple(int(c) for c in s)
def expected_freeze(a):
 if a.source not in Base.PW.VERT or len({a.a,a.b,a.c})!=3 or not all(v in Base.PW.NEIG[a.source] for v in (a.a,a.b,a.c)):
  raise ValueError("partition triple/source does not match the frozen local graph")
 if a.coefficient not in (-1,1):
  raise ValueError("partition epsilon coefficient must be +/-1")
 basis=Base.PW.basis_full_jhalf()
 if not 0<=a.column_index<len(basis):
  raise ValueError("partition column-index is outside the 32D logical basis")
 initial=basis[a.column_index]
 return dict(source=a.source,column_index=a.column_index,input_key=repr(initial),
  input_Ks=list(initial[1]),triple=[a.a,a.b,a.c],epsilon=a.coefficient,
  JMAX2=Base.JMAX2,TOL=Base.TOL,checkpoint_schema=Base.CHECKPOINT_SCHEMA,
  runtime=Base.runtime_fingerprint(),modules=Base.source_hashes())

def key_digest(keys):
 payload="\\n".join(repr(key) for key in keys).encode("utf-8")
 return hashlib.sha256(payload).hexdigest()

def run():
 a=cli();ijk=parse_key(a.channel);checkpoint=a.checkpoint_dir
 expected=expected_freeze(a)
 source,orig=Base.load_result(checkpoint,"kb",ijk,expected)
 freeze=orig["freeze"]
 rows=sorted(source.items(),key=lambda kv:repr(kv[0]))
 part_rows=rows[a.part_index::a.parts]
 subset=dict(part_rows)
 if a.mode=="part":
  LP.JMAX2=Base.JMAX2
  restore,_=LP.install_sine_cached_stack()
  try:
   if subset:
    out,d=RAW.KCOMP.C_K_component(subset,a.source,a.a,ijk[0],ijk[1],Base.JMAX2)
   else:out={};d={}
   diag={("CK_"+key):float(value) for key,value in d.items()
    if key in ("outer_complete_basis_leakage","internal_volume_sector_leakage","complete_charge_basis_leakage")}
   h=Base.maxdiag(orig["leakage_max"],diag)
   stem=checkpoint/f"ka_{a.channel}_part{a.part_index:03d}_of{a.parts:03d}"
   Base.persist_state(stem.with_suffix(".npz"),out)
   meta=dict(stage="ka_part",parts=a.parts,part_index=a.part_index,channel=list(ijk),
      kb_sha256=orig["sha256"],freeze=freeze,rows=len(subset),
      input_keys_sha256=key_digest([key for key,_ in part_rows]),support=len(out),
      norm=Base.norm(out),sha256=Base.digest(stem.with_suffix(".npz")),
      leakage_max=h)
   Base.write_atomic(stem.with_suffix(".json"),
     lambda f:f.write((json.dumps(meta,sort_keys=True,indent=2)+"\n").encode()))
   print("PASS_CHUNK",a.channel,a.part_index,a.parts,"input_keys",len(subset),
         "output_keys",len(out),"output_norm",Base.norm(out),flush=True)
  finally:restore()
  return 0
 # MERGE consumes ALL expected parts on one common source hash.
 total={};parts=[];diag=orig["leakage_max"];zero=0
 for part in range(a.parts):
  stem=checkpoint/f"ka_{a.channel}_part{part:03d}_of{a.parts:03d}"
  meta=json.loads(stem.with_suffix(".json").read_text())
  file=stem.with_suffix(".npz")
  expected_part_rows=rows[part::a.parts]
  if (meta.get("kb_sha256")!=orig["sha256"] or meta.get("freeze")!=freeze or
       meta.get("part_index")!=part or meta.get("parts")!=a.parts or
       meta.get("channel")!=list(ijk) or meta.get("sha256")!=Base.digest(file) or
       meta.get("rows")!=len(expected_part_rows) or
       meta.get("input_keys_sha256")!=key_digest([key for key,_ in expected_part_rows])):
   raise ValueError(f"ERROR: mismatched/missing checkpoint part {part}")
  state=Base.restore_state(file)
  if len(state)!=meta["support"] or not math.isclose(Base.norm(state),meta["norm"],abs_tol=1e-12,rel_tol=1e-12):
   raise ValueError(f"ERROR: corrupt part {part}")
  LP.add(total,state)
  diag=Base.maxdiag(diag,meta["leakage_max"])
  parts.append(meta)
 if sum(p["rows"] for p in parts)!=len(rows):
  raise ValueError("ERROR: missing input source keys")
 Base.save_result(checkpoint,"ka",ijk,total,freeze,diag,
          {"partitioned":True,"partition_count":a.parts,
           "source_kb_sha256":orig["sha256"],
           "note":"Floating-point order differs from monolithic; must check independent small-input/control."})
 print("PASS_MERGE",a.channel,"all",a.parts,"support",len(total),"norm",Base.norm(total))
 return 0
if __name__=="__main__":
 raise SystemExit(run())
