#!/usr/bin/env python3
"""Cheap K5 full-state input grid + completeness preflight. No fake amplitudes."""
import argparse,itertools,json,pathlib,sys
K=list(itertools.product((0,2),repeat=5))
PAIRS=[0,8,16,24]
assert [K[i][:2] for i in PAIRS]==[(0,0),(0,2),(2,0),(2,2)]
assert all(K[i][2:]==(0,0,0) for i in PAIRS)
def table(source):
    neigh=[v for v in range(5) if v!=source]
    ans={}
    for omit_index,omit in enumerate(neigh):
        base=tuple(v for v in neigh if v!=omit)
        for p in itertools.permutations(base):
            inversions=sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
            ans["".join(map(str,p))]=(-1)**(omit_index+inversions)
    assert len(ans)==24
    return ans
def check(root):
    missing=[];bad=[]
    for source in (0,1):
      for column in PAIRS:
        for tag,sign in table(source).items():
          stem=root/f"node{source}"/f"column{column}"/f"T_{tag}"
          jp=stem.with_suffix(".json");sp=stem.with_suffix(".npz")
          if not (jp.is_file() and sp.is_file()):
              missing.append(f"{source}:{column}:{tag}");continue
          d=json.loads(jp.read_text(encoding="utf-8"))
          if d.get("column_index")!=column or d.get("source_node")!=source or d.get("epsilon_coefficient")!=sign or not d.get("passed",False) or d.get("input_Ks")!=list(K[column]):
              bad.append(f"{source}:{column}:{tag}")
    return {"status":"READY" if not (missing or bad) else "BLOCKED",
            "sources":2,"columns_per_node":len(PAIRS),"triples_per_column":24,
            "required_forward_raw_packets":2*4*24,
            "note":"The 192 raw forward packets do not by themselves establish 192 direct adjoint columns. Hermitian completion requires independent adjoint action or full matrix-element identity.",
            "missing_count":len(missing),"bad_count":len(bad),
            "missing_examples":missing[:12],"bad_examples":bad[:12]}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=pathlib.Path,default=pathlib.Path("artifacts/lorentzian_columns"))
    args=p.parse_args()
    d=check(args.root)
    print(json.dumps(d,ensure_ascii=False,indent=2))
    return 0 if d["status"]=="READY" else 2
if __name__=="__main__":sys.exit(main())
