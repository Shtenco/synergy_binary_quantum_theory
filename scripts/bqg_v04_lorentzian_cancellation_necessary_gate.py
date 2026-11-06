#!/usr/bin/env python3
"""Rigorous necessary Lorentzian-cancellation gate from REAL j=2 W5 Euclidean K5 data.

A u := (H_E,0 u,..., H_E,4 u) in common direct-sum outgoing Hilbert.
B u := (H_L,0 u,..., H_L,4 u) MUST be independently computed on SAME W5,
normalization, spin cutoff, and full open output. No B measured here.

singular-value Weyl bound:
sigma_min(A+lambda B) >= sqrt(lambda_min(G_EE))-|lambda|*||B||.
If Lorentzian norm < bound, then NO physical common zero in THIS W5.
"""
import argparse,json,math
from pathlib import Path
import numpy as np
FILE=Path("research/node1_fullstate_gate_2026_10_31/V04_CROSS_BRANCH_REAL_EUCLIDEAN_W5_LOG_EVIDENCE.json")
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--lambda-abs",type=float)
 args=p.parse_args()
 data=json.loads(FILE.read_text())["data"]
 assert data["screening_conclusion"]=="NO_COMMON_EUCLIDEAN_ZERO_IN_W5"
 assert data["claim_scope"].startswith("Actual j=2 K5")
 G=np.array([[complex(x[0],x[1]) for x in row] for row in data["joint_euclidean_gram"]])
 assert G.shape==(5,5)
 assert np.linalg.norm(G-G.conj().T)<1e-10
 evals=np.linalg.eigvalsh(G)
 assert abs(float(evals[0])-data["min_eigenvalue"])<1e-8
 assert evals[0]>0
 smin=float(math.sqrt(evals[0]))
 result={
 "status":"RIGOROUS_NECESSARY_CONDITION_FROM_ACTUAL_EUCLIDEAN_W5",
 "Euclidean_W5_master_spectrum":[float(x) for x in evals],
 "sigma_min_actual_five_vertex_HE":smin,
 "necessary_condition_for_full_Lorentzian_kernel_in_W5":"abs(lambda_L)*||stack_v H_L,v restricted to W5|| >= sigma_min(A)",
 "necessary_lorentz_norm_times_coupling_min":smin,
 "joint_master_must_be_measured":"M(lambda)=G_EE+lambda*(G_EL+G_EL^dagger)+lambda^2*G_LL",
 "required_future_blocks":["same-space full Lorentzian output columns B on W5","cross Gram A^dagger B","Lorentzian Gram B^dagger B","independent D and HDA on common habitat"],
 "scientific_limit":"No physical projector or full HDA proved; a kernel outside W5 is not excluded."}
 if args.lambda_abs is not None:
  if args.lambda_abs<=0:raise RuntimeError("Need positive nonzero abs lambda")
  result["requested_abs_lambda"]=args.lambda_abs
  result["necessary_lower_bound_on_HL_operator_norm"]=smin/args.lambda_abs
 print("BQG_REAL_W5_LORENTZIAN_CANCELLATION_BOUND_BEGIN")
 print(json.dumps(result,indent=2))
 print("BQG_REAL_W5_LORENTZIAN_CANCELLATION_BOUND_END")
if __name__=="__main__":main()
