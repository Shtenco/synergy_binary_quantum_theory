#!/usr/bin/env python3
"""K5 q2 32D Euclidean source gap and Lorentzian parity lower bound.
Optional re-evaluation from the real committed fine_H0_full_images.npz.
This does NOT compute full physical Lorentzian HDA or a physical TT pole.
"""
from pathlib import Path
import math,json,argparse,numpy as np
REPORTED_GAP=3.97703147684868
def csr_packet(path):
 with np.load(path,allow_pickle=False) as z:
  fmt=z["format"].item()
  if isinstance(fmt,bytes):fmt=fmt.decode()
  assert fmt=="csr"
  shape=tuple(map(int,z["shape"]))
  idx,ptr,data=z["indices"],z["indptr"],z["data"]
  assert shape[1]==32 and len(ptr)==shape[0]+1 and ptr[-1]==len(idx)==len(data)
  H=np.zeros(shape,complex)
  for i in range(shape[0]):
   np.add.at(H[i],idx[ptr[i]:ptr[i+1]],data[ptr[i]:ptr[i+1]])
 return H
def check(npz=None):
 a=math.sqrt(REPORTED_GAP)
 pi=np.diag([1.,-1.])
 E=a*np.array([[0.,1.],[1.,0.]])
 L=a*np.eye(2)
 H=E+L
 assert np.allclose(pi@E@pi,-E)
 assert np.allclose(pi@L@pi,L)
 assert np.allclose(np.linalg.eigvalsh(H),[0,2*a])
 assert np.allclose(np.linalg.norm(H[:,0])**2,2*REPORTED_GAP)
 ker=np.array([1.,-1.])/math.sqrt(2)
 assert np.allclose(H@ker,0) and abs(ker[0]**2-.5)<1e-12
 theta=.2
 A=np.array([[a],[0.]])
 B=np.array([[math.cos(theta)],[math.sin(theta)]])
 f=float((A.T@(np.eye(2)-B@B.T)@A)[0,0])
 assert abs(f-REPORTED_GAP*math.sin(theta)**2)<1e-12
 result={"reported_real_Euclidean_H0_gap_squared":REPORTED_GAP,
 "code_only_full_Lorentzian_gap_squared_at_least":REPORTED_GAP,
 "code_only_full_Lorentzian_gap_norm_at_least":a,
 "toy_full_kernel_overlap_with_code":.5,
 "toy_schur_principal_angle_residual":f,
 "status":"VALID scoped mathematical theorem, NOT full BQG HDA or physical TT"}
 if npz:
  h=csr_packet(npz);vals=np.linalg.eigvalsh(h.conj().T@h)
  result["actual_packet_shape"]=list(h.shape)
  result["actual_packet_lambda_min"]=float(vals[0])
  result["actual_packet_rank"]=int(np.linalg.matrix_rank(h))
  assert abs(vals[0]-REPORTED_GAP)<5e-6
 return result
if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--npz",type=Path)
 a=p.parse_args()
 print(json.dumps(check(a.npz),indent=2))
