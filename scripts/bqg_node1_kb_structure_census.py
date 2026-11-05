#!/usr/bin/env python3
"""Cheap authentic KB phase and symmetry-block census, NOT invented symmetry."""
import collections,hashlib,json,math
from pathlib import Path
import numpy as np
p=Path("checkpoint/kb_000.npz")
m=json.loads(Path("checkpoint/kb_000.json").read_text(encoding="utf-8"))
assert hashlib.sha256(p.read_bytes()).hexdigest()==m["sha256"]=="f486e2250d4a5678a34438467b9e338d730ebf6e1d1ec956d2de7f8c2521b293"
with np.load(p,allow_pickle=False) as z:
    spins=z["spins"];others=z["Kother"];charges=z["charges"];amps=z["amp"]
assert len(amps)==120
spin_blocks=collections.Counter(tuple(int(x) for x in row) for row in spins)
environment=collections.Counter(tuple(int(x) for x in row) for row in others)
charge_blocks=collections.Counter(tuple(int(x) for x in row) for row in charges)
real_max=float(np.max(np.abs(amps.real)))
nrm=float(np.linalg.norm(amps))
assert real_max<1e-12
assert len(spin_blocks)==34 and len(environment)==5 and len(charge_blocks)==13
report=dict(status="PASS_GENUINE_FROZEN_KB_STRUCTURAL_CENSUS",
  original_KB_sha256=m["sha256"],source_node=m["freeze"]["source"],input_count=120,
  unique_spin_patterns=len(spin_blocks),unique_environment_K_patterns=len(environment),
  unique_charged_local_patterns=len(charge_blocks),
  environment_counts={str(k):v for k,v in sorted(environment.items())},
  max_abs_real_component=real_max,input_norm=nrm,
  complexity_note="Potential caching by common local spin/K/charge ingredients only; NO symmetry transports of uncomputed CK outputs have been proved.",
  physical_HDA="OPEN",TT="OPEN")
print("BQG_KB_STRUCTURAL_AUDIT_BEGIN")
print(json.dumps(report,indent=2))
print("BQG_KB_STRUCTURAL_AUDIT_END")
