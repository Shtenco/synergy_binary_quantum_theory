# BQG: exact G33 K5 Lorentzian micro-amplitude audit (2026-11-06)

Real upstream 4 CK outputs run https://github.com/Shtenco/binary_quantum_theory/actions/runs/42131670594 SUCCESS. Real independent 33-column Gram and coherent increment run https://github.com/Shtenco/binary_quantum_theory/actions/runs/42132361925 SUCCESS. Every input uses frozen genuine KB source SHA f486e2250d4a5678a34438467b9e338d730ebf6e1d1ec956d2de7f8c2521b293. All new outgoing NPZ evidence pinned in research/node1_fullstate_gate_2026_10_31/log_recovered_ka and verified by SHA in the G33 code.

## Measured actual G33

One physical-source microscopic *subamplitude*: node1, source column0, ordered T_(0,2,3), auxiliary 000, Jmax=7/2. NOT the full 8 aux sum, 24 orientation sum, or Hermitian H_L.

G33 rank 33/33, union 1407 outgoing covariant SU(2) keys; lambda_min=4.120535987166385e-7 and lambda_max=0.0014862687223758581. Gram condition ~3606.9791. Coherent output norm=0.07578186266888405, quadrature norm=0.07945770017102775, signed squared-norm interference=-0.0005706354069033425. Processed fraction of INPUT squared norm=0.754980618194833, 87 uncomputed input sources. True KB remainder input norm=0.06358518918942349, which does NOT certify an OUTPUT error.

## Crucial vector convergence result

G29 norm=0.07324866404521364. G33 norm increased only +3.4584%, but ||Psi33-Psi29||=0.019130397993351585, i.e. 25.2440% of ||Psi33||. Normalized real cross-overlap of new increment with prior Psi29 =0.004121884976388552, almost real-orthogonal. 87 new output keys; support count saturation does NOT certify vector convergence. A sufficient 10% full-output relative-error criterion is ||A Q87 psi||<=||Psi33||/11=0.006889260242625823; this tail has NOT been measured.

## Actual recoupling input structure, verified

Run https://github.com/Shtenco/binary_quantum_theory/actions/runs/42132264925 SUCCESS: 120 authentic KB source records, 34 distinct local edge-spin patterns, FIVE environment-intertwiner K patterns (counts 8,28,28,28,28), 13 local charged-label patterns, real component max 4.095565252728536e-17; entire source essentially purely imaginary. A possible optimization is memoizing repeated LOCAL SU(2) recoupling ingredients keyed by full physical local geometry. No symmetry transport of distinct sources is proved from these counts.

## Negative performance control

Run https://github.com/Shtenco/binary_quantum_theory/actions/runs/42130832847 SUCCESS: isolated caches 35.88715 s vs one shared in-process cache 35.85226 s on two authentic cheap sources [0,12]. Speedup ~0.1% only. Do not assume the existing lru_cache solves 87 high-spin cost.

## One-hit Gauss/charged closure only

Real one-hit projection control https://github.com/Shtenco/binary_quantum_theory/actions/runs/42129284121 PASS on two branches with pre- and post-projector squared norms ~0.5. This does NOT resolve max primitive CK complete_charge_basis_leakage=1.0 at the full HE_sine level. Need exact coherent POST-SUM charge leakage, with no unproven extra projection.

## Next physical gate

1. Complete the coherent charged-basis closure of actual HE_sine on same common graph-changing habitat BEFORE using projected amplitudes for full physical projector claims.
2. Construct independently defined D on full outgoing states and confirm full Lorentzian off-shell HDA defect; no HDA inference from Gram rank.
3. Compute all 8 auxiliary channels and 24 oriented source terms, adjoints, physical history/projector, connected gauge-invariant TT with positive massless pole. These are OPEN.

The kernel of a physical master constraint and TT pole do NOT follow from finite G33 positivity alone.
