# BQG: verified genuine G25 and coherent projection control — 2026-11-06

## Verified microscopic measurements

- Successful independent Actions run https://github.com/Shtenco/binary_quantum_theory/actions/runs/42130080910 (G25 v3).
- Four direct KA source jobs https://github.com/Shtenco/binary_quantum_theory/actions/runs/42128476913 (source 31,64,53,111), all success.
- Frozen KB parent SHA256: f486e2250d4a5678a34438467b9e338d730ebf6e1d1ec956d2de7f8c2521b293.
- All four raw NPZ outputs plus 21 previous true outputs are archived via full metadata, complex amplitudes and SHA256; every item is checked by G25 script before Gram computation.
- 25 of 120 nonzero KB source basis keys, one node1, one oriented triple (0,2,3), one auxiliary index channel 000, source column0, Jmax=7/2.
- Union of output SU(2)-keys: 1306.
- G25 numerical rank 25 of 25; smallest eigenvalue 4.1365110001389713e-7; largest 0.0014849597237780055. Gram condition number ~3589.8846, singular-value condition number ~59.9156.
- Norm of coherent sum AP25 psi=0.07151296436808603; incoherent quadrature=0.07594471320568569; coherent interference contribution to squared norm: -0.0006534953911827082.
- Input squared norm fraction 0.642142294353641; unprocessed 95 inputs squared-norm fraction 0.357857705646359; exact input residual norm 0.07684417184116164.
- Sufficient (NOT necessary) for 10% relative error to unknown full output: norm(delta95) <= norm(AP25 psi)/11 = 0.006501178578916912. No such estimate of delta95 has been established. Input 64.214% coverage is not output accuracy.

## Independent charge-projection control

Real one-hit SU2 job https://github.com/Shtenco/binary_quantum_theory/actions/runs/42129284121: Gauss input charge-projector leak 4.440892098500625e-16; genuine single holonomy hit generated 2 branches; coherent preprojected squared norm 0.5000000000000011, postprojected squared norm 0.5000000000000016, coherent post-sum leak zero within floating-point accuracy.

This is ONE hit, NOT the complete HE_sine or full HL. Previous complete_charge_basis_leakage maximum 1.0 over primitive CK branches is still NOT resolved for the full operator. Do not call this a full charge-closure proof.

## Correct next workload

Real next high-weight KB indices [15,39,40,80], weights [0.0183521,0.0162602,0.0134731,0.0134426] relative to total input squared norm. New 4 source jobs at https://github.com/Shtenco/binary_quantum_theory/actions/runs/42130349763. If all succeed and are independently reconstructed: 29/120 directions, input squared-weight coverage ~0.70367017546, 91 remaining basis directions, residual input norm ~0.06992676384. These are predicted INPUT quantities only.

Full 8 auxiliary paths, 24 oriented triple sum, Hermitian adjoint, node0 full output, independently constructed off-shell D and HDA, physical group-averaged projector and interacting massless TT remain OPEN.

## Important methodological NO-GO

No upper bound on unknown coherent delta95 follows from only the 25 known columns of A. A rank-one extension of A supported entirely on the unmeasured 95-dimensional subspace can reproduce all measured columns while changing delta95 arbitrarily. A physical error bound therefore requires additional structural operator bounds or the actual missing coherent output. This is an identifiability statement, not a physical no-go for BQG.
