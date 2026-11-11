# BQG — audit of prior full-HDA and TT-pole claims (2026-10-12)

**STATUS: no completed theory-specific physical Lorentzian HDA or interacting graviton TT 1PI pole.** This is a reconciliation of the actual preceding results, not a new no-go for all BCQG.

Repository `main` snapshot: `9a58b9cde993ef9e352942345c56439c29261473`. Main tree 591 entries; 122 branches inventoried. Current first lines of [THEORY_STATUS.md](../../THEORY_STATUS.md) explicitly say structural candidate closed=false, physicalization complete=false, experimentally confirmed=false. [physicalization_gates.json](../../physicalization_gates.json) has **9/9 required physical gates open**, including `PHYSICAL_PROJECTOR_HISTORY`, `CONNECTED_INTERBLOCK_HISTORY`, `PHYSICAL_TT_KERNEL`, `IR_SIX_VECTOR`. [theory_gates.json](../../theory_gates.json) has 55 structural gates and one core `no_go`.

## HDA: results that must be preserved

1. Standard ADM test and simplex/flux bridge establish exact/classical target identities `c=1/2, AB=1`, `E_l=3V grad(lambda_l)` and K5 lapse cochain `omega_ij=N_i M_j-N_j M_i` **inside their declared ansätze**.
2. [THREE_NODE_GRAPH_HDA_RESULT.md](../../THREE_NODE_GRAPH_HDA_RESULT.md): genuinely spin/graph-changing Euclidean Peter-Weyl+route finite model, 3 pairs, `Delta(1/64)=0.0252238079`, fitted exponent~1.00644. This is not ten full Lorentzian pairs.
3. Research branch [PETER_WEYL_TWO_NODE_SINE_HDA_RESULT.md](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/PETER_WEYL_TWO_NODE_SINE_HDA_RESULT.md): preregistered sine ordering Euclidean+route `Delta(1/64)=0.02003033878`. CI run [31855735615](https://github.com/Shtenco/binary_quantum_theory/actions/runs/31855735615) successful. The paper explicitly lists missing full Lorentzian.
4. Research branch [FULL_OPERATOR_FIRST_HDA_CERTIFICATE.md](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/FULL_OPERATOR_FIRST_HDA_CERTIFICATE.md): valid **conditional** `Delta_full <= Delta_route + C1 eps + C2 eps²` at fixed cutoff and specified `N=Nbar+eps n`, WKB `Omega~eps^-1` plus bounded geometry, and conditional diagonal cutoff `Jmax~eps^-1/8`. This is not a proof of the off-shell HDA over arbitrary fixed-physical-momentum/refinement habitats. The `EL/LE/LL` amplitude workers timed out; the successful JSON gate certifies inputs/identities, not their computed output.
5. [October 10 audit](../../research/audit_2026_10_10/FIVE_ROUTES_AUDIT.md) reports `cross/D=0.266767` unchanged over checked L=32 to 96 at fixed physical momentum, illustrating why normalization by a diverging D cannot certify full physical HDA.
6. [Real Euclidean two-hit](../../research/audit_2026_10_10/results.json): one source column/pair, `||[H0,H1]psi||=2.8794538147`, 514 output states, `||Q32[H0,H1]psi||=2.82637117724`. A `Dpsi` artificially restricted to 32D has relative residual >=0.98156503251; this says NOTHING against an independently derived D with the correct outgoing support.
7. Research branch [LORENTZIAN_DIRECT_HERMITIAN_COMMUTATOR_THEOREM.md](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/LORENTZIAN_DIRECT_HERMITIAN_COMMUTATOR_THEOREM.md) reorganizes the full 24-term raw plus adjoint sum into twelve termwise-Hermitian commutator/anticommutator words. Algebraically valid and potentially efficient; still requires full operator amplitudes and domain checks.
8. [QUANTUM_HDA_HABITAT_RESIDUAL_PREREGISTRATION.md](https://github.com/Shtenco/binary_quantum_theory/blob/research/physical-scalar-kernel-dynamics/QUANTUM_HDA_HABITAT_RESIDUAL_PREREGISTRATION.md) specifies actual frozen full-`EE/EL/LE/LL` production, independent D, full domain hashes and refinement. **No accepted completed physical packet found in audited sources**.

## TT: results that must be preserved

1. [TT_PROPAGATOR_FIRST_PASS.md](../../TT_PROPAGATOR_FIRST_PASS.md): exact two-polarization **reduced/free** pole `4sin²(omega/2)=(1/3)sum_i4sin²(k_i/2)`, Gaussian propagator with free normalization `Z_T`. Bare quartic isotropic `-1/45`, cubic `-1/12`. Not physical theory-specific connected `Gamma_TT`.
2. [TT_VACUUM_TWO_POINT_RESULT.md](../../TT_VACUUM_TWO_POINT_RESULT.md): exact reduced Gaussian covariance `C(k)=1/(Z_T Omega sqrt(Omega²+4))`, IR ~k^-1. The older inference of TT vacuum spectrum from observer smoothing `b^-2` is explicitly rejected.
3. [TT_REGGE_ZT_L6_RESULT.md](../../TT_REGGE_ZT_L6_RESULT.md): preregistered held-out Regge Hessian test `Z6 predicted 0.11876923193907167 vs measured 0.11876075461190198`, relative 0.00714%; geometric unit-action normalization only.
4. [leading TT reference](https://github.com/Shtenco/binary_quantum_theory/blob/physical-tt-kernel-leading-ir/PHYSICAL_TT_KERNEL_LEADING_IR_DERIVATION.md): conditional ADM action and massless 2-TT kernel, **not measured connected theory-specific physical kernel**.
5. [five-block production](https://github.com/Shtenco/binary_quantum_theory/blob/physical-tt-kernel-leading-ir/BQG_MULTIBLOCK_TT_WILSON_PRODUCTION_STATUS.md): actual `P` carrier 960x30, rank30, Gram error4.53e-17; full `HE/S` output contractions, Q/Schur and physical history remain uncomputed.
6. [S4 quartic TT basis](../../S4_TT_QUARTIC_COMPLETE_BASIS.md): exact 6-dimensional TT quotient and invertible six-observable extractor. Does not provide the six **physical** Wilson values.
7. [6→1→0 universality](https://github.com/Shtenco/binary_quantum_theory/blob/research/physical-pole-universality/TT_POLE_UNIVERSALITY_NO_GO.md): conditional on restored SO3 and Lorentz metric-only massless branch. Not a dynamical derivation of Lorentz restoration from BQG.

## Missing single decisive chain

`full H_E + lambda_L H_L` (common graph-changing domain and Hermitian convention)
→ independently derived flux/route `D_ij` acting on nonconstant off-shell probes
→ all 10 complete `[H_i,H_j]/iℏ-D_ij` streams / fixed-momentum normalized Gram / refinement
→ nonzero compatible physical rigging projector `P_BQG`
→ physical connected `Z[J_g] -> W[J_g] -> Gamma[g]`
→ positive interacting massless TT pole / Ward / two helicities
→ frozen 6 Wilson values and one common scale / blind comparison.

In particular [BCQG_IR_HDA_UNIVERSALITY_CLOSURE.md](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/BCQG_IR_HDA_UNIVERSALITY_CLOSURE.md) has **full HDA among its assumptions**. It cannot prove its own hypothesis. The physical tensor frequency cannot be renamed from a constraint-resolvent spectral parameter.

## Scientific safeguards

- Physical TT is not obtained by identifying a constraint resolvent `z` with frequency `omega`.
- Graph-changing output outside the original 32D logical code is not an HDA anomaly by itself.
- CI `success` of a conditional theorem script does not mean its physical input assumptions have been generated.
- A `0.02` WKB residual where the reference D grows as `1/eps` does not demonstrate vanishing HDA at fixed physical momentum.
- Ten K5 pairs are sufficient **on one K5 with nodal lapse vectors only**; continuum smooth lapse and arbitrary graph refinement need extra work.
- No physicalization gate should be switched to proved based on this audit.

**Audit breadth:** all branches enumerated; main status and relevant foundational/TT/HDA/physicalization source documents checked, advanced branches deeply compared, historical GitHub Action successes verified. No claim that all 122 branch file bodies or expensive historical amplitudes were re-executed.
