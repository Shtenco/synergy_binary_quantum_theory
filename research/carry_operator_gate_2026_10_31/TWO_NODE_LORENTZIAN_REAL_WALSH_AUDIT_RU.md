# REAL K5 Lorentzian adjacent-node diagonal block — 2026-10-31

**Status:** reconstructed **actual archived microscopic diagonal-neighbor matrix elements** at Jmax=7/2; **NOT** all outgoing Lorentzian columns, **NOT** full physical HDA, **NOT** a physical TT pole.

## Provenance and operator convention

Frozen upstream [Walsh amplitude JSON](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/verification_results/PETER_WEYL_LORENTZIAN_ENVTRACE_WALSH_NODE012.json), run 31836393368, source node 0, 24 oriented triples, Gauss-volume basis leakage 6.6944564e-16. The archived calculation measures only \(\langle K'_0,e|L_{raw,0}|K_0,e\rangle\) for diagonal neighboring environments, NOT \(\langle e'|L|e\rangle\) with \(e'\ne e\).

The already declared **candidate Hermitian completion**, \(\beta=\hbar=1\), is [LORENTZIAN_HERMITIAN_COMPLETION.md](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/LORENTZIAN_HERMITIAN_COMPLETION.md):
\[
H_{corr}=(16i/9)(L_{raw}-L_{raw}^\dagger)
=-\frac{32}{9}\,\frac{L_{raw}-L_{raw}^\dagger}{2i}.
\]

Fix K2=K3=K4=0. Reconstruct the actual four **projected** logical node0 amplitudes in basis K0,K1 = 00,02,20,22. The 10 raw Walsh coefficients reduce to:
\[
P_4 H_{phase,0} P_4=
0.33823953394649886\,Y_0I_1
-0.002352356009281853\,Y_0Z_1.
\]

Thus the normalized physical-candidate correction **on that compressed source block** is
\[
\boxed{P_4 H_{corr,0} P_4=-1.202629454031996\,Y_0I_1
+0.0083639324774466\,Y_0Z_1.}
\]

For neighbor Z1=+1, Y0 amplitude -1.1942655215545495; neighbor Z1=-1, -1.2109933865094424. Observed conditional difference =0.01672786495489298, relative split 1.39094090027567%. Eigen/singular moduli are each repeated twice: 1.1942655215545495,1.2109933865094424. Block rank 4, Frobenius norm2.405317076039484.

**New rigorous lower bound from real recorded matrix elements**: for every \(u=P_4u\) on the same domain and normalization,
\[
\boxed{\|H_{corr,0}u\|\ge\|P_4H_{corr,0}P_4u\|\ge1.1942655215545495\,\|u\|.}
\]
The **full outgoing vectors** themselves are not supplied in the archived diagonal packet. This finite cutoff source bound cannot be promoted to the physical continuum, gauge-invariant time or TT spectrum.

## Orientation sign vs logical Y sign

For the oriented four-simplex [01234], adjacent facets F0=[1234] and F1=-[0234] induce **opposite orientation** on their common face [234]. Exact fact \(\partial^2=0\). But this does **NOT** by itself prove a relative **Pauli-Y coefficient** for the two node operators: one must transport the logical intertwiner basis and face recouplings. The authentic **node1 Lorentzian raw block was not measured** in this packet.

## A mathematically rigorous data-identifiability obstruction

Knowing both neighbor-diagonal 2×2 blocks of any Hermitian 4×4 operator still leaves the complex 2×2 *offdiagonal neighbor-flipping block* unknown: **8 real free parameters**. E.g. \(H_{known}+tX_0X_1\) gives exactly the same measured diagonal-neighbor blocks for every real t. Furthermore its outgoing vectors outside P4 are wholly unmeasured. No unique complete 4-column output, relative Y0+Y1 vs Y0−Y1 sign, physical HDA, or Bell/TT observable is determined by the saved data.

**Conditional mathematical demonstration only:** if one assumes \(H_1=\eta(g_0 IY+g_1 ZY)\), where actual \(g_0=-1.2026294540,g_1=0.00836393248\), then
\[
-i[H_0,H_1]=2\eta g_0g_1(X_0Y_1-Y_0X_1).
\]
With externally stipulated \(\eta=-1\), coefficient +0.020117423097824168; this is **NOT measured microscopic HDA**.

## Next calculation determined by this result

1. **Lowest cost:** freeze the common K5 face orientation and SU2 recoupling map; run actual source1 24-term Lorentzian diagonal-neighbor contractions for K0=0,2, all other spectator nodes at K=0. Derive \(Y_1\), \(Z_0Y_1\), and orientation sign from the genuine source1 operator; do NOT extrapolate from source0.
2. Compute all **4 actual source columns per node** (8 total) including 24 forward raw and 24 adjoint contributions **without logical-output projection**. The older [full-state triple worker](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/scripts/peter_weyl_lorentzian_fullstate_triple_worker.py) currently hardcodes the first K5 state, so must be modified to accept the four K-pair source configurations and an independent adjoint execution. The [generic 24-way collector](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/scripts/peter_weyl_lorentzian_fullstate_24_collect_generic.py) already handles arbitrary source **node**, but it does NOT create missing columns or adjoint amplitudes.
3. Fail closed unless all complex outgoing keys and leakages are present; calculate full Gram / Schur overlap and then 10 independent full off-shell HDA commutators against separately constructed D, physically compatible projector/history and actual two-mode TT. No physicalization gate changed.

The locally downloaded package includes executable \`two_node_lorentzian_audit.py\`, frozen source coefficients, JSON output and a strict 8-column validator (which returns exit=2 BLOCKED on missing genuine data).