# BQG: conditional Lorentzian square/triangle loop phases, relational Bell/gluing and physical-projector no-go (2026-10-31)

**Status: exact symbolic mathematical results and reduced-model controls PASS. Full microscopic Lorentzian BQG, physical HDA and interacting TT pole remain OPEN.** No physicalization status changed.

## Verified but distinct upstream objects

1. Main [PETER_WEYL_LOGICAL_ANISOTROPY.md](../../PETER_WEYL_LOGICAL_ANISOTROPY.md): unbiased real Euclidean two-node RETURN matrix \(K_{\rm ret}=P(H_0^E+H_1^E)^2P\), NOT a physical Hamiltonian or Schur kernel. In its original pre-rotation computational basis, \(p=K_{01}=K_{02}=-0.09188432819822412\), \(q=K_{13}=K_{23}=0.2529732671985883\), and the two-bit \(00\leftrightarrow11\) edge has coefficient \(d=-1.7250844497317166\).
2. Legacy [Peter–Weyl Lorentzian environment-trace](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/PETER_WEYL_LORENTZIAN_ENVTRACE_RESULT.md): real raw \(L_{\rm raw}=icY\), \(c=1.3389293521464034\), onebody \(J_{\max}=7/2\), 24-oriented term orbit and 16 environments. The [phase certificate](https://github.com/Shtenco/binary_quantum_theory/blob/research/bcqg-core-candidate-v1/LORENTZIAN_COMMUTATOR_PHASE_CERTIFICATE.md) gives a **conditional** \(-i\) completion \(cY\) in its chosen nested-commutator convention. Final real scale, sign, two-node operator and common physical domain are not fixed.
3. Main [two-node reduced gluing](../../BQG_TWO_NODE_REDUCED_GLUING_RESULT.md): exact *explicitly modelled* \(H_g=\frac32 I-\frac34(X_AX_B+Z_AZ_B)\), unique zero \(\Phi^+=(|00\rangle+|11\rangle)/\sqrt2\); neither full gravitational constraint nor physical rigging/history.

**Do NOT sum object (1) with object (2) as though it were a measured physical Hamiltonian:** different cutoffs, \(H_E^2\) versus \(H_L\), unconstrained resolvent, and physical time/normalization not established.

## Gate 1 — Gauge-invariant closed-loop diagnostic (exact algebra)

For the purely DIAGNOSTIC hybrid \(H_{\rm test}=K_{\rm ret}+aY_A+bY_B\), orient the frozen Gray square \(0\to1\to3\to2\to0\) and set \(W_\square=H_{1,0}H_{3,1}H_{2,3}H_{0,2}\). It is invariant under arbitrary independent diagonal basis rephasings \(H\mapsto UHU^\dagger\).

\[
\boxed{W_\square=(p+ib)(q+ia)(q-ib)(p-ia),}
\qquad
\boxed{\Im W_\square=(q-p)(b-a)(pq-ab).}
\]

For symmetric local Lorentzian Y \(a=b\), the **Gray-square imaginary flux vanishes exactly**, even if the local Y matrix elements are imaginary. This DOES NOT imply a real gauge exists for the complete four-state graph: the two-bit diagonal edges support triangles. The triangle \(0\to1\to3\to0\) has \(W_\triangle=d(p+ib)(q+ia)\). For \(a=b=c\), \(W_\triangle=3.1327121118-0.3720777873i\), phase \(-0.118217964\) rad while \(W_\square=3.3442899112+0i\). For \(a=c,b=-c\), \(W_\square=2.9178823358-1.6340861925i\), phase \(-0.5105071283\) rad. These angle values are **not physical predictions**, only conditional hybrids using source-reported raw numbers.

The sign flip \((a,b)\mapsto(-a,-b)\) conjugates all loop products; orientation-unselected averages may erase imaginary phases. No exogenous timed Floquet protocol is inferred from a static Wilson phase.

## Gate 2 — Exact relational orientation selection from source reduced gluing

Let \(P_\pm=\frac12(I\pm Y_AY_B)\). Algebraically:
\[
[H_g,Y_AY_B]=0,\quad
(Y_A+Y_B)P_-=0,\quad (Y_A-Y_B)P_+=0,
\quad Y_AY_B|\Phi^+\rangle=-|\Phi^+\rangle.
\]

Thus symmetric onebody Y cannot act within the anti-aligned orientation sector. Antisymmetric Y does. Yet \(c\) and the relative sign of true neighbor Lorentzian couplings must be **derived** from orientation/gluing/ordering, not selected because a test produces chirality.

- \(H_g+c(Y_A+Y_B)\) has exact eigenvalues \(0,3,3/2\pm2c\); Bell \(\Phi^+\) remains an exact zero eigenvector, but ceases to be the lowest-energy vector for \(|c|>3/4\).
- \(H_g+c(Y_A-Y_B)\) restricted to \(\{\Phi^+,\Psi^-\}\) has exact matrix \(\left(\begin{smallmatrix}0&2ic\\-2ic&3\end{smallmatrix}\right)\), energies \((3\pm\sqrt{9+16c^2})/2\), remaining doubly-degenerate \(3/2\).

## Gate 3 — SAME reduced ground state supplies Bell and TWO shape-source histories, but not TT

For the *model* staggered H and a chosen raw-reference \(c=1.3389293521464034\), exact:
\[
\mathcal C=\frac3{\sqrt{9+16c^2}}=0.4887023577,\quad
\langle Y_A\rangle=-\frac{4c}{\sqrt{9+16c^2}}=-0.8724505749,
\quad\mathcal C^2+\langle Y_A\rangle^2=1.
\]
Pure-state Horodecki \(S_{\max}=2\sqrt{1+\mathcal C^2}=2.2260548011\), confirmed by direct correlation matrix. The same state, same \(H-E_0\), source vectors \(U=(X_A\Omega,Z_A\Omega)\) yield
\[
\operatorname{spec}C_0=\{1-|\langle Y_A\rangle|,1+|\langle Y_A\rangle|\}>0,
\quad C_t=e^{-t\Delta}C_0,\quad\Delta=\frac{\sqrt{9+16c^2}}2=3.0693529025.
\]
\(C_2-C_1 C_0^{-1}C_1=0\) within numerical \(3\times10^{-17}\). The response is GAPPED and NO spatial momentum or tensor gauge Ward identities were derived. These are shape excitations, NOT physical TT graviton.

## Gate 4 — Rank-one physical gluing projector FAILS both shape/Dirac probes

For the original pure gluing ground \(P_g=|\Phi^+\rangle\langle\Phi^+|\),
\(U=(X_A\Phi^+,Z_A\Phi^+)\):
\[
\boxed{U^\dagger U=I_2,\qquad U^\dagger P_g U=0_2.}
\]
A one-dimensional physical-algebra compression has no two independent Dirac-TT source modes, and the kinematic Bell violation is not a Bell test of two nontrivial commuting physical Dirac algebras. This does NOT rule out a larger fully graph-changing BQG Hilbert or relational boundary excitations.

## Automatic next full-microscopic calculation

**First measurement: real two-node Lorentzian SIGN and connected operator** on a common Peter–Weyl / Gauss habitat, with complete 24-term epsilon-oriented Hermitian columns, all outgoing states and fixed action normalization. Compute coefficients of \(Y_A+Y_B\), \(Y_A-Y_B\), and all correlated chiral channels; fix relative PL face signs *before* seeing matrix elements. Do not project away Q leakage or infer symmetry from a onebody environment trace.

**Second (only after 1):** derive a legitimate Schur/Feshbach return operator with independently fixed denominator from same \(\mathcal C^\dagger\mathcal C\), then evaluate gauge-invariant loop products and their refinement. A rank2-or-better physical TT source Gram is required; gluing rank1 is a NO-GO for two physical shape sources in this reduced model.

**Independent:** full Lorentzian off-shell 10-pair HDA, a nonzero compatible rigging/history, physical 2-polarization \(C_0,C_1,C_2\), massless positive TT pole and Ward identities. A nonzero triangle loop is NOT an HDA or TT proof.

Reproduction: user-facing verified package \`BQG_CHIRAL_RELATIONAL_20261031.zip\` with \`chiral_relational_gate.py\` (SymPy + NumPy), detailed report and JSON. This GitHub note is an abbreviated scientific evidence ledger.
