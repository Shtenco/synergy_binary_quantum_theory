# BQG: SO(5) graph-link to coherent carry, exact no-go and constructive chirality (2026-10-31)

**Scientific status:** finite **KINEMATIC OPERATOR THEOREM + restricted Hamiltonian no-go**. The exact microscopic full-Lorentzian \(H^E+\lambda_L H^L\), physical projector, independently measured HDA, true dynamical 2T and massless TT pole **remain OPEN**. Do not alter physicalization gates.

## Source provenance

Main:
- [SO(5) vector link representation](../../scripts/su2_quantum_link_vector5_gate.py): Hermitian \(U_a=M(a,4),\ a=0..3\), 4 active bi-doublet + one j=0 no-link state.
- [Graph-link Peter–Weyl bridge](../../scripts/q2_graphlink_peter_weyl_gate.py): exact \(P_g U_a P_0 U_b P_g=|a\rangle\langle b|\).
- [2T frontier](../../README_2T_FRONTIER.md): no actual microscopic \(Sp(2,R)\) embedding proved.

History branch [phase refinement and winding](https://github.com/Shtenco/binary_quantum_theory/blob/research/q2-history-phase-refinement-winding/Q2_HISTORY_PHASE_REFINEMENT_AND_WINDING.md) selects frozen oriented Gray order \(00\to01\to11\to10\to00\); first \(C_4\to C_8\) **history lift** exact but dynamical all-level carry selection OPEN.

## Theorem I: one no-link channel cannot implement unconditioned coherent carry

For \(P_0=|0\rangle\langle0|\), \(U_a=-i(|a\rangle\langle0|-|0\rangle\langle a|)\), any *single* Hermitian transporter Hamiltonian \(H(w)=\sum_{a=0}^3w_aU_a\) with real \(w\) has:

\[
P_g H(w)P_0 H(w)P_g=|w\rangle\langle w|,\quad \operatorname{rank}=1.
\]

The directed \(C_4\) shift \(S_4|a_k\rangle=|a_{k+1}\rangle\) is unitary rank4. By Eckart–Young:

\[
\boxed{\inf_{\operatorname{rank} B\le r}\frac{\|S_4-B\|_F}{\|S_4\|_F}=\sqrt{1-r/4};\quad r=1\Rightarrow\sqrt3/2=0.8660254037844386.}
\]

The frozen **undirected** Hamming adjacency \(X\otimes I+I\otimes X\) is rank2; its best rank1 relative Frobenius residual is \(1/\sqrt2\).

For \(\|w\|=1\) the full **actual-transporters restricted model** has exactly
\[
P_g e^{-itH(w)}P_g=I_4+(\cos t-1)|w\rangle\langle w|;
\]
three orthogonal dark directions remain unchanged. SciPy matrix exponentials confirm at times \(0,.13,.67,1.7,5\).

**Stronger Hermitian obstruction:** for every Hermitian \(H\), every \(P,Q\) with \(P+Q=I\),
\(PHQHP=(QHP)^\dagger(QHP)\succeq0\).
Since \(S_4\) is not Hermitian or positive, it cannot **equal a two-operator-insertion PSD return** even if \(\dim Q\) is increased. The closest arbitrary PSD Hermitian 4x4 matrix is a Frobenius distance \(\sqrt3\) (relative \(\sqrt3/2\)) away from \(S_4\). This is **not** a no-go for unitary time evolution \(\exp(-itH)\), different *ordered* elementary insertions, history memory or the full graph-changing BQG dynamics.

## Theorem II: real micro-transporters already contain a directed operator

The same \(U_a\), the already frozen Gray cycle \(a_k=(0,1,3,2)\), and the common no-link projector yield the exact Hermitian operators:

\[
K_\chi=\frac1{2i}\sum_{k=0}^3P_g(U_{a_{k+1}}P_0U_{a_k}-U_{a_k}P_0U_{a_{k+1}})P_g,
\]
\[
K_s=\frac12\sum_{k=0}^3P_g(U_{a_{k+1}}P_0U_{a_k}+U_{a_k}P_0U_{a_{k+1}})P_g.
\]

On the active 4D code,
\[
\boxed{K_\chi=(S_4-S_4^\dagger)/(2i),\quad K_s=(S_4+S_4^\dagger)/2,}
\]
\[
\boxed{S_4=K_s+iK_\chi,\quad [K_s,K_\chi]=0,\quad K_s^2+K_\chi^2=I_4,\quad K_sK_\chi=0.}
\]

\(\mathrm{spec}(K_\chi)=\mathrm{spec}(K_s)=\{-1,0,0,+1\}\); all identities are exactly reproduced using the original \(M(a,4)\) matrices. This is a **constructive operator available in BQG's declared kinematic algebra**, NOT a term measured or shown to be dynamically selected by its full Lorentzian constraint. \(K_\chi\) alone cannot generate the complete shift, since it fails to distinguish the two zero-eigenvalue states on which \(S_4\) is +1/-1.

## Theorem III: intermediate *history* must remember the edge

To implement an exact **discrete-time** 8-cycle, the four event states \(|m_k\rangle\) must remain distinguishable. An ideal 8D unitary \(T_8\) defined by \(|a_k\rangle\mapsto|m_k\rangle\mapsto|a_{k+1}\rangle\) has \(T_8^8=I\), \(\operatorname{order}(T_8)=8\), and \((T_8^2)_{active}=S_4\). If all four events are coherently erased to the identical normalized one-dimensional \(|0\rangle\), any path-uninformed decoder succeeds with average probability at most 1/4 for four equiprobable different orthogonal targets. This information-erasure bound is avoided by a separate path/event memory. Independent \(Z_4\times Z_2\) has maximum element order4, not8.

## Fast decisive next *physical* computation

1. Freeze an **independent** embedding of these four Gray-history channels into the real K5 Peter–Weyl / 5-block q4 Hilbert habitat. An SO(5) link identity alone is not this embedding. Specify full \(E_i+\lambda_LL_i\), 24 Lorentzian epsilon terms, gauge, common domains and output labels.
2. For two adjacent collective blocks, calculate the **four full outgoing vectors** \(QH_vP_{\mathrm{Gray}}\) on the *same* undropped source and outgoing basis. Stream amplitudes. Report rank, singular values, leakage, cutoff and refinement. If rank<4 on a claimed four orthogonal event-memory carrier, exact deterministic carry there is impossible. Rank4 is necessary, not sufficient.
3. Freeze orientations *before looking at data*, then compute ordered pairs and antisymmetric \([H_i,H_j]/(2i)\), compare their **projected** channel with \(K_\chi\) and orthogonal residual. Independently construct \(D\) and test 10 pair HDA relations; the chiral projection does **not** define \(D\).
4. Only if the real physical history/rigging limit and reflection-positive transfer exist, evaluate Bell with **two commuting Dirac algebras** and physical TT \(C_0,C_1,C_2\), positive massless pole and refinements. No C8 parameter or time-crystal subharmonic is silently relabelled physical frequency.

**Machine outcomes** for copied exact \(M(a,4)\) reference: full 16 matrix-unit identities discrepancy zero; rank1 unconditioned; S4 rank4; Kχ eigen ±1,0,0; Kχ reconstructed error zero; unitary T8 order8; closest rank1 and PSD error sqrt3. Local downloadable package also contains detailed Python test, JSON and Russian report.
