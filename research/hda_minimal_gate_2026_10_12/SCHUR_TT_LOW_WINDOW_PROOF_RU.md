# BQG: new exact Schur + low-window source proof (2026-10-12)

**Science status:** precise algebraic reduction and toy exact tests **PASS**; **full microscopic BQG normal projector/HDA/interacting TT not computed.** The existing five-block microscopic kinematic carrier P(960×30) has orthonormal rank30, but complete Hermitian Lorentzian H_E/H_L outgoing columns, complement Q closure and genuine physical history are still missing (see \`BQG_MULTIBLOCK_TT_WILSON_PRODUCTION_STATUS.md\` in branch \`physical-tt-kernel-leading-ir\`).

Let C denote the vertical stack of ALL actual constraint operators on the same Hilbert domain; M=C†C. In P/Q block form set X=C P, Y=C Q, A=X†X, B=X†Y, D=Y†Y. If D invertible, the **exact** reduced zero-constraint form is

\[
F=A-BD^{-1}B^\dagger=X^\dagger(I-P_{\mathrm{ran}Y})X\succeq0.
\]

For any trial retained vector u,

\[
u^\dagger Fu=\min_q\|Xu+Yq\|^2.
\]

When Y has null vectors, finite-dimensional Moore-Penrose inverse gives the same minimal residual; track Q-only zero states separately, use frozen tolerance and full outgoing basis. If E is a full orthonormal basis of ker F, K=D+ B†, the compressed true finite normal-zero projector is

\[
P P_{\ker M}P
=E[E^\dagger(I+K^\dagger K)E]^{-1}E^\dagger.
\]

For two sources U inside the P carrier their zero-mode Gram is \(U^\dagger E[E^\dagger(I+K^\dagger K)E]^{-1}E^\dagger U\). This does not assert Dirac-invariant TT sources or a physical Lorentzian history.

## Rational exact control

Take P dimension3, Q dimension2, \(D=\mathrm{diag}(2,3)\), \(K=\begin{pmatrix}1/3&0&1/2\\0&2/5&1/3\end{pmatrix}\), \(F=\mathrm{diag}(0,0,5)\), \(M=\begin{pmatrix}F+K^TDK&K^TD\\DK&D\end{pmatrix}\).

SymPy **exact** computations:
- positive M is rank3 with two zero modes; positive compressed A is rank3, \(\lambda_{\min}(A)=0.20026822778152817\) in numerical evaluation;
- completion by one Q direction changes effective rank 3→2; both Q directions make 2→1, revealing TWO exact normal zero modes;
- for two P-source vectors e1,e2, \(G_{\mathrm{toy}}^{zero}=\mathrm{diag}(9/10,25/29)\) is positive;
- deleting THREE constraint-output coordinates falsely creates an extra zero mode (truncated rank0 F vs true rank1).
The toy is not the BQG microscopic operator.

## Even cheaper sufficient two-source theorem

For any positive self-adjoint M, \(\delta>0\), \(P_\delta=\mathbf1_{[0,\delta]}(M)\), define \(K_0=U^\dagger U\), \(K_1=U^\dagger M U=\sum_v(H_vU)^\dagger(H_vU)\). Spectral calculus gives

\[
\boxed{U^\dagger P_\delta U\succeq K_0-\delta^{-1}K_1}.
\]

Thus \(\lambda_{\max}(K_0^{-1/2}K_1K_0^{-1/2})<\delta\) guarantees rank TWO in the **normal-master low-energy window**. It is sufficient only, cannot certify a true Dirac observable or massless TT propagation. In the toy \(K_1=\mathrm{diag}(2/9,12/25)\), \(\delta=1\) gives bound \(\mathrm{diag}(7/9,13/25)\). Actual zero projected Gram is \(\mathrm{diag}(9/10,25/29)\), and first positive eigenvalue1.878309>1, so the window equals exact zero projector in this toy.

## Automatic next calculation

1. Fast: construct two properly chosen collective TT source columns U from measured six-component metric carrier; apply frozen **complete** Hermitian \(H_v^E+\lambda_L H_v^L\) on every node, all unprojected output states. Compute 2x2 K0,K1, validate normalized inequality for preregistered \(\delta_n\).
2. If inconclusive: use actual 30 retained columns and an adaptively grown Q; stream complete complex outgoing amplitudes for five node constraints, and compute the true 30x30 Schur Gram without a dense ambient Hilbert matrix. Outgoing truncation can create false physical kernels, so fail closed on leakage.
3. Candidate normal master kernel then needs independent spatial-diffeomorphism D/HDA, Gauss, domains, rigging/refinement nonzero normalization, a true boundary-time physical source functional, connected TT and massless pole. Neither nullity F nor the cheap K1 certificate alone proves these.
4. No actual BQG Lorentzian amplitude source was supplied in this run; all actual PHYSICAL gates remain open.

Executable proof: \`research/hda_minimal_gate_2026_10_12/schur_tt_projector_symbolic_gate.py\`. The larger downloaded research package also includes a streaming Schur gateway, original exact SymPy controls, JSON and detailed report.
