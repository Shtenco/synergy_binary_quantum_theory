# BQG parity-protected q2 32D obstruction and next computation — 2026-10-12

**Status**: exact grading theorem using the saved numerical SU(2) Euclidean gap; actual full Lorentzian quantum HDA, physical projector and TT pole remain open.

Sources: main PETER_WEYL_LORENTZIAN_PARITY_FRONTIER.md; main research/audit_2026_10_10/FIVE_ROUTES_AUDIT.md; main research/audit_2026_10_10/fine_H0_full_images.npz.

For frozen Peter–Weyl grading Pi, Pi E Pi=−E, Pi L Pi=L. On the parity-pure original 32-dimensional code P32 and for any real lambda, outgoing images E_0 u and L_0 u are orthogonal. Thus:

\[
\|(E_0+\lambda L_0)u\|^2=\|E_0u\|^2+\lambda^2\|L_0u\|^2
\ge 3.97703147684868\|u\|^2.
\]

The derived norm lower bound is 1.9942496024441574 in the saved regulator normalization. Therefore the single node full normal constraint cannot annihilate any *nonzero source supported purely in P32* under the declared operator parity and numeric validity of stored H0. This does **NOT** exclude a full zero eigenvector with P32 overlap and additional components outside that subspace. The historical LDL certificate 3.89999968 applies to a rounded Gram and alone does not bound SU(2) recoupling error; if the raw Euclidean image error is epsilon the robust norm bound becomes max(0, 1.99424960244-epsilon).

A rigorous two-dimensional control makes the distinction explicit:
\[
E=a\sigma_x,\qquad L=aI,\qquad H=E+L,\quad a=\sqrt{3.97703147684868}.
\]
The pure code basis vector has nonzero full constraint norm squared 7.95406295369736, but the *full* normal constraint has zero vector (|even>-|odd>)/sqrt(2) with code weight 1/2.

On the complete parity even/odd carrier a Hermitian H has block form
\[
H=\begin{pmatrix}\lambda L_+&A\\A^\dagger&\lambda L_-\end{pmatrix}.
\]
For lambda nonzero and invertible L_-, the exact even-sector Schur zero equation is
\[
(\lambda L_+-\lambda^{-1}A L_-^{-1}A^\dagger)u=0.
\]
This is a **single-constraint** reduction, not proof of simultaneous five-node HDA or the physical projector.

For the vertical stack of all independent constraints C, with columns CP and CQ, the zero-energy Schur form is
\[
F=(CP)^\dagger(I-P_{\mathrm{ran}(CQ)})(CP)\succeq0.
\]
If sigma_min(CP)>=alpha and the minimum principal output-range angle between CP and CQ is theta_min, then the exact spectral bound is
\[
\lambda_{\min} F\ge\alpha^2\sin^2\theta_{\min}.
\]
A positive certified angle excludes an exact physical zero supported on the considered retained component; zero angle is merely a necessary possibility for a dressed physical zero.

**Automatic next computation:** stop searching for a pure P32 kernel. Use measured rank-30 collective q4 five-block metric carrier, freeze actual complete Hermitian H_E+lambda H_L ordering and the independent D, generate *all* outgoing images on the same full habitat for two declared collective TT probes (or all 30 if needed). Compute Schur principal angles with adaptively expanded Q, maintain complete support/leakage and test all 10 independent HDA pairs. Only afterwards is a normalized physical projector/history and C0/C1/C2 TT extraction meaningful.

Executable: parity_gap_source_gate.py; with --npz argument on cloned repo it independently computes 32-column SU(2) Gram from committed compressed CSR archive. The direct binary replay was not performed in this chat, so 3.97703147684868 is an upstream recorded result rather than a claimed fresh value.
