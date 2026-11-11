# Exact BQG TT source obstruction and collective remedy — 2026-10-12

Scientific scope **PROVED KINEMATIC; PHYSICAL TT/HDA OPEN**. Sources: main \`Q2_RELATIONAL_METRIC_SOURCE_GENERATING_FUNCTIONAL.md\` and the measured collective \`M_hq\` in branch \`research/bcqg-core-candidate-v1/BCQG_IR_HDA_UNIVERSALITY_CLOSURE.md\`. No microscopic physical-history data fabricated.

## Exact shape-source Gram

For traceless symmetric sources \(M_X,M_Z\) already frozen in the q2 relational bridge, use \(P_n=I-nn^T/q,\;q=x^2+y^2+z^2\) and
\(\mathrm{TT}_n(M)=P_nMP_n-\tfrac12P_n\mathrm{tr}(P_nM)\).
The 2x2 Frobenius Gram of the projected sources has determinant

\[
\boxed{\det G_{XZ}^{TT}(n)=\frac{3}{4q^3}F(x,y,z)^2,}
\]

\[
F=2(x^3+y^3+z^3)-x^2y-x^2z-xy^2-xz^2-y^2z-yz^2-3xyz.
\]

This is an exact SymPy symbolic factorization using the **actual BQG q2 shape-to-metric matrices**. At direction \((1,-1,0)\): \(TT(M_X)=0\), \(G=\operatorname{diag}(0,33/8)\): rank ONE. At \((0,0,1)\): \(G=\operatorname{diag}(3/2,2)\): rank TWO. At \((1,1,2)\), also rank ONE. Thus the original fixed **two** local metric shapes are not a globally complete **two-polarization TT source frame**.

## Rigorous all-direction six-dimensional repair

For the already measured \(\mathrm{rank}(M_{hq})=6\) collective metric map with squared singular values \(\{1/6\;(3\times),\,1/3\;(3\times)\}\), take any 2x6 orthonormal TT frame \(E_n E_n^\dagger=I_2\), and \(S_n=E_n M_{hq}\). By Rayleigh quotient:

\[
\boxed{\frac16 I_2\preceq S_nS_n^\dagger\preceq\frac13 I_2\quad\forall n\ne0.}
\]

Thus TWO kinematic TT sources survive for all momenta; the minimum determinant is \(\ge1/36\) in the declared metric normalization. Exact canonical A1/E/T2 numerical directions: n=(001) -> eigen {1/6,1/3}; (111) -> {5/18,5/18}; (1,-1,0) -> {5/24,1/3}. This **does not** imply two *physical Dirac* sources: the full common Lorentzian projector could kill both. The five-block \(960\times30\) real P carrier in \`physical-tt-kernel-leading-ir/BQG_MULTIBLOCK_TT_WILSON_PRODUCTION_STATUS.md\` is **kinematic**, not the physical master zero-sector.

## A second cheap blocker: Lorentzian local logical Y

Branch \`research/bcqg-core-candidate-v1/PETER_WEYL_LORENTZIAN_ENVTRACE_RESULT.md\` finds \(\;L_{raw,1body}=i\,cY\), \(c=1.3389293521464034\), with 16 diagonal environments traced. For its **conditional** odd Hermitian completion \(H_{odd}=cY\):

\[
[Y,X]=-2iZ,\quad [Y,Z]=2iX,\quad
\|[H_{odd},X]\|_{op}=\|[H_{odd},Z]\|_{op}=2c=2.677858704292807.
\]

So bare X/Z are NOT invariants of this isolated logical generator. Its 2x2 master is \(c^2I\) and has no zero sector if *mistakenly* treated as the entire scalar constraint. But this averaged block is NOT the full graph-changing \(H_E+\lambda_LH_L\); none of these local statements is a full physical no-go.

## Physical source survival certificate

For a genuine finite orthogonal physical projector \(P_{phys}\), source-excitation map \(U:\mathbb C^2\to\mathcal H\), and centered operators defined on the same physical vacuum:

\[
G^{phys}_0=U^\dagger P_{phys}U,\quad
\lambda_{\min}(G^{phys}_0)\ge
\lambda_{\min}(U^\dagger U)-\|(I-P_{phys})U\|^2.
\]

This is a rigorous sufficient rank-two condition, but the q4 kinematic \(1/6\) lower bound cannot be transplanted into \(\lambda_{\min}(U^\dagger U)\) unless the state-source embedding has independently been proven isometric.

## Independent cheap zero-kernel Schur screen

The collective source branch reports \(PGP=0\). On any selfadjoint graph-changing complete operator \(G=\left(\begin{smallmatrix}0&B\\B^\dagger&D\end{smallmatrix}\right)\) with invertible \(D\), the exact Schur zero-mode condition is
\(C_{eff}(0)=-BD^{-1}B^\dagger\).
If \(D\) is sign-definite and \(B^\dagger\) injective, no exact zero state lies in retained carrier. That is a necessary *screen* for an exact zero mode, not a measured result: the actual microscopic B,D remain uncomputed.

## Next automatic step

**Result:** 2 local X/Z source-only route **FAIL** (directional rank). Six-dimensional collective q4 metric route **PASS** (kinematic all-direction rank). Bare Dirac invariance **NOT CERTIFIED**; actual physical P/eta and boundary transfer **BLOCKED**. Therefore the next calculation is **not** another free TT or Wilson fit: reuse real five-block P and frozen hermitian E/S to compute complete action on its six q-directions; construct full common-domain constraint master and independent D; determine projected centered TT source Gram and physical history, then real C0,C1,C2. Never relabel constraint heat time as physical time.

Run \`python research/hda_minimal_gate_2026_10_12/tt_source_rank_frozen_gate.py\`. This file intentionally has no full-HDA/physical-TT PASS labels.
