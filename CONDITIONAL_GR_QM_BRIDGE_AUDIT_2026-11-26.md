# Conditional GR–QM bridge audit — 2026-11-26

This note isolates the valid implication chain from the audited binary/HDA, projector, and TT results. The assumptions are part of the theorem statement; they are not silently promoted to microscopic BQG conclusions.

## 1. Dimension inputs are independent

For a (q)-level carrier, the traceless Hermitian space has real dimension (D_q=q^2-1). If—and only if—the full adjoint carrier is identified with an isotropic spatial vector space whose frame symmetry is the full (SO(D_q)), equating the group dimensions gives the nontrivial solution (D_q=3,q=2). For (q=2), the adjoint action indeed realizes (SU(2)/\mathbb Z_2\simeq SO(3)).

Separately, the frozen graph-growth recurrence gives (d_g=\log_2(N_g/N_{g-1})\to3). These are two compatible finite/kinematic indicators of three dimensions, but no theorem in these results identifies the graph-growth exponent with the isotropic adjoint-space dimension or with physical spectral dimension.

## 2. HDA fixes the local ADM parameters

For the local ADM family

\[
H[N]=\int d^3x\,N\left[A\frac{\pi_{ab}\pi^{ab}-c\pi^2}{\sqrt q}-B\sqrt q(R-2\Lambda)\right],
\]

the documented classical bracket is

\[
\{H[N],H[M]\}=AB\left(D[\beta]+4(c-1/2)I[N,M]\right).
\]

Matching the standard hypersurface-deformation algebra for generic canonical data requires (c=1/2) and (AB=1). The ADM parameter gate was rerun at (L=7) and returned `passed: true`: the (c=1/2,AB=1) cases had standard-HDA residual about (4.65\times10^{-9}); the (c=0.4) case had defect (2.84594\times10^{-2}); cases with (AB=2) had defect (1/3). Varying (Lambda) from −10 to 10 left the bracket at the same spectral-error scale.

This closes parameter selection only inside the local ADM ansatz. It leaves the absolute Newton/canonical normalization (A/B) and (Lambda) free. It does not derive the microscopic BQG constraint family, its anomaly-free HDA, or these coefficients from a physical state.

## 3. From HDA to a graviton is conditional

If the spatial metric is the complete local configuration field and the (D+1) diffeomorphism constraints are independent, regular, and first class, then

\[
N_{\rm phys}^{\rm config}=\frac{D(D+1)}2-(D+1)=\frac{(D+1)(D-2)}2.
\]

At (D=3), this is two configuration modes. Identifying them with the two helicities of a massless spin-2 field additionally needs the Lorentzian gauge reduction and the correct tensor representation. A mode count alone does not prove a stable pole or positive norm.

If, further, the physical IR theory is a local, diffeomorphism-invariant, metric-only action with at most two derivatives, its leading action is the Einstein–Hilbert term plus a cosmological term, up to normalization and boundary terms. Neither the HDA bracket nor the degree count alone proves those action hypotheses. In particular, this conditional uniqueness statement should not be advertised as a replacement for the full Hojman–Kuchař–Teitelboim theorem.

## 4. Quantum physicalization is a separate bridge

For fixed finite constraints (C_A) and positive constraint metric (G),

\[
M_G=C^\dagger(G\otimes I)C\ge0,
\qquad\ker M_G=\bigcap_A\ker C_A.
\]

The finite master-projector gate was rerun and passed. This gives an exact finite-dimensional common-kernel theorem. To define a quantum-graviton observable one still needs a nontrivial, regulator-stable physical inner product/history; a source functional (Z[J]); connected correlators and a physical 1PI kernel; and then a TT pole with positive residue and the required gauge-reduced helicities. These objects do not follow from the finite projector identity.

The reduced TT lattice transfer separately has (c_{\rm red}=1/\sqrt3) in its own lattice units and a bare isotropic quartic coefficient −1/45. The repository's map (c_{TT}=v_0a_*/\tau_*) converts units algebraically; it does not derive (a_*/\tau_*), establish a physical photon kernel, or carry the reduced pole through the constrained-history RG limit.

## 5. Shared photon/graviton cone and information mode

With one emergent Lorentz metric and a physical Maxwell sector using the same metric/Hodge structure, both massless inverse kernels depend at leading order on the same invariant (s=-\omega^2+c_{\rm metric}^2k^2). Their leading cones therefore coincide under those assumptions. The physical common metric and Maxwell sector remain open inputs in the current BQG ledgers.

The branch's Pancharatnam gate was executed from its fetched source. Link-covariance error was \(2.48\times10^{-16}\), holonomy-invariance error zero, Bargmann-holonomy error \(2.22\times10^{-16}\), and the numerical Hopf Chern number was \(1.0000000000000013\). The link/holonomy checks use the sampled spinors; the Chern integration directly integrates the analytic \(F=\tfrac12\sin\theta\,d\theta\wedge d\phi\), rather than estimating curvature from plaquette data. This verifies the finite ray/phase identities and checks the known analytic Chern integral, but those are separate controls. Its canonical-Maxwell gate also passed: on the finite 16-cell, \(d_1d_0=0\), \(\mathrm{rank}\,d_0=7\), \(\mathrm{rank}\,d_1=17\), \(b_1=0\), and the unit-Hodge curl spectrum was \(0^{\times7},4^{\times6},6^{\times8},8^{\times3}\). However, the gate's reported \(Z_A\) cancellation is implemented by simplifying \((1/Z_A)Z_A=1\) after assuming the canonical equations; it does not derive those equations from BQG or compute a physical stiffness. The source explicitly leaves deconfinement, final Hodge weights, continuum photon poles, and \(Z_A\) open. A small-holonomy Wilson action expands to a quadratic curvature term only conditional on a positive blocked stiffness and geometric Hodge weights. The branch cites Polyakov's compact electrodynamics result in 2+1 dimensions and Guth's weak-coupling nonconfining result for four-dimensional Villain U(1) as dimensional compatibility evidence; neither selects the Coulomb basin in BQG.

The branch's symbolic QGT gate was also executed and passed. For the stated Bloch spinor it returns
\[
Q_{ab}=\begin{pmatrix}\tfrac14&i\sin\theta/4\\-i\sin\theta/4&\sin^2\theta/4\end{pmatrix},
\quad
\operatorname{Re}Q=g_{FS},
\quad
F_{\theta\phi}=2\operatorname{Im}Q_{\theta\phi}=\tfrac12\sin\theta,
\quad
\frac1{2\pi}\int_{S^2}F=1.
\]
This is an exact common ray-space origin for distinguishability geometry and compact Berry phase. It does not provide a derived map from \(g_{FS}\) to the spacetime metric, nor does it make the Berry connection a propagating photon.

The current T75 source/artifact was also rechecked by rerunning the joint-quotient generator with seed 20261118 and four samples per fixture. The path\(_3\) fixture has \(12V=36\), \(3V+3E=15\), and quotient dimension 21; cycle\(_4\) has \(12V=48\), \(3V+3E=24\), and quotient dimension 24. Independence of the stated Gauss and edge-singlet generators is proved by the fixture's unpaired boundary half-links. The rerun found no additional QFI nulls in eight samples; minimum sampled eigenvalues were 0.04557 and 0.01617, with maximum joint annihilation residual \\(1.25\\times10^{-15}\\). Those eigenvalue values and residuals remain floating sample diagnostics without interval certificates. The subsequent exact minor certificates in Section 7 prove generic absence of additional tangent nulls for these two fixed fixtures, but do not give an interval spectral gap or uniform refinement lower bound. This QFI quotient is not yet the spacetime metric: a physical geometric map and a certified lower singular bound are still needed.

Two separate finite controls supply adjacent pieces, but their composition is not established. The local logical-shape gate maps a two-coordinate qubit singlet shape \((X,Z)\) to a rank-two, trace-free tangent of one reconstructed tetrahedron metric. The C8 source gate separately verifies a rank-two relational projector control and computes \(Z(j_X,j_Z)=\cosh\sqrt{j_X^2+j_Z^2}\), \(C_{metric}=BB^T\), with Frobenius-component response eigenvalue \(9/2\) and tangent pseudoinverse eigenvalue \(2/9\). Those numerical response eigenvalues depend on the declared source/component normalization. The C8 system step \(R=J\) and the T75 random-intertwiner path/cycle quotient are not shown to be the same physical Hilbert space or the same history. No refinement-stable map from the T75 quotient to the local shape doublet is supplied, so the gates cannot yet be multiplied into a physical spacetime metric correlator.

An infoton candidate requires a gauge-invariant physical operator and positive physical spectral pole. If it exists in the same history as the graviton, a common retarded inverse kernel has the block form

\[
\mathbb K^R=\begin{pmatrix}K_T^R&K_{TI}^R\\K_{IT}^R&K_I^R\end{pmatrix},
\qquad\det\mathbb K^R=0.
\]

The determinant locates candidate normal modes only; positive spectral weights and regulator-stable residues are also required. (K_{TI}) and the (hII) vertex must come from the same physical effective action. The current finite Floquet example has externally prescribed modulation and null derived coupling, so it establishes no physical infoton or interaction.

## 6. Combined conditional theorem

The defensible one-pass statement is:

\[
\begin{aligned}
&\text{binary carrier + full isotropy}\Rightarrow D=3,\\
&\text{local ADM ansatz + standard first-class HDA}\Rightarrow c=1/2,\ AB=1,\\
&\text{regular diffeomorphism constraints}\Rightarrow 2\text{ configuration modes},\\
&\text{local metric-only two-derivative IR hypotheses}\Rightarrow\text{EH leading dynamics},\\
&\text{physical projector/history + positive TT 1PI pole}\Rightarrow\text{quantized massless spin-2 excitations},\\
&\text{same Lorentz metric + physical Maxwell sector}\Rightarrow c_g=c_\gamma,\\
&\text{physical }Q_I+\text{positive pole + common }\Gamma\Rightarrow\text{candidate infoton and derived coupling}.
\end{aligned}
\]

The first lines are conditional kinematic/classical theorems and finite controls. The decisive microscopic antecedents in the final lines—physical history, TT pole, Maxwell kernel, and infoton operator/pole—are not established by the current evidence. Thus this is a valid conditional bridge map, not a completed proof of quantum binary gravity.

## 7. Constructive finite T75 rank witnesses

A targeted search found explicit local intertwiner coefficient vectors with full numerical rank in the two audited fixtures. In the local logical basis \((L_0,L_1)\), choose the unnormalized node coefficients

\[
\begin{array}{c|c}
\text{path}_3 &(1,-1+2i),\;(-1+2i,2-2i),\;(-1,1-i)\\
\text{cycle}_4 &(i,-2i),\;(1,2i),\;(2,-1-i),\;(1,1+i)
\end{array}
\]

The numerical gate normalized each node before contraction; the exact minor uses the listed vectors as unnormalized ray representatives, which changes a common scale but not rank. The floating implementation returned quotient dimensions 21 and 24, ranks 21 and 24, and smallest singular values 0.212005 and 0.0584274 respectively (rank threshold (10^{-9})). Thus the earlier absence-of-extra-null result is no longer supported only by random sampling: these are reproducible full-rank finite witnesses with explicit coefficients.

The exact certificate is reproducible with `python scripts/q2_relational_exact_tangent_minor_certificate.py`; it reconstructs the stated coefficients, logical singlets, normalized edge epsilon, and projective tangent formula. The resulting determinants in \\(\\mathbb Q(\\sqrt2,\\sqrt3,i)\\) are

\\[
\\det M_{path_3}=\\frac{10187584608944228345425+5881804716378939721075\\sqrt3}{231812806445087701493115518976}>0,
\\]

\\[
\\det M_{cycle_4}=-\\frac{688337770980009774477524933597}{832230283493566508763382283981944794232021271445504}-\\frac{33117666447391300525114547657\\sqrt3}{69352523624463875730281856998495399519335105953792}<0.
\\]

These are exact algebraic nonzero-minor certificates. The tangent used for the exact determinant is \\(d\\phi-\\phi\\langle\\phi,d\\phi\\rangle/\\langle\\phi,\\phi\\rangle\\), which differs from the normalized-ray tangent only by a nonzero common scale; the node vectors in the table may likewise be normalized without changing rank. Together with the exact independent Gauss and edge-singlet kernel directions (15 for \\(path_3\\), 24 for \\(cycle_4\\)), they prove quotient ranks 21 and 24 at these witnesses. The projective tangent entries are rational real-analytic functions of the node-ray coordinates wherever the postselection norm is nonzero. Clearing powers of that norm from each minor gives a real-analytic numerator on the connected product of node intertwiner rays. Since that numerator is nonzero at the exact witness, its zero set has empty interior; its nonzero set is open dense. The exact symmetry kernel count supplies the matching rank upper bound, so quotient rank is generically 21 and 24 on the two fixed graph fixtures. Exceptional states can still have extra stabilizers.

This exact finite generic-rank result is not a uniform lower singular-value bound under graph refinement, a proof for other graphs, or a derivation of a spacetime metric. The T75 result remains a finite two-graph theorem/control and does not by itself establish the continuum BQG metric.


For reproducibility, pivoted QR on the realified full tangent first selected these zero-based row/column indices. The exact script freezes the indices and verifies the corresponding determinants:

| Fixture | Rows | Columns | Numerical minor |
|---|---|---|---|
| `path_3` | 307,85,52,45,203,476,492,75,459,311,362,59,277,55,179,396,54,42,213,317,314 | 0,1,2,30,31,32,18,19,20,4,3,5,33,34,35,6,13,8,23,22,15 | sign +; log10(abs(det)) = −19.28954; condition ≈ 15.73 |
| `cycle_4` | 102,331,347,43,294,23,103,152,356,411,286,70,145,401,326,480,283,155,11,26,357,281,282,89 | 3,4,5,13,12,14,34,33,35,9,10,11,45,47,46,24,25,26,37,21,23,15,16,17 | sign −; log10(abs(det)) = −23.74258; condition ≈ 64.63 |

The row convention is first all real boundary-amplitude coordinates, then all imaginary coordinates; tangent columns follow node, leg, and Pauli-axis order in `projective_tangent_jacobian`. QR originally selected the indices in floating point, but the exact script now freezes and verifies the listed minor determinants independently; a different BLAS implementation may select other indices without affecting this existing certificate.


## 8. Exact local quantum-ray to spatial-shape metric

The four-spin singlet logical qubit has intrinsic shape Bloch coordinates \\(X,Z\\); its oriented-volume coordinate is \\(Y=\\pm\\sqrt{1-X^2-Z^2}\\). The pure-qubit Fubini–Study line element is one quarter of the unit Bloch-sphere line element. At either regular tetrahedral branch \\(X=Z=0, Y=\\pm1\\),

\\[
 ds^2_{FS}=\\frac14(dX^2+dZ^2).
\\]

The exact local reconstruction gate gives the background edge metric \\(g_0\\) and trace-free metric derivatives \\(M_X=\\partial_Xg|_0\\), \\(M_Z=\\partial_Zg|_0\\). Its exact DeWitt-type tangent inner product is

\\[
\\operatorname{Tr}(g_0^{-1}M_Ag_0^{-1}M_B)=\\frac32\\delta_{AB},
\\qquad A,B\\in\\{X,Z\\}.
\\]

Therefore, at the isotropic tetrahedron and with these stated normalizations, the pullback of the spatial shape tangent inner product is exactly \\(6g_{FS}\\). This supplies a local, analytic identification between distinguishability of logical quantum rays and intrinsic tetrahedral metric-shape response. The factor is convention-dependent, and the identity is a tangent-space result at the regular shape. It does not construct a graph-to-continuum metric field, physical time evolution, an Einstein kinetic operator, a graviton pole, or a photon/infoton sector. The missing T75-to-shape map remains explicit: the T75 frame coordinates and node-shape coordinates remain distinct input sectors. Section 9 tests them together through the same finite boundary contraction, but this does not identify a physical projection from frame quotient to shape or derive a continuum history map.


## 9. Combined finite frame-and-shape response

Extend the same graph boundary map by allowing both half-link frame variations and two real projective variations of each node intertwiner ray:

\\[
(\\delta\\omega_{frame},\\delta X_1,\\delta Z_1,\\ldots,\\delta X_V,\\delta Z_V)
\\longmapsto\\delta[\\mathcal C_\\Gamma(\\otimes_v|\\iota_v\\rangle)].
\\]

For each logical ray \\(c=(a,b)\\), the orthogonal tangent vector \\(u=(-\\bar b,\\bar a)\\) and \\(iu\\) span its two real ray directions. The shape derivative is computed by replacing the corresponding node tensor in the multilinear contraction, then projecting off the boundary ray. For the path_3 exact T75 witness and the cycle_4 nondegenerate witness recorded in scripts/q2_relational_frame_shape_response_gate.py, every node has nonzero oriented-volume Bloch coordinate \\(Y\\), so \\(X,Z\\) are valid local shape coordinates. The response gate first found full rank numerically; an exact residue certificate now upgrades the path case, while the cycle case remains numerical:

| Fixture | Frame quotient | Node-shape directions | Combined rank / dimension | Certificate |
|---|---:|---:|---:|---|
| \\(path_3\\) | 21 | 6 | 27 / 27 | Exact modular minor: determinant residue 53 mod 97 |
| \\(cycle_4\\) | 24 | 8 | 32 / 32 | Floating SVD; minimum singular value 0.292778 |

For \\(path_3\\), the exact modular certificate proves that the combined boundary-ray differential has rank 27 at the listed witness. The exact frame kernel has dimension 15, so this is the maximal possible rank of its 42 frame-plus-shape input coordinates. The determinant is a real-analytic minor after clearing the postselection-norm denominator; being nonzero at one point, it proves maximal rank on an open dense subset of this fixed path graph's node-ray space. Since all three node volume coordinates are nonzero, \\(X,Z\\) are valid shape coordinates there, and Section 8 maps each such tangent with rank two to the local tetrahedral metric shape.

For \\(cycle_4\\), the same exact multilinear shape derivative and quotient frame calculation give full floating rank at the displayed nondegenerate witness. Its 32×32 exact modular determinant has not yet been certified. These are finite graph results and give no graph-refinement lower bound. Neither result constructs a physical Hamiltonian history, time metric, continuum spacetime field, or graviton, photon, or infoton pole. Frame and shape remain distinct input sectors even where their boundary responses are jointly injective.
