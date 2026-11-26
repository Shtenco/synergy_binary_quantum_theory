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