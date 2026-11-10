# Formula ledger — article 1046730

Every formula in the article is to receive a stable ID. This is the controlling checklist; no formula is silently skipped. Status vocabulary: `NOT_STARTED`, `TRANSCRIBED`, `RECOMPUTED`, `PROVED`, `REFUTED`, `CONDITIONAL`, `UNRESOLVED`.

| ID | Article location / expression | Audit question | Status |
|---|---|---|---|
| G-01 | §Graph model: (pK=N^{-1/3}), (p=K^{-1}N^{-1/3}) | Is the scaling derived from graph ensemble/dynamics or assumed? Check dimensions, expected shortcut count, connectivity and limits. | TRANSCRIBED |
| G-02 | (K=2d_s=6) | Is spectral dimension equal to local coordination/2? Test on counterexample graphs. | TRANSCRIBED |
| G-03 | (R(N)=\pi[1/3-(K-\ln K)/\ln N]) | Definition vs derived order parameter; derive critical point. | TRANSCRIBED |
| G-04 | (G(N)=(R(N)-1)^2) | Tautological minimization and uniqueness/domain. | TRANSCRIBED |
| G-05 | (ln N_*= (K-\ln K)/(1/3-1/\pi)) | Recompute numerical value and sensitivity to K/dimension. | TRANSCRIBED |
| Z-01 | (zeta_G(s)=\sum_{i=1}^{N-1}\lambda_i^{-s}) | Definition; zero mode, connectedness, normalization. | TRANSCRIBED |
| Z-02 | (K(t)=\mathrm{Tr}\,e^{-t\Delta}), (K(0)=N) | Finite graph heat trace and Mellin transform domain. | TRANSCRIBED |
| Z-03 | (ln N_\zeta=K^{3/2+\zeta(2)}) | Derive from spectrum or refute. | TRANSCRIBED |
| Z-04 | (lambda_n=(2\pi/L)^2|n|^2) | Boundary conditions, torus normalization, finite volume. | TRANSCRIBED |
| Z-05 | (sum_{n\ne0}|n|^{-2}=K\zeta(2)) in 3D | Convergence; shell-counting and cutoff dependence. | TRANSCRIBED |
| Z-06 | (ho(\lambda)\sim\lambda^{1/2}), (mathcal N\sim(2K)^{3/2}) | Weyl-law normalization and relation to graph eigenvalue cutoff. | TRANSCRIBED |
| Z-07 | (zeta_G(0)_{reg}\sim N(2K)^{3/2}) | Meaning of regularization; compare finite-graph identity. | TRANSCRIBED |
| Z-08 | (p_n=p/n^2), (zeta_{nl}(1)\sim p\ln K\zeta(4)) then (p\ln K\zeta(2)) | Recompute series; claimed replacement of zeta(4) by zeta(2); graph ensemble meaning. | TRANSCRIBED |
| Z-09 | (exp[\zeta(2)\ln K]=K^{\zeta(2)}) | Algebraically true, but check whether preceding exponent is justified. | TRANSCRIBED |
| Z-10 | (ln N_\zeta=K^{3/2}K^{\zeta(2)}) | Product/power algebra and interpretation. | TRANSCRIBED |
| J-01 | (j_{uv}=\mathrm{Im}(\Psi_u^*S_{uv}\Psi_v)) | Scalar/spinor consistency, orientation, Hermiticity, units. | TRANSCRIBED |
| J-02 | (sum_{e}j_e\to J^\mu(x)) | Coarse-graining map, scaling and tensorial transformation. | TRANSCRIBED |
| J-03 | (J^\mu=\Psi^\dagger\Gamma^\mu\Psi) | Derivation of (Gamma^\mu), representation and current conservation. | TRANSCRIBED |
| J-04 | (partial_\mu J^\mu=0) and (partial_t\rho+\nabla\cdot J=0) | Derive from discrete unitary dynamics; specify density and edge divergence. | TRANSCRIBED |
| J-05 | Integral continuity equation | Boundary orientation, measure and assumptions. | TRANSCRIBED |
| J-06 | (g_{\mu\nu}\sim\langle\partial_\mu J^\alpha\partial_\nu J_\alpha\rangle) | Tensor signature, nondegeneracy, covariance, units, dynamics. | TRANSCRIBED |
| J-07 | (L_{eff}\sim N^{1/3}) vs (L_{eff}\sim\ln N) | Graph-distance theorem and regime of validity. | TRANSCRIBED |
| J-08 | Madelung velocity and quantum Hamilton–Jacobi / Q | Derive from Schrödinger equation; sign and coefficient of Q. | TRANSCRIBED |
| T-01 | Fixed point (Psi^*=\mathcal T[Psi^*]) unique/stable | Existence, uniqueness, contraction/compactness hypotheses. | TRANSCRIBED |
| T-02 | (Psi_v(t+t_P)=\mathcal T[\Psi_v(t)]) | Well-defined state space and dimensions. | TRANSCRIBED |
| T-03 | (\mathcal T=\mathcal T_{loc}+\mathcal T_{nonloc}+\mathcal T_{meta}) | Operator domains, linearity/nonlinearity, normalization. | TRANSCRIBED |
| T-04 | Local neighbor sum | Adjoint pairing and unitary/contractive conditions. | TRANSCRIBED |
| T-05 | Nonlocal sum plus memory kernel | Kernel causality, summability, stability and normalization. | TRANSCRIBED |
| T-06 | Meta term (-\eta\sum\delta\mathcal M/\delta\Psi^*) | Variational derivative convention and units. | TRANSCRIBED |
| T-07 | Node current balance (sum j=0) | Is this a dynamical consequence or an imposed constraint? | TRANSCRIBED |
| T-08 | (mathcal T[\Psi]=\Psi+t_P\mathcal D(J)+O(t_P^2)) | Continuum expansion and convergence. | TRANSCRIBED |
| T-09 | (|\mathcal T[\Psi]|=|\Psi|), norm (=N) | Define norm and prove unitarity. | TRANSCRIBED |
| S-01 | Random graph axioms / local (mathbb Z^3) / nonlocal (p(r)) | Fully specify graph probability law and thermodynamic limit. | TRANSCRIBED |
| S-02 | (mathcal H_N\cong\mathbb C^{2N}) | Tensor-product vs direct-sum Hilbert space and many-body dimension. | TRANSCRIBED |
| S-03 | (S_{uv}\in SL(2,\mathbb C)) | Compatibility with unitary quantum evolution; link orientation. | TRANSCRIBED |
| S-04 | (U_N=\exp(-iH_N)) | Missing (hbar), time interval, Hermiticity of H. | TRANSCRIBED |
| S-05 | Graph Hamiltonian edge sum | Hermiticity, double counting, dimensions, gauge covariance. | TRANSCRIBED |
| M-01 | Particle masses from (sqrt{\lambda_i}) | Units and physical mass operator; eigenmode identification. | TRANSCRIBED |
| M-02 | Particle-specific exponent tuples / algebraic basis | Count choices and demonstrate derivation before seeing data. | TRANSCRIBED |
| C-01 | All fundamental constants and particle masses/lifetimes | Reproduce each individual formula from source code; units and dependencies. | IN_PROGRESS |
| A-01 | 55×6 exponent matrix, rank and nullspace identities | Recompute exact rational rank/kernel; distinguish identities from physical laws. | IN_PROGRESS |
| V-01 | Free-energy functional and optimum (x_*) | Recompute all terms; parameter provenance; global vs local minimum. | TRANSCRIBED |
| V-02 | (k=xN^{2/3}), (p=k/N), (mu=\ln k/\ln N) | Check built-in exponent and definitions. | TRANSCRIBED |
| V-03 | (L=\ln N/\ln(K+\lambda k)) | Derive graph-distance scaling and check positivity/domain. | TRANSCRIBED |
| V-04 | Second derivative / claimed critical curvature | Correct nonuniform-grid derivative and stationarity. | TRANSCRIBED |

The ledger is intentionally not a summary of selected headline equations. The complete transcription must continue through every displayed formula in the article, including each constant/mass/lifetime formula and each formula embedded in images. New IDs are appended, never reused.
