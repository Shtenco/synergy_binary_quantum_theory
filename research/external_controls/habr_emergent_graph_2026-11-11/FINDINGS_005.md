# Findings 005 — exact finite-graph spectral identities and TSCO/current consistency

## 1. Finite-graph spectral zeta identities

For a finite connected graph with combinatorial Laplacian Δ and eigenvalues 0=λ0<λ1≤...≤λN−1, define the nonzero-spectrum zeta function exactly as the article does:

    ζ_Δ(s) = Σ_{i=1}^{N−1} λ_i^(−s).

Then at s=0 no analytic continuation is needed because the sum is finite:

    ζ_Δ(0) = Σ_{i=1}^{N−1} 1 = N−1.

The article later states ζ_Δ(0) ~ −ln N. That is incompatible with this finite-graph definition for large N. A different regularized/infinite-volume zeta may have a different value, but it must be separately defined and cannot be silently identified with the finite-graph sum.

Differentiate the finite sum:

    ζ'_Δ(0) = −Σ_{i=1}^{N−1} ln λ_i.

Therefore:

    −ζ'_Δ(0) = ln(Π_{i=1}^{N−1} λ_i) = ln det'(Δ).

For a connected unweighted graph, the matrix-tree theorem gives:

    Π_{i=1}^{N−1} λ_i = N τ(G),

where τ(G) is the number of spanning trees. Hence:

    −ζ'_Δ(0) = ln N + ln τ(G),

not ln N in general. The article's equation lnN = −ζ'_Δ(0) would require τ(G)=1 (or a separately justified normalization subtracting ln τ(G)); this is not true for a generic connected graph and is not established for the proposed small-world graph.

## 2. Heat-trace statement is correct but does not imply the claimed N formula

For a finite graph, K(t)=Tr(exp(−tΔ)) satisfies K(0)=N because exp(0)=I. This elementary identity is correct. But the finite heat trace at t=0 gives N directly; it does not derive the article's expression lnNζ=K^(3/2+ζ(2)). Any route from a continuum Weyl law to a finite graph count must specify the cutoff, eigenvalue normalization, boundary conditions, and the order of limits.

## 3. The graph degree K=6 does not by itself derive spectral dimension 3

K=6 is the coordination number of the nearest-neighbor cubic lattice in three dimensions. But the identity K=2d is specific to a hypercubic nearest-neighbor lattice. A two-dimensional triangular lattice also has degree 6 while its large-scale spectral dimension is 2. Therefore degree alone does not determine spectral dimension. The article's local Z^3 assumption already inserts the three-dimensional lattice structure; it is not an emergent result of the subsequent K=6 substitution.

## 4. TSCO fixed-point theorem: hypotheses missing

The article asserts existence, uniqueness, and stability of Ψ*=T[Ψ*]. Those conclusions do not follow from writing a map T. A fixed-point theorem needs a specified complete state space, a well-defined map, and conditions such as a contraction (or another explicit existence/uniqueness theorem). A local derivative bound |T'|<1 is not enough unless the norm, domain, differentiability, and global/appropriate local hypotheses are specified.

More importantly, the displayed nonlocal TSCO term contains a history sum Σ_{τ<t} K(t−τ) Ψ_v(τ). It is therefore not a map of the instantaneous Ψ(t) alone. To formulate a fixed-point theorem, the state must be enlarged to include the required history (or the memory kernel must be eliminated by an equivalent Markovian auxiliary-state construction). The article does not do this.

## 5. Current conservation: local continuity is not zero divergence

The edge current is j_uv=Im(Ψ_u† S_uv Ψ_v). For an oriented edge, antisymmetry j_vu=−j_uv follows if the reverse link obeys S_vu=S_uv† (with consistent orientation and inner product). This condition is not stated as an axiom in the displayed construction.

For a time-dependent node density ρ_v=||Ψ_v||², the discrete continuity equation has the form:

    dρ_v/dt + Σ_{u~v} j_{vu} = 0

up to the chosen sign convention and inclusion of nonlocal edges. The article instead writes Σ j_{vu}=0 at each node. That is a divergence-free stationary-current constraint, not the general continuity equation. Global norm conservation under unitary evolution does not imply that the instantaneous current divergence vanishes at every node; it implies the summed/global continuity relation.

## 6. SL(2,C) links and unitary evolution require a compatibility condition

The article assigns S_uv∈SL(2,C) to links and separately asserts a unitary global update U=exp(−iH). Membership in SL(2,C) does not imply unitarity: SL(2,C) is noncompact and its generic elements are not unitary. To obtain a Hermitian H and unitary U, the link variables and reverse-edge pairing must satisfy explicit adjointness/inner-product conditions. A possible construction may exist, but it is not implied by S_uv∈SL(2,C) alone.

The article's U=exp(−iH) also suppresses ℏ and the time interval. In ordinary units, U(Δt)=exp(−iHΔt/ℏ); if natural units and a unit time step are intended, those conventions must be fixed.

## 7. Metric-from-current formula is not yet a metric construction

The article proposes g_{μν} ∼ <∂_μ J^α ∂_ν J_α>. This expression requires a definition of coordinates/derivatives, an index-lowering metric for J_α, a tensorial transformation law, and a proof that the resulting bilinear form is nondegenerate with Lorentzian signature. Using J_α already presupposes a metric or an explicitly defined internal bilinear form; this must be distinguished from the emergent spacetime metric to avoid circularity.

## 8. Mass formula changes within the article

The article first states m_i=√λ_i (spectral-mass section), later states m_n∝λ_n (spectral axioms section), and later returns to m∼√λ_m for composite modes. These are not the same scaling. The local dispersion λ_n∼|n|²/N^(2/3) implies √λ_n∼|n|/N^(1/3), whereas λ_n itself scales as |n|²/N^(2/3). The choice must be fixed and derived from a physical dispersion relation before matching masses.

## 9. Current verdict

- Finite-graph ζ(0)=N−1: exact theorem; article's ζ(0)~−lnN is incompatible with the stated finite sum.
- −ζ'(0)=ln det'(Δ)=ln[Nτ(G)]: exact; article's −ζ'(0)=lnN omits the spanning-tree factor.
- Heat trace K(0)=N: PASS, but it does not derive the claimed spectral formula.
- K=6⇒d_s=3: REFUTED as a general implication; only follows after choosing the cubic-lattice local structure.
- TSCO existence/uniqueness/stability: UNPROVED; state space and contraction hypotheses missing.
- Local current antisymmetry and continuity: CONDITIONAL; reverse-link adjointness and time-dependent continuity equation required.
- SL(2,C) link variables imply unitary update: FALSE without extra compatibility conditions.
- Metric-from-current: UNPROVED and potentially circular as written.
- Mass spectral scaling: INTERNAL NOTATION/SCALING INCONSISTENCY.