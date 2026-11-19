# Q2 continuum-limit research program (19 November 2026)

## Status and objective

The current G0--G3 computations establish finite-level combinatorial and spectral facts. They do **not** establish a continuum limit. In particular, the G3 heat-trace/effective-dimension data do not yet exhibit a demonstrated scale-independent plateau, and the current G3 Feshbach report evaluates four deterministic coarse probes rather than the complete 512 x 512 self-energy matrix.

The objective is to construct and test a mathematically specified sequence of measured spaces, Hilbert spaces, closed quadratic forms, and self-adjoint operators whose limit is independent of arbitrary discretization choices. A continuum claim is accepted only after analytic hypotheses and numerical falsification tests are both met.

## Literature anchors

1. K. Kuwae and T. Shioya, *Convergence of spectral structures: a functional analytic theory and its applications to spectral geometry*, Communications in Analysis and Geometry 11 (2003), 599--673. DOI: 10.4310/CAG.2003.v11.n4.a1. Framework includes Mosco convergence of quadratic forms, spectral structures, Riemannian manifolds, locally finite graphs, and metric spaces with Dirichlet forms.
   - https://www.intlpress.com/site/pub/files/_fulltext/journals/cag/2003/0011/0004/CAG-2003-0011-0004-a001.pdf
2. I. Babuška and J. E. Osborn, *Finite Element-Galerkin Approximation of the Eigenvalues and Eigenvectors of Selfadjoint Problems*, Mathematics of Computation 52 (1989), 275--297. DOI: 10.1090/S0025-5718-1989-0962210-8.
   - https://www.ams.org/journals/mcom/1989-52-186/S0025-5718-1989-0962210-8/
3. J. Kigami, *Analysis on Fractals* (Cambridge University Press, 2001), especially the construction of limiting Dirichlet/resistance forms from compatible finite networks.
   - https://doi.org/10.1017/CBO9780511470943
4. J. Ambjørn, D. N. Coumbe, J. Gizbert-Studnicki, J. Jurkiewicz, *Searching for a continuum limit in causal dynamical triangulation quantum gravity*, Physical Review D 93, 104032 (2016). DOI: 10.1103/PhysRevD.93.104032.
   - https://doi.org/10.1103/PhysRevD.93.104032
5. M. Reuter and F. Saueressig, *Quantum Einstein Gravity*, New Journal of Physics 14 (2012), 055022. DOI: 10.1088/1367-2630/14/5/055022.
   - https://inspirehep.net/literature/1088697
6. G. Dusson, I. M. Sigal, B. Stamm, *The Feshbach--Schur map and perturbation theory*, 2021; and *Analysis of the Feshbach--Schur method for the Fourier spectral discretizations of Schrödinger operators*, Mathematics of Computation 92 (2023), 217--249.
   - https://ems.press/books/ecr/187/3518
   - https://doi.org/10.1090/MCOM/3774
7. G. Elek, *Weak convergence of finite graphs, integrated density of states and a Cheeger type inequality*, Journal of Combinatorial Theory, Series B 98 (2008), 62--68. DOI: 10.1016/j.jctb.2007.03.004.
   - https://doi.org/10.1016/j.jctb.2007.03.004

8. J. Ambjørn, J. Gizbert-Studnicki, A. Görlich, D. Németh, *Is lattice quantum gravity asymptotically safe? Making contact between causal dynamical triangulations and the functional renormalization group*, Physical Review D 110, 126006 (2024). DOI: 10.1103/PhysRevD.110.126006. This is a particularly relevant current bridge between CDT finite-volume scaling and FRG fixed-point language; the authors explicitly state that the numerical precision allows a putative UV fixed point but does not prove its existence.
   - https://doi.org/10.1103/PhysRevD.110.126006

9. N. Gigli, A. Mondino, G. Savaré, *Convergence of pointed non-compact metric measure spaces and stability of Ricci curvature bounds and heat flows*, Proceedings of the London Mathematical Society 111 (2015), 1071--1129. DOI: 10.1112/plms/pdv047. This supplies a measured-space framework, Γ/Mosco convergence of Cheeger energies, heat-flow stability, and spectral convergence under additional hypotheses.
   - https://doi.org/10.1112/plms/pdv047
10. L. Ambrosio and S. Honda, *New stability results for sequences of metric measure spaces with uniform Ricci bounds from below* (2016). Relevant for extrinsic embeddings into a common metric space, Mosco convergence, Sobolev compactness, and stability of Cheeger constants.
   - https://arxiv.org/abs/1605.07908

11. N. García Trillos, M. Gerlach, M. Hein, D. Slepčev, *Error estimates for spectral convergence of the graph Laplacian on random geometric graphs towards the Laplace--Beltrami operator* (2019). Provides quantitative eigenvalue/eigenvector convergence rates for appropriately scaled geometric graphs sampled from manifolds. It is not directly the same deterministic PL-dual construction, but it is a valuable benchmark for what a graph-to-manifold theorem must quantify.
   - https://arxiv.org/abs/1801.10108
12. J. Calder and N. García Trillos, *Improved spectral convergence rates for graph Laplacians on ε-graphs and k-NN graphs*, Applied and Computational Harmonic Analysis 60 (2022), 123--175. Useful for explicit scaling regimes and convergence-rate diagnostics.
   - https://doi.org/10.1016/j.acha.2022.02.004
13. M. Belkin and P. Niyogi, *Convergence of Laplacian Eigenmaps* (2008). Shows convergence of suitably constructed graph-Laplacian eigenvectors to Laplace--Beltrami eigenfunctions under manifold sampling assumptions.
   - https://www.aoml.noaa.gov/ftp/pub/phod/pub/thacker/SSAmanifoldPapers/BelkinNiyogi2008.pdf

14. J. Dodziuk, *Finite-Difference Approach to the Hodge Theory of Harmonic Forms*, American Journal of Mathematics 98 (1976), 79--104. For suitable inner products on cochains of smooth Riemannian triangulations, this constructs a combinatorial Hodge theory approximating the smooth one; the paper also proves convergence results for eigenvalues of the Laplacian on functions under subdivision. This is especially relevant to a barycentrically refined PL complex.
   - https://doi.org/10.2307/2373615
15. M. Wardetzky, S. Mathur, F. Kälberer, E. Grinspun, *Discrete Laplace Operators: No Free Lunch*, Symposium on Geometry Processing (2007), DOI: 10.2312/SGP/SGP07/033-037. The paper proves that on general triangular meshes no single discrete Laplacian satisfies all the natural desiderata simultaneously. It is a warning to state explicitly which properties our geometric dual and primal FEM operators preserve; its exact theorem is for surface meshes, so extension to tetrahedral 3-complexes must be justified separately.
   - https://doi.org/10.2312/SGP/SGP07/033-037

9. N. García Trillos, M. Gerlach, M. Hein, D. Slepčev, *Error estimates for spectral convergence of the graph Laplacian on random geometric graphs towards the Laplace--Beltrami operator* (2018). arXiv:1801.10108. Gives quantitative eigenvalue/eigenvector convergence rates under sampling and bandwidth hypotheses, with transport-based comparison of empirical and continuum measures.
   - https://arxiv.org/abs/1801.10108
10. N. García Trillos and D. Slepčev, *Continuum limit of total variation on point clouds* (2016), and related work on graph Dirichlet energies and spectral convergence. The key transferable technique is the TL^p transport topology plus Gamma-convergence and compactness; the hypotheses are for random geometric graphs approximating a known manifold/domain, not automatically satisfied by our combinatorial PL dual graphs.
   - https://www.math.cmu.edu/users/slepcev/spectral_gamma10.pdf
11. M. Keller, D. Lenz, R. Wojciechowski, *Graphs and Discrete Dirichlet Spaces* (2021). A systematic framework for graph Laplacians, regular Dirichlet forms, intrinsic metrics, heat semigroups, recurrence and stochastic completeness.
   - https://www.math.uni-potsdam.de/professuren/graphentheorie/research/ansicht/graphs-and-discrete-dirichlet-spaces

These works provide mathematical methods and comparison points; they do not validate SBQT's physical interpretation.

## Required mathematical object at every refinement level

For each level n, define explicitly:

- a finite measured space (X_n, d_n, mu_n), including the metric construction and normalization of mu_n;
- a Hilbert space H_n = L^2(X_n, mu_n);
- a dense domain D(E_n) and a nonnegative, closed, symmetric quadratic form E_n;
- the self-adjoint operator A_n associated with E_n;
- maps/comparison operators between H_n and a candidate limit Hilbert space H, or a common ambient space;
- the physical length scale h_n and every normalization factor in A_n.

For a weighted graph, a candidate form has the structure
E_n(f) = c_n/2 sum_{(i,j) in E_n} w_{ij}^{(n)} |f_i - f_j|^2,
with the scalar c_n, edge weights, vertex measure, and inner product derived from the model's geometry. They must not be chosen after inspecting the desired spectral dimension. The associated operator is the operator represented by this form in L^2(X_n, mu_n), not an unlabelled adjacency matrix.

The current G3 coarse matrix U^T L U is a finite fibre-compressed operator. Its construction is useful, but it is not by itself a derivation of the continuum normalization c_n or a proof that the sequence of forms has a limit.

## Critical source-level audit: the current graph operator is not yet a geometric Laplacian

The current exporter starts from the boundary of the four-dimensional cross-polytope: 16 tetrahedral facets. It repeatedly applies barycentric subdivision, which creates 24 child tetrahedra per parent tetrahedron. Thus the generation-3 count 221,184 is exactly 16 x 24^3.

There is a useful canonical *reference embedding*: assign the eight cross-polytope vertices to the points ±e_i in R^4 and place each new barycentric vertex at the Euclidean barycentre of its parent face. Each initial facet is a regular tetrahedron of side sqrt(2) and 3-volume 1/3. The total piecewise-flat boundary 3-volume is therefore 16/3, and barycentric subdivision preserves it. At generation 3, every child tetrahedron has volume 1/(3 x 24^3), with total volume 16/3.

By contrast, the current chunked exporter records measure_total = n_child / 24^3 = 16. This is a combinatorial counting normalization (one unit per initial tetrahedron), not the Euclidean PL 3-volume. The constant factor of three may cancel in some mass-normalized operators, but the two measures must not be given the same geometric interpretation. The continuum program must persist both separately and state which one is used by each operator.

A second distinction is equally important: the current fine graph is the **unweighted dual adjacency graph** of the tetrahedra. Its unit-edge combinatorial Laplacian is a legitimate graph operator, but it is not automatically a finite-volume or finite-element approximation to the Laplace--Beltrami operator of the reference PL metric. A geometric finite-volume candidate would use face-area / dual-distance conductances and cell-volume masses; a primal FEM candidate would use the stiffness and mass matrices obtained by integrating piecewise-linear basis functions on the tetrahedra. These are separate operators and must be compared rather than silently identified.

Therefore the next geometry task is to export, at every generation, (i) the barycentric coordinates of all PL vertices, (ii) each tetrahedron's intrinsic 3-volume and diameter, (iii) shared-face areas and cell-to-face altitudes, (iv) the maximum mesh diameter h_n and mesh-quality bounds, and (v) both combinatorial and geometric measures. Then construct at least one geometry-weighted dual finite-volume operator and one primal FEM operator. Only after this can we test whether the unweighted dual Laplacian, the geometric finite-volume operator, and the FEM operator share the same rescaled low-energy limit.

This is a source-level mathematical distinction, not a claim that the reference embedding is already the physical spacetime metric of SBQT.

## Shape-regularity obstruction measured on G0--G3

The geometric audit artifact (workflow 44446384790) adds a dimensionless tetrahedral quality
q = 12 (3 V)^(2/3) / sum_{edges} l_e^2,
which equals 1 for a regular tetrahedron and tends to 0 for degenerate cells. It also records the maximum diameter-to-inradius ratio. The measured sequence is:

| Level | Minimum q | Median q | Maximum diameter / inradius |
|---|---:|---:|---:|
| G0 | 1.0000 | 1.0000 | 4.90 |
| G1 | 0.5968 | 0.5968 | 10.35 |
| G2 | 0.1954 | 0.4250 | 29.77 |
| G3 | 0.0531 | 0.2844 | 88.44 |

Thus mesh diameter decreases, but the worst-cell shape quality deteriorates markedly. Uniform shape regularity is not established. Standard finite-element error estimates that assume a uniformly shape-regular family cannot be invoked for this barycentric sequence without a separate argument. This is a real obstruction to a quick FEM convergence claim.

The next geometry gate must compare at least two alternatives: (i) retain barycentric subdivision and prove estimates that tolerate its measured anisotropy, and (ii) construct a quality-controlled refinement family (for example, a compatible red-refinement or mesh-optimization scheme that preserves the same PL complex and causal/fibre structure). Any optimized family must preserve topology, orientation, face pairing, and pushforward measure, and its results must be compared against the unoptimized reference family. We must not optimize cells in a way that silently changes the model's combinatorics.

## Continuum acceptance theorem to target

The primary target is Mosco convergence of E_n to a closed nonnegative form E on H, with explicit identification of the Hilbert-space convergence/comparison maps. The proof obligations are:

1. **Liminf:** whenever u_n converges weakly to u under the chosen varying-Hilbert-space notion, E(u) <= liminf E_n(u_n).
2. **Recovery:** for every u in D(E), construct u_n converging strongly to u with limsup E_n(u_n) <= E(u).
3. **Compactness / tightness:** establish the hypotheses needed for compact resolvent or, if the limit is noncompact, state the correct local or measure-theoretic spectral mode of convergence.
4. **Measure and geometry:** prove convergence/tightness of mu_n and control metric distortion, volume growth, and degenerating bottlenecks.
5. **Operator consequence:** derive strong resolvent convergence; claim eigenvalue/eigenfunction convergence only when the additional compactness and multiplicity hypotheses hold.
6. **Heat semigroup:** derive convergence of e^{-t A_n} on the appropriate domain of t. Heat-trace convergence requires additional uniform spectral-tail/trace-class control; it does not follow from strong resolvent convergence alone.

A numerical fit to a power law is not a substitute for these conditions.

## Refinement protocol: eliminate false convergence

A single nested sequence is insufficient. For each target scale, generate at least three independently constructed refinement families (different admissible subdivision orders, relabellings, and local refinement choices) with matched macroscopic geometry and measure. Where the model admits it, add controlled perturbations of edge weights and fibre assignments.

At each level record h_n, vertex/edge counts, measure total, degree distribution, fibre-size distribution, graph diameter, bottleneck/isoperimetric diagnostics, low spectrum with residuals, spectral counting function, heat trace, and resource use. Preserve raw operator and measure artifacts with hashes.

Compare observables only after the same dimensional scaling is derived and applied to all sequences. Report uncertainty and residuals; never select the time window or exponent after seeing which one gives the preferred dimension.

## Required convergence diagnostics

### A. Quadratic-form and resolvent layer

- Construct explicit restriction/prolongation or comparison maps.
- Test energy liminf and recovery inequalities on growing deterministic function families, including constants, localized bumps, low eigenmodes, and random smooth-proxy fields.
- Compare resolvents (A_n + z)^(-1) on matched probes for multiple positive z.
- Track the operator residual and the error under independent refinement families.

Finite probe tests are falsification tools, not a proof of Mosco convergence. The eventual proof must identify a dense core and uniform estimates.

### B. Spectrum and heat kernel

- Compute low eigenpairs with residuals and multiplicity-aware subspace angles.
- Track eigenvalue counting functions and heat traces over a scale window that expands with refinement.
- Separate ultraviolet, crossover, and infrared regimes.
- Demonstrate that the candidate scaling window moves with h_n while the rescaled curve stabilizes.
- Quantify finite-size corrections and compare at least two independent extrapolation ansätze.
- Test robustness to changing the heat-time grid and differentiation method for effective spectral dimension.

A plateau is accepted only if it persists across refinement families, expands or remains stable in rescaled coordinates, and its uncertainty shrinks under refinement. A local flat segment on one level is not enough.

### C. Geometry and measure

- Verify that total measure conservation is accompanied by convergence of measures on test functions, not merely equality of total mass.
- Measure volume growth and ball-volume ratios at matched rescaled radii.
- Track diameter, injectivity/bottleneck proxies, and isoperimetric/Cheeger quantities.
- Check that no vanishing-measure appendages or exceptional fibres dominate the heat trace.

### D. Feshbach / coarse-graining layer

For P_n the coarse projection and Q_n = I - P_n, define the energy-dependent Schur complement with an explicit spectral parameter z:
F_n(z) = P_n(A_n-z)P_n - P_n A_n Q_n [Q_n(A_n-z)Q_n]^{-1} Q_n A_n P_n,
on a domain where the inverse exists.

Required:
- certify or bound the Q-block spectral gap and resolvent norm over a common z-window;
- calculate the full coarse self-energy matrix when feasible, otherwise report matrix-free operator actions on a systematically expanding, reproducible probe basis with error bounds;
- compare independent probe bases, not only four seeded vectors;
- test symmetry, positivity/monotonicity where theoretically required, and finite-difference agreement for the self-energy derivative;
- establish convergence of F_n(z) after the same coarse-space identification maps;
- show that truncation and solver errors are below the refinement error.

Four positive quadratic-form probes are a useful initial smoke test, not a positivity theorem for the full matrix and not a continuum closure.

### E. Renormalization and universality

If the graph ensemble or action has tunable couplings, identify dimensionless couplings and a finite-size scaling prescription before fitting. Search for critical surfaces/fixed points through multiple volumes and independent bare-coupling trajectories. Record relevant/irrelevant directions and stability under regulator/parameterization changes. A stable fixed point and universality across discretization families are stronger evidence than apparent convergence of one deterministic graph sequence.

## Immediate SBQT work queue

1. **Freeze definitions:** write down the G_n geometry, h_n, mu_n, comparison maps, E_n, and A_n. Explicitly distinguish the full 221,184-vertex G3 graph Laplacian, the 512-fibre Galerkin compression, and the Feshbach effective operator.
2. **Audit normalization:** derive c_n and weights from the intended continuum energy/dimension. Add dimensional-analysis tests and prohibit post-hoc per-level normalization.
3. **Build independent refinement families:** produce at least three G0--G3 families with matching macroscopic invariants; compare spectra and heat traces without pooling unlike operators.
4. **Add form-convergence harness:** implement liminf/recovery probe suites and resolvent-action comparisons; persist all probe vectors/seeds and errors.
5. **Expand Feshbach:** replace the four-vector-only report with deterministic basis blocks and randomized cross-checks; compute the complete 512-dimensional coarse self-energy if memory/time permit. If not, produce rigorous operator-norm or quadratic-form bounds, not a false PASS.
6. **Extend levels:** G4 and beyond only after the normalization and exporter are stable; retain bounded-memory export and verify geometry/measure at each generation.
7. **Independent reproduction:** run the same analysis from a second implementation or independently reconstructed artifacts.
8. **Physical bridge remains separate:** even a proven spectral continuum limit does not prove Lorentzian causal propagation, Einstein dynamics, or physical QFI conductance. Those require independent operator and observable derivations.

## New finite-form implementation and regression result

The research branch now contains a finite graph Dirichlet-form module, with explicit positive vertex measure, positive edge conductances, scale factor, stiffness matrix, energy, mass norm, associated generator, and mass-weighted resolvent action. It also validates nonnegative constant-preserving prolongation maps and computes candidate recovery-sequence probe diagnostics. The report deliberately labels these as finite probes, not a Mosco proof.

The dedicated Q2 graph/export workflow run 44445748662 completed successfully on commit d453900eab2ba4d67d644d6e087b89c988918520. Its 23 pytest invocations passed 150 tests in total; the new Dirichlet-form suite contributed 5/5 passing tests. The run also completed the existing graph, QFI, G3 exporter, and Feshbach-control suites and uploaded its artifact.

This is a tested numerical foundation for the next stage. It does not yet provide actual SBQT geometric conductances, independent refinement families, a liminf theorem, or a recovery theorem on a dense core.

## New finite-level result: complete G3 Feshbach matrices (19 November 2026)

Dedicated workflow run 44442960180 completed successfully on the research branch and uploaded artifact 14493520542. The artifact contains both complete 512 x 512 matrices and the manifest. It is the full-column reconstruction, not the four-probe matrix-free report.

- 512/512 coarse basis columns solved for both the correction and self-energy slope.
- 230,632 total projected-CG iterations; solve time 1,374.05 s; peak RSS 148,156,416 bytes.
- Maximum first/second solve residuals: 1.4085e-8 and 1.4962e-8.
- Correction: symmetry error 1.6228e-9, minimum eigenvalue 2.3506e-15, maximum 2.0440, trace 348.2500.
- Slope: symmetry error 1.1520e-7, minimum eigenvalue -1.3963e-15, maximum 27.8169, trace 514.2515.
- No eigenvalues below the workflow's negative tolerance in either matrix.

An independent artifact-side check reconstructed the expected constant coarse vector q_i = sqrt(|F_i|), normalized it, and found alignment 1.0 with the numerical null eigenvector of each matrix. The residual norms were 3.17e-9 for the correction and 3.94e-8 for the slope. The remaining 511 eigenvalues are positive at the numerical resolution of this artifact.

This is a strong finite G3 Schur/Feshbach result. It does not supply the missing geometric conductances, a refinement sequence of these matrices, Mosco convergence, or a physical QFI interpretation. The matrices are generated from the unweighted dual graph. Artifact run: https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44442960180

## New finite-level result: complete G3 Feshbach matrices (19 November 2026)

The dedicated GitHub Actions run completed successfully:
- Run: https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44442960180
- Artifact: `q2-g3-full-feshbach-matrices` (artifact ID 14493520542).
- Actual graph: 221,184 vertices, 442,368 edges, 512 nonempty coarse fibres, one connected component.
- All 512 coarse basis columns were solved for both the correction and self-energy-slope matrices (not four probes).
- 230,632 total projected-CG iterations; 1,374.05 seconds solve time; max first/second solve residuals 1.41e-8 / 1.50e-8.
- Correction matrix: symmetry max-entry error 1.62e-9; smallest eigenvalue 2.35e-15; largest 2.0440; trace 348.250.
- Self-energy slope: symmetry max-entry error 1.15e-7; smallest eigenvalue -1.40e-15 (roundoff); largest 27.8169; trace 514.252.
- Both matrices are 512 x 512 and numerically positive semidefinite to the stated tolerances. Independent artifact inspection also verified the normalized constant coarse mode (proportional to square-root fibre sizes) is annihilated to relative residual about 1.6e-9 / 1.4e-9.

Independent matrix-action cross-check on the persisted G3 graph (seed 20261119, one normalized random coarse vector, independently rebuilt sparse Laplacian): projected-CG converged on both solves with residuals 1.25e-9 and 8.89e-10. Relative discrepancy between direct operator action and stored matrix action was 8.16e-9 for the correction and 1.40e-7 for the slope. This closes the earlier four-probe limitation for the **finite G3 geometric graph**. It does not establish convergence of these matrices across refinement levels, nor identify them with physical QFI conductance. The solve residuals are numerical stopping diagnostics, not rigorous operator-norm error bounds. The next test is matrix-level refinement stability across independently generated G0--G3 families, using common coarse-space identifications and fixed normalization.

## Current decision rule

Status labels must distinguish:

- FINITE_LEVEL_PASS: finite graph/operator invariants pass;
- NUMERICAL_SCALING_EVIDENCE: multi-level data show a stable trend but theorem assumptions remain unproved;
- MOSCO_CONVERGENCE_PROVED: liminf and recovery plus the required Hilbert-space/measure hypotheses are proved;
- SPECTRAL_CONVERGENCE_PROVED: operator/spectral consequences are established with compactness or the correct noncompact substitute;
- CONTINUUM_PHYSICAL_LIMIT: only after the mathematical limit and independent physical interpretation gates both pass.

At the present G0--G3 stage, the correct label is finite-level evidence only. No stable continuum limit is claimed.

## Immediate references for the next implementation pass

- Kuwae--Shioya (2003): varying Hilbert spaces, Mosco convergence, spectral structures.
- Babuška--Osborn (1989): Galerkin eigenvalue/eigenvector error estimates.
- Kigami (2001): compatible network energies and limiting resistance forms.
- Ambjørn et al. (2016): lattice-spacing determination and finite-size scaling in CDT.
- Dusson--Sigal--Stamm (2021, 2023): Feshbach--Schur estimates and discretization error.
- Elek (2008): weak graph convergence and integrated density of states.

The first implementation milestone is not another plot. It is a fully specified, scale-normalized sequence of energy forms with comparison maps and a test suite capable of detecting failure of liminf, recovery, and resolvent convergence.
