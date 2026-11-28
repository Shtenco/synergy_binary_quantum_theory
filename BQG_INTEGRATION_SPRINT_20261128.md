# BQG integration sprint — W5 → physical projector → HDA

Date: 2026-11-28  
Branch: research/bqg-physical-projector-hda-gates-20261128  
Base: research/bqg-fs-symplectic-reduction-20261128

## Objective

Close the logical interfaces between the W5 Lorentzian operator, finite physical projector, and microscopic hypersurface-deformation algebra. No new physical hypotheses are introduced.

## Work package A — W5 operator (first blocker)

The source audit establishes that current thresholded closure cannot be treated as a linear operator. The fixed-basis assembler removes magnitude pruning but explicitly states that it does not yet evaluate BQG KKV columns.

Required sequence:

1. Freeze W5 input seeds, source-node convention, graph orientation, epsilon ordering, link-spin cutoff, and orthonormal basis hash.
2. Generate a conservative basis from exact SU(2) triangle/parity rules.
3. Evaluate all 24 epsilon terms for every required source node and basis column; keep exact algebraic zeros separate from unresolved numerical intervals.
4. Accumulate all paths and epsilon terms before any display threshold.
5. Report forward leakage and independently ordered adjoint leakage.
6. Check homogeneity/additivity on superpositions, not merely on matrix assembly after the fact.
7. Increase cutoff and precision independently; retain a machine-readable error budget.

Stop condition: a state-limit guard, partial epsilon sum, unresolved coefficient interval, or omitted output key means W5_OPERATOR=OPEN, not PASS.

## Work package B — physical projector

Use one frozen set of constraint matrices from work package A and a declared common output fiber (or explicit direct-sum fibers with a positive block pairing). Then:

1. establish exact/interval Hermiticity of the constraint-label metric;
2. construct M=C†(G⊗I_K)C;
3. certify positivity and ker(M)=intersection ker(C_A);
4. compute the full non-vacuum kernel dimension, not only vacuum inclusion;
5. construct the orthogonal projector and certify P²=P=P†;
6. test stability under cutoff/refinement and precision changes;
7. report induced projective physical dimension and Fubini–Study rank.

The helper scripts/bqg_finite_master_kernel_gate.py is an exact finite reference implementation for this interface. It accepts supplied matrices; it does not yet load the full W5 operator.

## Work package C — HDA

Use exactly the same frozen constraint definitions as work package B. Evaluate DD, DH and HH independently, retaining the inverse-metric structure function in HH. Record off-shell, dual-habitat and on-shell results separately.

The acceptance criteria are specified in BQG_HDA_INTEGRATION_GATE.md. A successful master-kernel calculation does not imply HDA closure.

## Integrated acceptance

The combined physical consistency gate is PASS only if:

- W5 is a fixed-basis linear operator with complete forward/adjoint coverage;
- the physical kernel is non-vacuum and regulator-stable;
- the physical projector is exact or has a certified error bound;
- the same operators satisfy the declared HDA representation;
- the physical subspace is invariant under required observables/evolution;
- all conclusions survive the declared refinement limit.

Until then the integrated result is OPEN. A scientifically valid no-go is recorded as NO_GO and is not counted as a successful theory-closure gate.

## Known facts at sprint start

- The W5 threshold audit documents nonlinearity of hard amplitude pruning.
- The fixed-basis assembler retains arbitrarily small nonzero coefficients but does not evaluate physical W5 columns.
- The q=2 Fubini–Study calculation is a positive control, not the BQG physical projector.
- Exact vacuum inclusion has been shown for a recovered generic master-row construction, but its identity with the historical branch-sum operator and non-vacuum rank remain open.
- The W5 raw 24-term run was previously interrupted after one epsilon term; the full coherent column is unknown.
- The current ADM HDA calculation is classical and ansatz-conditional; microscopic BQG DD/DH/HH operator brackets remain open.
- A first staged W5 KKV attempt reached the charged H_E path but aborted because the highest-weight cache wrapper asserted J,M preservation and encountered an output label J2=2 from an input J2=0. This is not yet evidence of a physical anomaly: the assertion belonged to the optimization wrapper, not an independently verified operator identity.
- The highest-weight magnetic-sector shortcut has been disabled in the sine-ordered W5 path. The exact full magnetic labels now remain in H_E cache keys; the rerun must determine whether the underlying unsimplified operator satisfies the expected charge-sector selection rules.
- Current W5 run 8 evaluates only seed K=0 at source node 0 with Jmax2=10 and max_states=1. It is a diagnostic column, not the full W5 operator. The finite physical-projector CI passed again on the current tested head; microscopic HDA remains OPEN.

## Numerical audit of staged W5 run 8

The workflow label "exact-amplitude" is not a claim of exact arithmetic. The implementation in `k5_peter_weyl_safe_hda_column.py` evaluates Clebsch–Gordan coefficients to Python floats and stores branch tensors/amplitudes in NumPy complex floating-point arrays. The closure runner additionally uses the declared numerical support threshold `amp_tol=1e-12`. Therefore its output, if produced, is a finite-cutoff numerical witness, not an exact symbolic amplitude certificate.

The staged parameters `--K-values 0 --nodes 0 --Jmax2 10 --max-states 1 --max-depth 0` restrict the test to one W5 seed, one source node, and a spin wall J≤5. The state/depth limits do not establish closure of the reachable basis. Any output must be labelled partial unless it explicitly reports a completed adjoint-closed basis; a zero observed leakage at a finite spin wall is not by itself a proof of zero leakage beyond that wall.

The current field `final_nonscalar_J_weight_fraction` is not a fraction: it is hard-coded to 0 when no rejected coefficient was observed and 1 otherwise. It must not be interpreted quantitatively. The underlying `scalar_projection_rejected_norm` is a coefficient-norm diagnostic, but the implementation does not yet report a denominator from the total pre-projection covariant coefficient norm. A subsequent diagnostic patch must add that denominator and separate J, M, and charged-label rejection categories before any quantitative charge-leakage conclusion is drawn.

The finite projector helper now accepts an optional exact input Gram matrix Q. With Q supplied, it constructs the Q-self-adjoint master operator Q⁻¹C†(G⊗I)C and the Q-orthogonal projector N(N†QN)⁻¹N†Q; omitting Q retains the Euclidean/orthonormal-coordinate default. The W5 basis Gram matrix still has to be derived and supplied from the actual W5 inner product before this metric-aware path can certify the physical projector.

These are audit limitations, not findings of a physical inconsistency. W5, the full physical projector, and microscopic HDA remain OPEN.



## Exact Gram audit of the five W5 seed intertwiners

The five seed columns use doubled link spin 4 (j=2) on all ten edges and differ only by the source-node recoupling channel K2 in {0,2,4,6,8}. The normalized four-valent intertwiner is constructed by coupling legs (0,1) and (2,3) to the same K and then coupling K⊗K to the singlet, with normalization 1/sqrt(K2+1). Exact SymPy Clebsch–Gordan arithmetic gives the 5×5 intertwiner Gram matrix Q_seed=I_5. Since the remaining four node intertwiners are identical normalized K=0 states and the edge Peter–Weyl basis is assumed Haar-normalized, the full five-seed spin-network Gram matrix is also I_5.

A reproducible gate was added:
- `scripts/bqg_w5_seed_gram_gate.py`
- `tests/test_bqg_w5_seed_gram_gate.py`

The finite-projector CI workflow now runs both the finite master-kernel tests and the exact W5 seed Gram gate. The newly wired GitHub Actions workflow run 46178728858 completed successfully on head 657985b1f41fe1fc3f67c88c637a3c6b20a2efd8; both the combined test suite and standalone Gram certificate passed. This establishes the metric only on the declared five-dimensional seed span. It does not establish the Gram matrix of a larger dynamically generated closure basis, where duplicate/noncanonical states and basis normalization must be checked explicitly.

The metric-aware finite projector can therefore be supplied Q=I_5 for the five seed columns once the corresponding full constraint maps on that same five-dimensional input span and common output fiber are assembled. It still cannot certify the physical kernel until those W5 constraint matrices exist.


## Full K5 seed-network Gram factorization — 2026-11-28

The seed Gram gate has been strengthened from a local four-valent intertwiner calculation to an explicit factorized full-network certificate. The five seeds share all ten link labels (j=2); node 0 carries the five normalized K channels; nodes 1–4 carry the same normalized K=0 intertwiner in every seed. Under normalized Haar measure and normalized Peter–Weyl matrix elements, the link integrations identify endpoint magnetic indices by unitary contractions, so the spin-network inner product factorizes into the node intertwiner Gram factors. Node 0 contributes I_5 and each fixed node contributes [1], hence Q_full = I_5 and rank 5.

The implementation now exposes `exact_full_seed_gram()` and reports the factorization assumptions in the JSON certificate. A dedicated test checks the full-network 5×5 Gram result and rank. This is still a seed-span result: it does not calculate the Gram matrix of the dynamically generated adjoint-closure states, where basis collisions, distinct link labels, and sector-dependent normalization must be treated explicitly.

The expanded test suite was committed as `a00b275c4b4a0b647413a05ac5eecedb45a07bef`. GitHub Actions runs `46179668056` (finite-projector + W5 Gram) and `46179667988` (W5 first-column workflow) were queued at the time of this update; do not report the expanded suite as passing until the former completes successfully.
