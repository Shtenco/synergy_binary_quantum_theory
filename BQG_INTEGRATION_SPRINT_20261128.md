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
