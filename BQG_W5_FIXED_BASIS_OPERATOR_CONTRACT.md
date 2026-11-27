# Fixed-basis linear operator contract for W5

This is the next infrastructure gate after the threshold audit.

## Key design decision

Do not discover the operator domain by thresholded output support. First define a finite ordered basis using exact representation-theoretic selection rules. Then compute every column of the operator on every basis vector. The matrix is a linear map by construction; coefficients outside the declared basis are explicit leakage, not silently dropped.

The utility in \`scripts/bqg_w5_fixed_basis_operator.py\` implements this contract:
- no magnitude-tolerance parameter exists;
- every nonzero complex coefficient, however small, is retained;
- exact zero is omitted from sparse storage;
- output outside the fixed basis is returned as leakage;
- the conjugate-transpose matrix is constructed independently of any physical interpretation;
- homogeneity, additivity, and adjoint residuals can be measured on the assembled matrix.

Unit tests use a small complex matrix and include a coefficient of order \(10^{-200}\) to ensure it is not pruned.

## How to use it for BQG

1. Freeze the W5 input set \(K=(0,2,4,6,8)\), source-node set, link cutoff and a deterministic basis ordering.
2. Generate a conservative candidate basis using exact SU(2) fusion/triangle/parity rules, not numerical amplitude size.
3. Run the forward KKV column and the independently reversed-order adjoint column for every basis vector.
4. Assemble both matrices with this fixed basis utility.
5. Any output key outside the candidate basis is a structural closure failure. Expand the candidate basis by exact selection rules and rebuild; do not append only outputs above a tolerance.
6. Check matrix adjoint residual, cutoff stability and numerical error bounds.
7. Only after these gates, form the Hermitian constraint and inspect its non-vacuum kernel.

## Important limitation

This infrastructure does not yet generate the full W5 candidate basis and does not evaluate the BQG KKV columns. It therefore makes no claim about W5 rank, non-vacuum kernel, physical projector, HDA, TT modes or Einstein infrared dynamics. It removes one source of logical circularity from the future implementation.
