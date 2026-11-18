# Graph topology dependence of finite vertex-fusion Gram maps (2026-11-19)

## Question

The existing independent-local-Haar pushforward result uses connected path blocks. A path is a tree. Does the same isometric contraction structure survive when the fused block branches or contains a loop?

We construct the finite linear map
\[
C:\bigotimes_{v\in V}\mathcal H^{\mathrm{inv}}_{4,j}
\longrightarrow \mathcal H^{\mathrm{inv}}_{4|V|-2|E|,j}
\]
by contracting every internal edge of one graph block with the normalized invariant spin-\(j\) singlet. Its Gram operator \(G=C^\dagger C\) controls both the kernel of fusion and the norm of the resulting block tensor.

## Results

| Graph | Spin | Input dimension | Boundary legs | Boundary singlet dimension | Rank of \(C\) | Rank defect | \(\sigma_{\max}\) | Smallest positive \(\sigma\) | \(\kappa(C|_{\mathrm{supp}})\) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| star_4 (tree) | \(1/2\) | 16 | 10 | 42 | 16 | 0 | 0.125000 | 0.125000 | 1.000000 |
| cycle_4 (one loop) | \(1/2\) | 16 | 8 | 14 | 12 | 4 | 0.125000 | 0.0589256 | 2.121320 |
| cycle_4 (one loop) | \(1\) | 81 | 8 | 91 | 61 | 20 | 0.0370370 | 0.00979908 | 3.779645 |

For the four-node star, \(C^\dagger C=2^{-6}I_{16}\) at spin \(1/2\), up to a maximum absolute residual of \(1.39\times10^{-17}\). This is the tree result expected from leaf elimination: each normalized singlet contraction contributes one factor \(d^{-2}\) to the double-layer Gram map. The maximum diagonal fused-block Gauss residual over all product-basis columns is \(8.31\times10^{-16}\).

For the four-cycle, the Gram matrix is not proportional to identity and has an exact numerical kernel in the tested finite sectors:
- spin \(1/2\): rank \(12/16\), four annihilated input directions;
- spin \(1\): rank \(61/81\), twenty annihilated input directions.

The fused tensors nevertheless satisfy the diagonal boundary Gauss law to residuals below \(4.5\times10^{-15}\). Thus the rank loss is not a failure of the invariant boundary constraint; it is a property of the loop contraction map itself.

## Consequence for normalized ensembles

For a tree block, the scalar Gram identity makes \(\|C\psi\|\) independent of the normalized microscopic input \(\psi\). The normalized map is therefore a linear isometry up to a fixed scale, which justifies the simple covariance pushforward used in the chain Haar calculation.

For a loop block, \(C^\dagger C\) is non-scalar and rank-deficient. The normalized state map
\[
|\psi\rangle\mapsto \frac{C|\psi\rangle}{\|C\psi\|}
\]
is input-dependent and nonlinear; an ensemble pushforward must account for the factor \(1/\|C\psi\|^2\) and the kernel support. It is incorrect to reuse the tree formula \(C(I/D)C^\dagger\) as the covariance of normalized outputs without evaluating that reweighting.

This gives a concrete reason that graph topology matters for a future RG construction: tree blocks and loop blocks do not induce the same measure transformation, even before introducing a physical BQG history measure.

## Verification and scope

The regression suite reconstructs the complete finite contraction matrix from an orthonormal local intertwiner basis, checks the Gram spectrum, rank, boundary singlet dimension, and diagonal boundary Gauss residual. The CI workflow exports separate JSON controls for the star and cycle fixtures at the listed spins.

This is a finite kinematic result for the named graph fixtures and spins. It does not derive the physical BQG measure, sum over spins or internal labels, define a physical boundary trace, or establish a continuum RG flow or fixed point.
