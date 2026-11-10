# Findings 017 — exact Bernoulli graph control, 190 independent realizations

## 1. Exact ensemble definition

This control replaces the previous Poissonized shortcut sampler with the intended independent Bernoulli law. For a periodic cubic lattice of side (L), (N=L^3), retain the six deterministic nearest-neighbour edges. For every other unordered pair ({u,v}), let (r(u,v)) be periodic Euclidean distance and include the edge independently with
[
q_{uv}=\frac{p}{r(u,v)^2},\qquad p=\frac{1}{6L}.
]
No Poisson approximation is used.

Efficient exact sampling groups pairs by displacement class. If a class contains (M) eligible pairs with common probability (q), draw (K\sim\mathrm{Binomial}(M,q)), then choose a uniform (K)-subset. This is exactly equivalent in distribution to (M) independent Bernoulli trials. Opposite displacement classes are represented once; self-inverse classes use (M=N/2).

Implementation: `exact_bernoulli_gap.py`.

## 2. Numerical results

Sparse combinatorial Laplacian, λ₂ computed by symmetric Lanczos/ARPACK (`eigsh`, smallest-magnitude eigenvalues). All sampled graphs in the reported batches were connected.

| L | N | realizations | mean λ₂ | sample SD | CV | mean degree |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 512 | 20 | 0.951126727 | 0.032443856 | 3.4111% | 6.9678 |
| 10 | 1,000 | 20 | 0.779016677 | 0.027049980 | 3.4723% | 7.0056 |
| 12 | 1,728 | 20 | 0.694949108 | 0.018431898 | 2.6523% | 7.0666 |
| 14 | 2,744 | 20 | 0.650590924 | 0.015728806 | 2.4176% | 7.0946 |
| 16 | 4,096 | 20 | 0.617000356 | 0.015128739 | 2.4519% | 7.1219 |
| 20 | 8,000 | 20 | 0.574706789 | 0.010051839 | 1.7490% | 7.1572 |
| 24 | 13,824 | 20 | 0.553933989 | 0.007278427 | 1.3139% | 7.1726 |
| 28 | 21,952 | 20 | 0.541104641 | 0.006253781 | 1.1557% | 7.1874 |
| 32 | 32,768 | 20 | 0.534760183 | 0.005626131 | 1.0521% | 7.2018 |
| 40 | 64,000 | 10 | 0.526079408 | 0.003480680 | 0.6616% | 7.2144 |

Total: 190 independent graph realizations. For L=40, ten more attempts at L=48 were started but the eigensolver batch was interrupted by the execution time limit; no L=48 statistic is claimed.

Approximate Student-(t) 95% confidence intervals for the batch means (descriptive, conditional on this sampling protocol):
- L=8: [0.935943, 0.966311]
- L=16: [0.609920, 0.624081]
- L=24: [0.550528, 0.557340]
- L=32: [0.532127, 0.537393]
- L=40: [0.523590, 0.528569]

## 3. Interpretation

1. The exact-Bernoulli sample means continue to decrease across the tested finite sizes, while the coefficient of variation of λ₂ decreases from about 3.4% to 0.66%. This is evidence of increasing concentration over this finite range, not a proof of concentration as (N\to\infty).
2. The observed gap is approaching a nonzero-looking plateau near 0.5 over the largest sizes, but the current range does not establish a positive limiting gap. A slow decay, crossover, or positive limit remains to be discriminated by larger (L), more seeds, and concentration bounds.
3. The exact-Bernoulli means are close to the previous Poissonized-control means at overlapping sizes (e.g. L=16, 24, 32), but the samplers are not identical and should remain separate in the archive.
4. This random-graph λ₂ is not the same object as λ₂ of the expected Laplacian. The expected-operator coordinate-mode limit near 0.72835 cannot be substituted for the expectation or typical value of the random graph gap; the exact-Bernoulli sample means are materially lower.
5. These results test only the specified cubic-torus plus inverse-square Bernoulli ensemble. They do not derive the cubic background, the inverse-square kernel, or the coefficient (p=1/(6L)) from BQG.

## 4. Reproducibility and limitations

- Exact independent Bernoulli sampling is used per unordered pair.
- Fixed deterministic seed schedules; every realization is independently seeded.
- λ₂ numerical accuracy depends on ARPACK convergence; the script raises on solver failure rather than silently replacing the result.
- The reported batches have no failed connectivity cases.
- No finite-size scaling law is selected as final. In particular, the prior ((\ln N)^{-2}) claim remains withdrawn.
- The full CSV can be regenerated with the script; no claim is made that an interrupted L=48 batch completed.

## 5. Verdict

- Exact Bernoulli sampler: **IMPLEMENTED**.
- Finite-size spectral-gap sample: **190 realizations, L=8…40**.
- Decreasing sample CV: **OBSERVED**, not an asymptotic concentration theorem.
- Positive limiting spectral gap: **PLAUSIBLE BUT UNPROVED**.
- Expected-operator spectrum = random-graph expected gap: **NOT VALID**.
