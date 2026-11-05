# BQG: G25 deep pass and post-sum charged closure (2026-11-06)

## Verified real controls

G21 genuine node1 KA auxiliary 000: 21/21 rank, 1202 distinct outgoing states, coherent norm 0.0677323644561008; remaining 99 of 120 source states with norm 0.0872045595018736. Corrected direct two-input SU(2) coherent linearity is PASS in run 42113437487: difference relative to direct norm 8.790476459045389e-16, all 36/36 output keys match after normalizing nested archive key to genuine flat output key.

## Actual G25 run

Run 42128476913 computes genuine held-out full KA sources 31,64,53,111 (node1, column0, T023, aux 000, Jmax=7/2, same frozen KB SHA). Do NOT claim 25/25 until all four complete and full complex outputs pass SHA and common source-basis metadata.

If all four succeed: 25/120 inputs, 64.214229435% of squared input weight; residual 95 input states with input norm 0.07684417184116164. These inputs DO NOT bound outgoing error without a norm on a fixed residual ray.

## Charged-basis leak: new exact criterion (not yet evaluated)

The existing complete_charge_basis_leakage=1.0 comes from maximum over single primitive unsummed branches. This is NOT a measured global post-summed gauge leak.

For a common orthogonal projector P on one fixed full Hilbert sector and the unprojected coherent operator result R=sum_alpha c_alpha phi_alpha:

  ||(I-P)R||^2 = ||R||^2 - ||PR||^2

Hence define the actual post-summed relative leak as sqrt(max(0,||R||^2-||PR||^2))/||R||. Different spin sectors must remain orthogonal, and intermediate output charge labels must share one inner-product convention.

One can avoid huge dense state tensors by factorizing inner products of genuine product spin-network branches:

  ||R||^2 = sum_{alpha,beta: same spins} conj(c_alpha)c_beta * product_{v=0..4} <X_alpha_v,X_beta_v>.

Naively this O(M^2) branch Gram is expensive; group spins, hash local tensor sectors, reuse SU(2) local contractions. Only then may the 1.0 primitive leakage be pronounced canceled or surviving on the actual coherent output.

## Prioritized gates

1. Verified G25 (rank, spectrum, overlap and coherent output norm).
2. Direct actual residual of 95 unsolved input sources A Q95 psiKB; check required ||delta|| <= ||Psi25||/11 for at most 10% error against full output.
3. Physical charge closure after ALL primitive contributions combined, not branchwise max.
4. Other seven auxiliary channels, 24 oriented terms, conjugate outgoing, source0 reference, ten independent off-shell HDA constraints, physical projector, TT pole.

Full Lorentzian HDA, interacting physical TT pole, and physical projector remain OPEN.
