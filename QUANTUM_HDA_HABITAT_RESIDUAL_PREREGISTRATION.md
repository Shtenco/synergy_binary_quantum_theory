# Quantum HDA habitat residual — production preregistration

## Status

This is a **fail-closed microscopic gate**.  It does not certify the existing target-side construction by itself.  It certifies only an explicit quantum packet on a declared graph-changing habitat.

The production projector is not emitted until this gate and the constraint-master gate both pass on compatible hashes.

## Frozen operator decomposition

Write

\[
H_v = H_{E,v}+\lambda H_{L,v},\qquad \lambda=1+\beta^2.
\]

For every ordered constraint pair \((k,l)\) and every declared boundary/generated column \(i\), the quantum commutator packet is decomposed before any numerical cancellation:

\[
HH_{EE}^{kl,i}=[H_{E,k},H_{E,l}]|\psi_i\rangle,
\]

\[
HH_{EL}^{kl,i}=[H_{E,k},H_{L,l}]|\psi_i\rangle,
\]

\[
HH_{LE}^{kl,i}=[H_{L,k},H_{E,l}]|\psi_i\rangle,
\]

\[
HH_{LL}^{kl,i}=[H_{L,k},H_{L,l}]|\psi_i\rangle.
\]

The independently constructed target may also be represented as a polynomial in the same frozen convention,

\[
D^{target}=D_0+\lambda D_1+\lambda^2D_2.
\]

The residual components are therefore

\[
r_0=HH_{EE}-D_0,
\]

\[
r_1=HH_{EL}+HH_{LE}-D_1,
\]

\[
r_2=HH_{LL}-D_2,
\]

and

\[
\Delta(\lambda)=r_0+\lambda r_1+\lambda^2 r_2.
\]

No post-hoc choice of \(\beta\) or \(\lambda\) is allowed.

## Frozen norm diagnostic

For the aggregate pair/column residual, the certificate records the full quartic polynomial

\[
\|\Delta(\lambda)\|^2
=c_0+c_1\lambda+c_2\lambda^2+c_3\lambda^3+c_4\lambda^4,
\]

with

\[
c_0=\langle r_0,r_0\rangle,
\qquad
c_1=2\,\mathrm{Re}\langle r_0,r_1\rangle,
\]

\[
c_2=\langle r_1,r_1\rangle+2\,\mathrm{Re}\langle r_0,r_2\rangle,
\]

\[
c_3=2\,\mathrm{Re}\langle r_1,r_2\rangle,
\qquad
c_4=\langle r_2,r_2\rangle.
\]

The gate never minimizes this polynomial in \(\lambda\).  It may evaluate it only at a **preregistered microscopic/convention value** supplied by the packet.

If \(r_0=r_1=r_2=0\) within the frozen residual tolerance, the result is stronger: quantum HDA closure is beta-independent on the declared finite domain.

## Mandatory common metadata

A physics packet must declare nonempty:

- `habitat_hash`
- `domain_hash`
- `constraint_packet_hash`
- `convention_hash`
- `jmax`
- complete expected constraint-pair list
- complete expected column list
- cutoff leakage
- habitat leakage
- recoupling error
- frozen tolerances

Every expected pair × column record must occur exactly once.  Missing, duplicate, unexpected, or hash-inconsistent records fail closed.

## Refinement test

When `require_refinement_certificate=true`, a single finite residual is not enough.  The preregistered fit is

\[
\|\Delta(r)\|=A r^p+B,
\]

with a deterministic grid \(p\in[0.10,6.00]\) in steps of 0.01 and linear least-squares solution for \(A,B\) at each grid point.  The best SSE point is reported.

The refinement gate requires

\[
|B|\le B_{tol},\qquad p\ge p_{min}.
\]

A nonzero stable intercept is an anomaly/blocker, not a quantity to absorb into a fitted coupling.

## Certification rule

`quantum_habitat_residual_certified=true` iff all of the following hold:

1. common hashes are present and record-consistent;
2. the declared pair × column domain is complete;
3. cutoff/habitat leakage and recoupling error satisfy frozen tolerances;
4. either all three component residuals close independently, or the total residual closes at a preregistered \(\lambda\);
5. the refinement test passes whenever refinement is required;
6. no fitted \(\lambda\) is used for certification.

The following are explicitly **insufficient**:

- an existing `Dtarget`/classical target certificate by itself;
- small Euclidean-only \(HH_{EE}-D_0\) when Lorentzian terms are required;
- logical-sector projection only when outgoing graph-changing support exists;
- a best-fit value of \(\beta\) or \(\lambda\);
- a residual computed with mismatched habitat/domain/convention hashes;
- one finite-regulator point when the production preregistration requires refinement.

## Production dependency

The legal front-end remains

\[
\{C_A\}
\rightarrow M_{\rm BQG}
\rightarrow P_{\rm BQG}
\rightarrow Z[J]
\rightarrow W[J]
\rightarrow \Gamma_{\rm BQG}.
\]

The HDA gate is a prerequisite for `P_BQG`; it is not itself evidence for dark matter, dark energy, lensing, or a scalar pole.

Only after the physical projector/history exists may scalar and FLRW observables be derived from the same physical \(\Gamma_{\rm BQG}\).
