# q=2 graph-link → Peter-Weyl representation bridge

## 1. Four-state-only obstruction

The frozen q=2 route cell has four basis states. The existing exact two-qubit SO(5) spinor quantum-link representation also has dimension four, and the q=2 Hamming adjacency embeds exactly in its operator algebra:

\[
A_{Q_2}=X\otimes I+I\otimes X
=2(M_{12}+M_{34}).
\]

However under the endpoint gauge algebra

\[
SU(2)_L\times SU(2)_R
\]

the four-state spinor decomposes as

\[
\mathbf4=(\mathbf2,\mathbf1)\oplus(\mathbf1,\mathbf2),
\]

not as the Peter-Weyl spin-1/2 link bi-doublet \((\mathbf2,\mathbf2)\). The left and right spin-1/2 Casimir projectors are complementary and have zero overlap.

Therefore a strict four-state-only identification is ruled out.

## 2. The missing state is already present: graph absence

The graph-changing cylindrical Hilbert space already contains a natural state in which a link is absent / reduced to j=0. No new fundamental local label is required.

Take

```text
4 active q=2 labels + 1 no-link state.
```

The existing five-state SO(5) vector quantum-link representation decomposes exactly as

\[
\mathbf5=(\mathbf2,\mathbf2)\oplus(\mathbf1,\mathbf1).
\]

The rank-4 active projector satisfies

\[
C_L=C_R=\frac34
\]

on the active sector, so it is precisely \(j_L=j_R=1/2\). The fifth state is an endpoint gauge singlet and is identified with the no-link state.

Every fundamental transporter component has only off-diagonal active/singlet blocks:

\[
P_g U_a P_g=0,
\qquad
P_0U_aP_0=0,
\]

and toggles active link ↔ no-link with unit norm.

## 3. Exact graph-changing factorization of the frozen q=2 adjacency

![SU(2) и Peter–Weyl representation bridge](docs/figures/05-su2-peter-weyl.svg)


For the four transporter components \(U_a=M_{a4}\),

\[
P_gU_aP_0U_bP_g=|a\rangle\langle b|.
\]

Thus every matrix unit on the four active q=2 labels is a two-step graph-changing excursion through the no-link state.

Using the **already frozen** q=2 Hamming edge set,

\[
A_{Q_2}
=
\sum_{a\sim_{Q_2}b}
P_gU_aP_0U_bP_g.
\]

So the microscopic q=2 adjacency has an exact implementation as

```text
active q=2 label
 -> no-link / j=0
 -> Hamming-neighbor active q=2 label.
```

Machine gate: `scripts/q2_graphlink_peter_weyl_gate.py`.

## 4. Higher-j representation growth by symmetric strand blocking

One active micro-link carries \((1/2,1/2)\). If a coarse link contains n indistinguishable active strands and the endpoint blocking is fully symmetric,

\[
\mathrm{Sym}^n(\mathbb C^2)_L
\otimes
\mathrm{Sym}^n(\mathbb C^2)_R
=
(j=n/2,j=n/2),
\]

with dimension

\[
(n+1)^2=(2j+1)^2.
\]

Allowing graph occupancy \(n=0,1,\ldots,N\) therefore supplies exactly one diagonal sector for every

\[
j=0,\frac12,1,\ldots,\frac N2,
\]

and

\[
\sum_{n=0}^N(n+1)^2
=
\sum_{j\le N/2}(2j+1)^2,
\]

which is the dimension of the Peter-Weyl tower truncated at \(J_{max}=N/2\).

`scripts/q2_symmetric_block_peter_weyl_growth_gate.py` explicitly constructs the symmetric SU(2) generators, verifies their commutators and Casimirs, and checks the counting identity.

This representation-growth theorem is **conditional on symmetric blocking**. The microscopic Hamiltonian has not yet been proved to dynamically select that symmetric coarse sector or its occupancy weights.

## 5. Claim boundary

The combined result substantially narrows `MICRO_TO_QGEOM`:

```text
q=2 active labels
 + graph-changing no-link state
 -> exact (2,2)+(1,1) quantum-link representation
 -> exact two-step graph-change realization of q=2 adjacency
 -> conditional symmetric-block growth of the full Peter-Weyl j tower.
```

Still open:

1. derive/uniquely select the SO(5) transporter Hamiltonian coefficients from the frozen microscopic action rather than choose a compatible completion;
2. show dynamical preference for the active geometric sector and symmetric blocking under coarse graining;
3. derive the physical occupancy weights / effective action over j;
4. connect that dynamics to the same semiclassical B-field and HDA regime without post-hoc projection.
