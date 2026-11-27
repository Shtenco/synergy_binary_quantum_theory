# W5 SU(2) structural basis contract

## Purpose

`scripts/bqg_w5_su2_structural_basis.py` enumerates admissible representation
labels for a fixed ten-link K5 spin assignment. Its support is defined only by
SU(2) tensor-product rules; no amplitude cutoff, numerical coefficient, CG
overlap, or eigensolver output is consulted.

## Conventions

All spins are doubled integers. For two representations,

`s1 tensor s2 = |s1-s2|, |s1-s2|+2, ..., s1+s2`.

At a four-valent node, pair channels are
`K12 in s1 tensor s2` and `K34 in s3 tensor s4`.

* Gauss-singlet labels obey `K12=K34=K`.
* An explicit charged source label obeys
  `|K12-K34| <= J2 <= K12+K34` and
  `K12+K34+J2 = 0 mod 2`, with all `M2=-J2,-J2+2,...,J2`.

The covariant K5 key is `(spins, labels)`: four nodes carry `("G",K)`;
the selected source carries `("C",J2,M2,K12,K34)`. This matches the label
conventions used by the existing covariant K-leg implementation.

## W5 seed

The five frozen seeds have ten doubled link spins equal to 4 and node labels
`K=(0,2,4,6,8)`. Their source-covariant embedding uses `J2=M2=0` and
`K12=K34=K` at the source node.

## Scope and non-claims

This is an exact representation-theoretic *ambient basis* for each supplied
spin sector. It is not the operator-reachable support of the full 24-term
epsilon KKV map. The operator may occupy a strict subspace because of
intertwiners, orientation, holonomy-index contractions, volume blocks and
coherent cancellations. Conversely, this generator does not decide whether
a coefficient is nonzero.

The full W5 reachable-basis problem still requires a word-aware support
propagator through the actual forward and independently reversed-adjoint
operator words, with every intermediate link-spin sector represented. No
physical projector, HDA closure, or nonzero physical kernel is implied.
