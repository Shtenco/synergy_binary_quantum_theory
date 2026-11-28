# BQG HDA integration gate — fail-closed specification

Date: 2026-11-28  
Status: **OPEN — no microscopic BQG HDA certificate is present in this gate.**

## Purpose

This document separates five logically distinct results:

1. classical ADM closure;
2. finite-regulator operator commutators;
3. continuum/dual or habitat closure;
4. recovery of the metric-dependent structure function;
5. preservation of the physical constraint space.

A classical ADM derivation, a finite toy matrix, or a master-constraint kernel does not by itself close the microscopic BQG HDA.

## Algebra to be recovered

With one frozen sign and smearing convention, the classical target is

\[
\{D[\vec N],D[\vec M]\}=D[[\vec N,\vec M]],
\]

\[
\{D[\vec N],H[M]\}=H[\mathcal L_{\vec N}M],
\]

\[
\{H[N],H[M]\}
=D[q^{ab}(N\partial_bM-M\partial_bN)].
\]

The last bracket is an algebroid relation: its structure function depends on the inverse spatial metric. It must not be silently replaced by a constant Lie-algebra structure coefficient.

## Required source objects

Before evaluating a bracket, freeze and hash:

- the actual BQG graph/representation carrier and inner product;
- the domain of each regulated constraint;
- lapse and shift discretization/interpolation;
- inverse-metric or inverse-volume prescription;
- operator ordering and adjoint convention;
- boundary terms and graph-changing output spaces;
- the map identifying output fibers for different constraints;
- regulator and refinement map.

If two constraints have different output fibers, their cross-pairing must be specified explicitly. A bare numeric matrix of coefficients is insufficient.

## Required bracket certificates

For each ordered pair of smearings and each tested basis state, store:

- independently evaluated left-hand commutator;
- independently evaluated right-hand constraint combination;
- residual before projection to the constraint surface;
- residual after any declared dual/habitat map;
- leakage outside the declared domain/output sector;
- precision and propagated numerical error bounds;
- regulator, graph, source-state and smearing hashes.

The DD, DH and HH brackets must be tested separately. HH is the decisive metric-dependent test.

## Off-shell and on-shell distinction

An identity that holds only after imposing constraints is not an off-shell representation of HDA. Report separately:

- OFF_SHELL_CLOSURE;
- DUAL_HABITAT_CLOSURE;
- ON_SHELL_ANOMALY_FREENESS.

Do not promote one status to another.

A master constraint can select a common kernel even when the original constraints do not form a strongly represented operator algebra. Therefore master-kernel existence and HDA closure are independent gates.

## Fail-closed acceptance criteria

The microscopic HDA gate remains OPEN unless all of the following are met:

1. all three brackets are evaluated from the frozen BQG operators;
2. HH reproduces the inverse-metric structure function with the declared ordering;
3. residuals vanish exactly or have certified error bounds converging under refinement;
4. no unexplained central/anomalous term remains;
5. domain and graph-changing leakage are controlled;
6. closure is shown on the declared off-shell domain, or the weaker dual/on-shell result is explicitly labelled;
7. the same operators preserve the declared physical space;
8. the classical/semiclassical limit recovers ADM HDA, not merely an Abelianized or gauge-fixed algebra.

## Current evidence boundary

The repository's ADM coefficient derivation is a classical consistency check inside an ADM ansatz. It does not establish microscopic BQG constraint commutators. The finite q=2 Fubini–Study test is a kinematic positive control. W5 structural-support tests do not evaluate amplitudes and do not establish HDA. The finite master-kernel helper certifies only explicitly supplied finite matrices.

The next physical HDA implementation must import the same frozen constraint definitions used by the physical-projector gate. If the two gates use different operators, the integrated BQG claim fails.
