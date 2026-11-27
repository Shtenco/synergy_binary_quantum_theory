# QAG explicit source-measurement dilation

**Date:** 2026-11-27  
**Scope:** finite-dimensional measurement-theory extension. This is not a gravity model.

## Purpose

The recoil model previously applied source postselection as an abstract rank-one projection. This module makes the measurement interaction explicit by appending a meter qubit, preparing it in |0>, correlating its computational label with the source measurement basis, and only then reading the meter.

For an orthonormal source basis {|b0>, |b1>}, define

\[ U_{SM}=|b_0\rangle\langle b_0|\otimes I_M+|b_1\rangle\langle b_1|\otimes X_M. \]

This is unitary because the source projectors are orthogonal and sum to identity, and I and X are unitary. With meter input |0>,

\[ |b_k\rangle|0\rangle_M\longmapsto |b_k\rangle|k\rangle_M. \]

Reading the meter in the computational basis realizes the projective measurement instrument.

## Implementation

- qag_measurement_dilation.py constructs the 4x4 unitary and applies it to a source register entangled with any existing environment tensor.
- The meter remains an explicit final tensor axis. Outcome slices are unnormalized conditional states; their squared norms are outcome probabilities.
- The implementation validates basis unitarity, source-axis dimension, and finite amplitudes.

The module does not assume that the apparatus is energetically free. Its meter is an abstract degenerate qubit with no calibrated Hamiltonian, spatial coordinate, or momentum. This is an explicit dilation of the measurement channel, not yet a mechanical detector model.

## Regression coverage

tests/test_qag_measurement_dilation.py verifies exact unitarity, Born-rule probabilities, retention of an environment entangled with the source, orthogonal meter outcome sectors, and rejection of malformed inputs.

## Remaining physical closure

A mechanical apparatus still requires a pointer coordinate and conjugate momentum, a specified interaction Hamiltonian, and a conserved total generator for source + probe + apparatus. The present dilation formalizes measurement information transfer but does not establish mechanical momentum or energy balance. Gravity-specific coupling and physical QAG/BQG observables remain separate unresolved gates.