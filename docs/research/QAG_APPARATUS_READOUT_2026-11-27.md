# QAG apparatus readout dilation — 2026-11-27

## Added operation

The relational interaction model produces a joint source/apparatus position
state \(\Psi(x_s,x_a)\). This increment adds an explicit N-level meter
initialized in \(|0\rangle_M\), with a modular controlled-addition unitary
\[
U_{AM}|x_a,m\rangle=|x_a,m+x_a\pmod N\rangle.
\]
For the initialized meter,
\[
\sum_{x_s,x_a}\Psi(x_s,x_a)|x_s,x_a,0\rangle
\longmapsto
\sum_{x_s,x_a}\Psi(x_s,x_a)|x_s,x_a,x_a\rangle.
\]
The meter outcome therefore records the apparatus position. Conditioning on
outcome \(m\) returns the normalized source-apparatus slice at
\(x_a=m\), with its Born probability. All N outcomes are retained and their
probabilities sum to one.

## What this closes

The previous apparatus model generated source/apparatus entanglement but did
not explicitly model the transfer of the apparatus position record to a
separate readout register. This dilation now supplies that finite-dimensional
measurement-information transfer and conditional readout.

## What remains open

- The meter is an abstract N-level register with no calibrated energy
  spectrum, spatial support, mass, or amplification dynamics.
- The readout is an ideal projective measurement after unitary premeasurement;
  environmental decoherence and robust classical records are not modeled.
- No physical detector support, external frame, switching work, or full
  mechanical energy budget is included.
- The finite ring has modular translations, not an unbounded canonical
  momentum. The readout operation is not evidence of a laboratory momentum
  balance or gravitational effect.
- This does not establish antigravity, propulsion, or a new gravitational
  interaction.

The next closure should model the meter's free Hamiltonian, measurement
switching, environmental decoherence, and the conserved generator of the
complete source-apparatus-meter-support system. It must then compare the
full unconditional ensemble with the postselected weak-value description.

## Symmetry audit: the readout imports an origin

The controlled-addition meter starts in \(|0\rangle_M\) and records the
absolute apparatus site label. This is not invariant under a simultaneous
translation of apparatus and meter:
\[
[U_{AM},T_A\otimes T_M]\ne0.
\]
The meter's initialized zero state therefore supplies an implicit spatial
reference frame. This is not a numerical defect; it is a physical resource
hidden in the ideal readout convention. A translation-covariant readout would
need an explicit reference-frame/support register and a relational record,
rather than silently treating the label zero as an absolute origin. The
regression suite now tests and preserves this non-covariance as an explicit
limitation.
