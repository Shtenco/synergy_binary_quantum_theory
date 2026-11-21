# Foundational audit: starting qubit versus emergent complex/history structure

Status: **logical scope clarification / foundational open fork. No existing finite geometry or constraint result is invalidated by this audit.**

## 1. The issue

The canonical BQG presentation starts with a quantum two-level carrier

\[
|\psi\rangle=\alpha|0\rangle+\beta|1\rangle,
\qquad
|\alpha|^2+|\beta|^2=1,
\]

with the usual complex Hilbert-space interpretation.

Later, independently interesting q=2 history/arithmetic results derive a real quarter-turn operator

\[
J^2=-I
\]

from the oriented `C4/C8` structure and use it to represent continuous phase rotations

\[
W(\theta)=e^{-\theta J}
\leftrightarrow e^{i\theta}.
\]

The same branch also identifies a unique positive quadratic phase-weight precursor under stated symmetry assumptions.

These statements are compatible only if their logical scope is stated carefully.

## 2. What may not be claimed simultaneously

If complex amplitudes and ordinary qubit Hilbert space are part of the microscopic starting ansatz, then the later derivation of `J`, complex phase and a quadratic norm **cannot be advertised as a derivation of complex quantum mechanics itself from pre-quantum binary distinctions**.

Doing so would be circular:

```text
assume complex qubit quantum mechanics
-> derive binary geometry/history
-> recover an internal complex structure
-> claim complex quantum mechanics was derived from binary structure.
```

The recovered `J` can still be a genuine structural result: it can identify why the q=2 orientation/history carrier naturally supports the same real algebra as multiplication by `i`. It simply does not remove the original quantum-mechanical assumption.

## 3. Two scientifically consistent programmes

### Programme A — quantum gravity from qubits

Declare explicitly:

```text
microscopic distinguishability + quantum two-level Hilbert carrier = starting ansatz.
```

Then the theory asks whether qubit relations generate

```text
3D geometry
-> constraints / GR structure
-> physical histories
-> gravity / gauge effective actions.
```

In this programme:

- superposition is inherited quantum kinematics;
- complex amplitudes are not derived;
- the Born rule is not derived;
- C4/C8 `J` is an **internal geometric/history realization of complex structure**, not the origin of complex numbers in quantum theory;
- light interference still requires a derived photon/Maxwell physical sector, but the quadratic quantum probability rule may be treated as part of the assumed quantum framework.

This is the smallest change to the current canonical project and is fully compatible with the existing SU(2)/Peter-Weyl machinery.

### Programme B — quantum mechanics from pre-quantum binary structure

If the stronger claim is desired, the microscopic starting point must be reformulated before complex Hilbert structure is used.

A possible research target would be

```text
real binary state/transition carrier
-> oriented C4/C8 algebra
-> J^2=-I
-> emergent complex scalar multiplication
-> positive inner product / norm
-> composition/tensor structure
-> unitary dynamics
-> measurement/probability rule
-> only then complex qubits and SU(2).
```

This is a much larger foundational programme. Existing SU(2), spin-network and Peter-Weyl results may serve as downstream targets, but they cannot be used upstream to prove the quantum axioms from which they were built.

## 4. Consequence for the interference claim

The observed two-path formula

\[
P
=|A_1+A_2|^2
=|A_1|^2+|A_2|^2
+2\operatorname{Re}(A_1A_2^*)
\]

requires at least

1. coherent linear amplitude composition;
2. complex/phase structure or an equivalent real two-dimensional representation;
3. a quadratic probability rule.

The current BQG repository has strong structural material for item 2 through the real `J` history/orientation algebra, and standard finite reference controls for coherent composition/decoherence/Sorkin `I3`.

But under Programme A, items 1 and 3 are inherited from ordinary quantum mechanics rather than derived from BQG.

Under Programme B, all three must be reconstructed without importing the complex qubit/Born framework at the start.

Therefore the strongest current fail-closed statement is:

\[
\boxed{
\text{BQG supplies a binary/history realization naturally compatible with quantum phase interference,}
}
\]

not yet

\[
\boxed{
\text{BQG derives the full quantum interference probability law from purely classical binary distinctions.}
}
\]

## 5. Photon physics remains a separate dynamical gate

Even a complete derivation of complex phase would not by itself derive light.

A physical photon requires the independent open gate

\[
\Gamma_{AA}^{(2)}(\omega,k)
\]

with

- gauge/Ward consistency;
- a massless transverse pole;
- positive physical residue;
- deconfined long-range phase;
- Maxwell stiffness `Z_A`;
- a common IR causal cone with the gravitational sector.

Thus the logical chain for optical interference is

```text
quantum amplitude/phase structure
+ physical photon sector
+ preparation/path/detection boundary amplitude
+ probability rule
-> observed interference.
```

No single C4 phase identity substitutes for the Maxwell and measurement parts of this chain.

## 6. Born-weight precursor status

A symmetry argument selecting a positive quadratic norm such as

\[
|z|^2=x^2+y^2
\]

is a useful precursor but is not a measurement theorem.

A full Born-rule derivation would additionally have to explain why outcome probabilities for arbitrary projective/POVM measurements are given by the Hilbert-space quadratic measure and why the rule composes consistently for subsystems.

Until that is done, the repository should retain the wording `Born-weight precursor` rather than `Born rule derived`.

## 7. Recommended canonical choice

For the current physicalization programme the scientifically economical choice is **Programme A**:

> BQG is a candidate quantum-gravity theory whose starting microscopic carriers are qubits; its q=2 history structure supplies a nontrivial internal geometric realization of orientation and complex phase, while a derivation of quantum mechanics itself is not currently claimed.

This preserves every current geometry/constraint result and removes the foundational circularity from the public claim surface.

Programme B can remain a separately labelled foundational research branch if desired.

## 8. New falsification/status rule

Until Programme B is completed, any statement of the form

```text
complex quantum mechanics / Born rule is derived from binary distinctions
```

must be rejected by the canonical audit.

Allowed statements include

```text
q=2 oriented history derives a real J with J^2=-I;
this supplies an internal representation of complex phase;
the current candidate gravity theory nevertheless starts from quantum two-level carriers.
```

This scope clarification is especially important before comparing BQG interference predictions with experiment.
