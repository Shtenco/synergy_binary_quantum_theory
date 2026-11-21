# q=2 orientation reversal and complex conjugation share one label reflection

**Status:** exact common-covariance theorem on the frozen q=2 label carrier; dynamical identification of the geometry and phase sectors remains open.

This note sharpens a connection that previously looked only analogical:

```text
orientation reversal <-> complex conjugation.
```

The new result is that, on the frozen q=2 four-label carrier, one can choose a single explicit reflection `R` that reverses both:

1. the real two-dimensional complex-structure block `J` of the oriented C4 shift;
2. the oriented Walsh tetrahedral frame built from the same four labels.

Because the existing geometry-qubit result identifies oriented volume with the logical `Y_L` pseudoscalar, the same reflection has the common sign rule

\[
\boxed{
R:\quad J\mapsto-J,\qquad Q\mapsto-Q,\qquad Y_L\mapsto-Y_L.
}
\]

This is an exact representation/covariance result. It does **not** yet prove that one microscopic Hamiltonian dynamically locks the two sectors together.

---

# 1. The four q=2 labels

Take the four frozen route labels in Gray-cycle order:

```text
0: 00
1: 01
2: 11
3: 10
```

Adjacent entries differ in exactly one bit, including the closing edge `10 -> 00`, so this is the q=2 square `C4`.

Define the oriented cyclic shift

\[
S|k\rangle=|k+1\pmod4\rangle.
\]

Define the reflection

\[
R:k\mapsto-k\pmod4.
\]

Explicitly it fixes labels `00` and `11` and swaps

```text
01 <-> 10.
```

Therefore

\[
\boxed{RSR^{-1}=S^{-1}}.
\]

This is the ordinary dihedral relation for the square.

---

# 2. The phase representation

On the real Fourier pair

\[
u_c=(1,0,-1,0),
\]

\[
u_s=(0,1,0,-1),
\]

the shift acts as

\[
Su_c=u_s,
\]

\[
Su_s=-u_c.
\]

Hence in the basis `(u_c,u_s)`:

\[
S\mapsto J=
\begin{pmatrix}
0&-1\\
1&0
\end{pmatrix},
\qquad J^2=-I.
\]

The reflection acts as

\[
Ru_c=u_c,
\]

\[
Ru_s=-u_s,
\]

so

\[
R_2=
\begin{pmatrix}
1&0\\
0&-1
\end{pmatrix}.
\]

Therefore

\[
\boxed{R_2JR_2^{-1}=-J}.
\]

If `J` is denoted spectrally by `i`, then this is exactly

\[
\boxed{i\mapsto-i}.
\]

So complex conjugation is implemented by an ordinary reflection of the underlying real phase plane.

---

# 3. The same labels define a tetrahedral frame

The three nontrivial Walsh characters of

\[
\mathbb Z_2^2
\]

map each binary label to an integer vector in `R^3`:

```text
00 -> ( 1, 1, 1)
01 -> (-1, 1,-1)
11 -> (-1,-1, 1)
10 -> ( 1,-1,-1)
```

After normalization by `sqrt(3)` these are the four unit face normals of a regular tetrahedron.

Using the Gray-cycle ordering above, the oriented determinant of the three edge differences is

\[
\boxed{\det=-16}
\]

for the unnormalized integer frame.

Apply the **same reflection**

```text
01 <-> 10.
```

The determinant becomes

\[
\boxed{+16}.
\]

Therefore

\[
\boxed{\text{tetrahedral orientation sign flips}.}
\]

The reflection is an odd permutation of the four face labels.

---

# 4. Connection to the logical geometry-qubit Y channel

The already established exact geometry-qubit bridge has

\[
\boxed{
Q=J_1\cdot(J_2\times J_3)=\frac{\sqrt3}{4}Y_L.
}
\]

`Q` is an oriented triple product. An odd permutation of the face labels changes its sign:

\[
Q\mapsto-Q.
\]

Therefore

\[
\boxed{Y_L\mapsto-Y_L}.
\]

The same q=2 label reflection thus gives the sign table

| carrier | object | action of R |
|---|---|---|
| C4 phase plane | `J` | `J -> -J` |
| complex spectral notation | `i` | `i -> -i` |
| Walsh tetrahedron | oriented determinant | `det -> -det` |
| logical geometry qubit | `Q ~ Y_L` | `Y_L -> -Y_L` |

This is the central exact result.

---

# 5. What this does and does not mean

It **does** mean:

> The geometry orientation pseudoscalar and complex conjugation can be represented as the same sign character of one explicit q=2 label reflection.

It does **not** yet mean:

> The physical quantum phase of nature has been derived from the geometry-qubit `Y` operator.

Why not?

Because the two representations currently live in different declared roles:

```text
Walsh / face-qubit sector -> spatial geometry
oriented C4 spectral sector -> finite phase/complex structure.
```

A physical identification requires a **joint microscopic operator** or history amplitude showing that the same dynamical reflection acts on both and that their signs cannot be chosen independently.

---

# 6. New candidate microscopic principle

The new theorem suggests a very economical possibility:

```text
q=2 labels
-> one orientation character
   |-> geometry orientation Y
   |-> phase complex structure J
```

or symbolically

\[
\boxed{
\chi_{\rm orient}(R)=-1
\Rightarrow
\begin{cases}
Y_L\to-Y_L,\\
J\to-J.
\end{cases}}
\]

If future dynamics enforces this shared character, then complex conjugation would no longer be an independently inserted algebraic operation. It would be the phase-space image of the same microscopic orientation reversal already present in the spatial quantum geometry.

This is now a precise falsifiable target.

---

# 7. Relation to parity and time reversal

One must not jump directly from the finite reflection `R` to physical `P`, `T` or `CPT`.

The current theorem only concerns:

- a four-label q=2 reflection;
- tetrahedral orientation;
- the sign of a real complex-structure operator.

To identify it with spacetime parity/time reversal one would need the full Lorentzian/history representation and its antiunitary/unitary implementation.

So the legal ladder is

```text
q2 reflection
-> orientation sign / conjugation sign          EXACT
-> physical spatial parity                       OPEN
-> physical time reversal / antiunitarity        OPEN
-> CPT interpretation                            OPEN.
```

---

# 8. Executable gate

Run

```bash
python scripts/q2_orientation_complex_conjugation_gate.py
```

or

```bash
python scripts/q2_orientation_complex_conjugation_gate.py \
  --json verification_results/q2_orientation_complex_conjugation_gate.json
```

The gate checks with exact integer arithmetic:

- `R S R = S^-1`;
- the real Fourier `J` block;
- `R J R = -J`;
- exact Walsh tetrahedron integer vectors;
- determinant `-16 -> +16`;
- odd permutation parity;
- common sign `-1` for phase complex structure, tetrahedral orientation and logical-Y pseudoscalar covariance.

---

# 9. Strongest current statement

The q=2 arithmetic/geometry layer now supports the exact common-covariance chain

\[
\boxed{
\text{one binary-label reflection}
\to
\begin{cases}
\text{tetrahedral orientation reversal},\\
\text{complex-structure reversal}
\end{cases}
}
\]

with both represented by the same sign character.

The next killer gate is dynamical: construct one microscopic/history operator that acts simultaneously on the geometry and phase carriers and test whether the relative sign is forced rather than chosen.
