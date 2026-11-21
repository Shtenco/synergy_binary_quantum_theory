# Why a phase plane naturally produces |z|^2

**Status:** exact uniqueness theorem for positive quadratic phase-invariant weights; not a complete derivation of the quantum Born measurement rule.

The modular/arithmetic programme now has a real two-dimensional phase carrier with

\[
J=
\begin{pmatrix}
0&-1\\
1&0
\end{pmatrix},
\qquad J^2=-I.
\]

After the Archimedean lift, a complex amplitude can be written equivalently as

\[
z=a+bi
\]

or as the real vector

\[
v=(a,b).
\]

This note asks:

> If an observable weight is quadratic in the amplitude and must not depend on the arbitrary quarter-turn phase orientation, what scalar is allowed?

The answer is unique up to one overall factor.

---

# 1. General quadratic weight

Let

\[
Q(v)=v^TAv
\]

with a real symmetric matrix

\[
A=
\begin{pmatrix}
\alpha&\beta\\
\beta&\gamma
\end{pmatrix}.
\]

Require invariance under the exact q=2 quarter-turn phase operation:

\[
Q(Jv)=Q(v)
\]

for every `v`.

This is equivalent to

\[
\boxed{J^TAJ=A}.
\]

Direct multiplication gives

\[
J^TAJ=
\begin{pmatrix}
\gamma&-\beta\\
-\beta&\alpha
\end{pmatrix}.
\]

Therefore

\[
\alpha=\gamma,
\qquad
\beta=-\beta.
\]

Over `R` or `Q`, characteristic is not two, so

\[
\beta=0.
\]

Hence

\[
\boxed{A=\lambda I}.
\]

and

\[
\boxed{Q(a,b)=\lambda(a^2+b^2)}.
\]

In complex notation:

\[
\boxed{Q(z)=\lambda|z|^2}.
\]

Thus the squared modulus is not an arbitrary decorative choice inside the class of quadratic phase-invariant weights: it is the unique one up to normalization.

---

# 2. Positivity fixes the sign

If `Q` is to be a nonnegative weight,

\[
Q(v)\ge0
\]

for all `v`, then

\[
\boxed{\lambda\ge0}.
\]

If we also normalize one unit amplitude by

\[
Q(1,0)=1,
\]

then

\[
\boxed{\lambda=1}.
\]

Therefore

\[
\boxed{Q(z)=|z|^2}.
\]

Again: this is a theorem about a **quadratic phase-invariant positive scalar**. Calling it the full Born rule would require additional measurement/composition assumptions.

---

# 3. Why the Archimedean step matters

In a finite Gaussian modular ring,

\[
z\bar z=a^2+b^2\pmod N
\]

is an exact algebraic norm residue.

But a residue has no intrinsic ordinary order:

```text
3 mod 5 is not intrinsically "more positive" than 4 mod 5.
```

So at the finite modular stage one has

\[
\boxed{\text{norm algebra}}
\]

but not yet

\[
\boxed{\text{positive real weight}}.
\]

After the ordinary rational/Archimedean lift,

\[
a,b\in\mathbb R,
\]

and now

\[
\boxed{a^2+b^2\ge0}.
\]

Therefore the arithmetic ladder separates two logically different ingredients:

```text
finite cyclic/complex structure -> tells us which quadratic invariant exists
Archimedean ordered continuum    -> lets that invariant become a positive magnitude.
```

This is a potentially important conceptual bridge between finite phase arithmetic and ordinary probabilities/intensities.

---

# 4. Interference appears automatically from the same quadratic scalar

For two amplitude vectors

\[
v=(a,b),\qquad w=(c,d),
\]

we have

\[
Q(v+w)=\lambda\|v+w\|^2.
\]

Expanding,

\[
\boxed{
Q(v+w)=Q(v)+Q(w)+2\lambda\,v\cdot w.
}
\]

In complex notation

\[
v\cdot w=\operatorname{Re}(z\bar w),
\]

so

\[
\boxed{
|z+w|^2
=|z|^2+|w|^2+2\operatorname{Re}(z\bar w).
}
\]

The familiar interference cross-term is therefore simply the polarization identity of the unique quadratic phase-invariant real norm.

This does not yet tell us which alternatives in a physical experiment must have their amplitudes added before weighting. That is part of the measurement/history interpretation, not the algebra alone.

---

# 5. Relation to the matrix decomplexification theorem

Recall

\[
\Phi(z)=
\begin{pmatrix}
a&-b\\
b&a
\end{pmatrix}.
\]

Then

\[
\boxed{
\Phi(z)^T\Phi(z)
=(a^2+b^2)I.
}
\]

and

\[
\boxed{
det\,\Phi(z)=a^2+b^2}.
\]

Thus three viewpoints coincide:

```text
complex squared modulus |z|^2
real Euclidean norm     a^2+b^2
matrix determinant      det Phi(z).
```

At the modular level the same identity holds modulo `N`; after Archimedean lifting the scalar becomes an ordinary nonnegative magnitude.

---

# 6. Why only a quarter-turn is already enough

One might expect that full continuous `U(1)` symmetry is needed to force the Euclidean norm.

Surprisingly, for a symmetric quadratic form in two dimensions, invariance under the single discrete generator

\[
J^2=-I
\]

already gives

\[
A=\lambda I.
\]

So the finite q=2 `C4` phase symmetry is already strong enough to determine the quadratic metric on the amplitude plane up to scale.

This is especially relevant to the current programme because full physical `U(1)` dynamics has not yet been derived from finite refinement, while the order-four phase block is exact.

---

# 7. What would still be needed for a real Born-rule derivation

The current result does **not** prove all of quantum measurement theory.

To reach a physical probability rule one would still need to derive or justify, in the model:

1. which physical alternatives correspond to linear amplitude components;
2. when amplitudes add coherently before weighting;
3. when weights add incoherently for exclusive alternatives;
4. normalization across a complete outcome set;
5. compatibility with tensor products/composite systems;
6. the physical state/projector/history interpretation.

The legal statement today is therefore

\[
\boxed{
\text{q=2 phase symmetry}
+\text{quadraticity}
+\text{positivity}
+\text{normalization}
\Rightarrow |z|^2\text{ as the unique local scalar weight}.
}
\]

This is a **Born-weight precursor**, not the complete Born rule.

---

# 8. Why this is interesting for the information-graph programme

The arithmetic story now contains a remarkably coherent sequence:

\[
\boxed{
\begin{aligned}
&\text{q=2 C4 oriented shift}\\
&\to J^2=-I\\
&\to\text{complex phase representation}\\
&\to\text{ordinary modular/integer arithmetic}\\
&\to\mathbb Q\\
&\to\text{Archimedean completion }\mathbb R\\
&\to\mathbb R[J]\cong\mathbb C\\
&\to\text{unique quadratic phase scalar }|z|^2.
\end{aligned}}
\]

This suggests that neither complex numbers nor the squared norm need to be primitive microscopic instructions. Both can arise from a real finite orientation/rotation structure plus the ordinary-number continuum lift.

The next physical question is whether the actual history/projector dynamics uses this scalar as its observable/measure weight.

---

# 9. Executable gate

Run

```bash
python scripts/phase_invariant_quadratic_weight_gate.py
```

or

```bash
python scripts/phase_invariant_quadratic_weight_gate.py \
  --json verification_results/phase_invariant_quadratic_weight_gate.json
```

It checks exactly:

- the coefficient equations `J^T A J=A`;
- exhaustive rational-grid uniqueness of `A=lambda I`;
- phase invariance under the four `mu4` rotations;
- positivity/normalization;
- the exact interference polarization identity.
