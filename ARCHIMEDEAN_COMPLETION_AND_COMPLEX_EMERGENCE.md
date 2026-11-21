# From modular arithmetic to ordinary real/complex arithmetic

**Status:** exact arithmetic controls + standard number-theory theorems; physical selection of the Archimedean place and scale remains open.

The previous two notes established

```text
finite complex/modular arithmetic
<-> ordinary modular matrix/scalar arithmetic
-> winding / CRT lift
-> Z
-> bounded rational reconstruction
-> Q.
```

This note asks the next question:

> Once we have ordinary rational arithmetic, what extra structure distinguishes the familiar real number line used by macroscopic physics from p-adic/profinite arithmetic?

The answer is very sharp.

---

# 1. Modular refinement and real completion are different directions

Congruences such as

\[
x\equiv y\pmod{p^n}
\]

say that `x-y` is divisible by a large power of `p`.

This naturally defines the p-adic notion of closeness

\[
|x-y|_p=p^{-v_p(x-y)}.
\]

Completing the rationals in this metric gives

\[
\boxed{\mathbb Q_p}.
\]

By contrast, ordinary physical magnitudes use the Archimedean absolute value

\[
|x|_\infty=|x|
\]

and its completion gives

\[
\boxed{\mathbb R}.
\]

Thus

```text
more and more modular precision -> p-adic / profinite direction
more and more Archimedean metric precision -> real direction.
```

They are not the same continuum limit.

---

# 2. Ostrowski's theorem identifies the missing choice

A standard theorem of number theory states that every nontrivial absolute value on `Q`, up to equivalence, is one of:

\[
|\cdot|_\infty
\]

or

\[
|\cdot|_p
\]

for a prime `p`.

Therefore once the arithmetic layer has reconstructed

\[
\mathbb Q,
\]

the possible completions are not an arbitrary zoo.

They split into

```text
one Archimedean place      -> R
one finite place per prime -> Q_p.
```

This makes the role of the ordinary continuum explicit:

\[
\boxed{
\text{ordinary real arithmetic}
=\text{rational arithmetic}
+\text{Archimedean size/topology}
+\text{completion}.
}
\]

Finite modular arithmetic alone cannot choose the Archimedean place because every finite ring has no compatible translation-invariant total order.

---

# 3. An exact finite illustration: the same rational sequence behaves oppositely

Take

\[
x_n=1+p+p^2+\cdots+p^{n-1}.
\]

In ordinary arithmetic

\[
x_n\to+\infty.
\]

Indeed the increment is

\[
x_{n+1}-x_n=p^n,
\]

which becomes **larger and larger** in the ordinary absolute value.

But p-adically

\[
|x_{n+1}-x_n|_p=p^{-n}\to0.
\]

So the same sequence is p-adic Cauchy and converges to

\[
\boxed{\frac1{1-p}}
\]

inside `Q_p`.

This is a very concrete warning:

> `continuum limit` has no meaning until the notion of distance/size has been declared or derived.

The executable gate checks this exactly for `p=2,3,5,7` without floating point.

---

# 4. A finite control of Q -> R completion

To illustrate the Archimedean completion with exact rational arithmetic, define

\[
a_n=\left\lfloor2^n\sqrt2\right\rfloor.
\]

Then

\[
\frac{a_n}{2^n}
\le\sqrt2
<\frac{a_n+1}{2^n}
\]

and the width is exactly

\[
\boxed{2^{-n}}.
\]

These nested rational intervals shrink to one Archimedean real point.

The gate constructs `a_n` using integer square roots only and verifies the nested brackets through `n=24`.

This is a finite certificate of the mechanism

\[
\boxed{
\mathbb Q
\xrightarrow{\text{Archimedean Cauchy completion}}
\mathbb R.
}
\]

The theorem that the full completion is `R` is standard mathematics; the gate is a reproducible exact positive control, not a replacement for the theorem.

---

# 5. The product formula links the finite and infinite notions of size

For every nonzero rational number,

\[
\boxed{
|x|_\infty\prod_p|x|_p=1.
}
\]

Only finitely many factors differ from one for a fixed rational `x`.

So the Archimedean and p-adic places are not unrelated decorations: they obey one exact global constraint.

For example, prime factors that make a rational number p-adically small are compensated globally by its Archimedean and other-place sizes.

The executable gate verifies the product formula exactly with rational arithmetic for several nontrivial fractions.

This suggests a useful architecture for the information-graph theory:

```text
finite modular data -> finite places / divisibility information
history + ordinary lift -> Q
Archimedean place -> macroscopic ordered magnitude
all places together -> global arithmetic consistency.
```

No claim is made that a physical spacetime degree of freedom literally stores all adeles or all p-adic numbers.

---

# 6. The real numbers do not require primitive irrational symbols

Once `Q` and the Archimedean metric are present, irrational numbers can appear as limits of rational information.

Thus a microscopic theory does not need to store `sqrt(2)` or `pi` as primitive infinite strings.

It can instead provide increasingly accurate rational observables whose Cauchy class is the real number.

Conceptually:

\[
\boxed{
\text{finite arithmetic at every finite resolution}
+\text{compatible refinement}
\to\text{real continuum value}.
}
\]

This is closely analogous to the repository's spacetime coarse-graining picture:

```text
finite/discrete description at each scale
+ compatible refinement/convergence
-> smooth continuum object.
```

The analogy is structural, not a proof that the same RG map controls both geometry and arithmetic.

---

# 7. Complex numbers can then emerge from a real oriented two-plane

Recall the already-derived ordinary real/integer operator

\[
J=
\begin{pmatrix}
0&-1\\
1&0
\end{pmatrix},
\qquad
J^2=-I.
\]

Over the rationals define

\[
\boxed{
\mathbb Q[J]=\{aI+bJ:a,b\in\mathbb Q\}.
}
\]

Multiplication gives

\[
(aI+bJ)(cI+dJ)
=(ac-bd)I+(ad+bc)J.
\]

This is exactly Gaussian rational arithmetic:

\[
\boxed{\mathbb Q[J]\cong\mathbb Q(i)}.
\]

Now perform the same Archimedean completion on the two real coefficients. The result is

\[
\boxed{
\mathbb R[J]\cong\mathbb C.
}
\]

Therefore one mathematically clean emergence route is

\[
\boxed{
\mathbb Z_M
\to\mathbb Z
\to\mathbb Q
\to\mathbb R
\xrightarrow{J^2=-I}
\mathbb C.
}
\]

The complex field need not be a primitive microscopic datatype. It may be the compact/spectral language of an oriented real two-dimensional structure after the ordinary-number completion has been selected.

---

# 8. Or one may introduce J before completion — the diagram commutes

There is another order:

\[
\mathbb Q\to\mathbb Q(i)\to\mathbb C.
\]

Because `J` has rational/integer matrix entries, introducing the complex structure before or after Archimedean coefficient completion yields the same final algebra:

```text
Q -----------------> R
|                     |
add J                 add J
|                     |
v                     v
Q(i) --------------> C
```

This is useful physically: the orientation/phase structure and the ordinary continuum scale can be conceptually separated.

The executable gate exhaustively checks the `Q[J]` addition/multiplication identities on a finite rational test grid.

---

# 9. Complex conjugation remains orientation reversal

On the real two-plane,

\[
J\to-J
\]

under reversal of the chosen orientation.

Therefore

\[
aI+bJ\to aI-bJ
\]

which is exactly

\[
a+bi\to a-bi.
\]

Hence the candidate relation becomes

\[
\boxed{
\text{orientation reversal}
\longleftrightarrow
\text{complex conjugation}.
}
\]

The algebraic statement is exact.

The theory-specific claim that the **same** q=2 geometry orientation operator dynamically selects this `J` is still open and requires a joint operator gate.

---

# 10. Why the one common physical scale is unavoidable, not embarrassing

After winding/history reconstruction we may have an ordinary integer count

\[
N\in\mathbb Z.
\]

Suppose a physical quantity `Q` is additive under concatenation:

\[
Q(N+M)=Q(N)+Q(M).
\]

On the integers this alone implies

\[
\boxed{Q(N)=sN}
\]

for one constant

\[
s=Q(1).
\]

So a free global conversion factor is exactly what arithmetic predicts when a dimensionless count is promoted to a dimensionful observable.

Arithmetic can determine

```text
how counts add
how residues compose
how winding unwraps
how ratios reconstruct
```

but it cannot determine whether one elementary count means

```text
one metre
one second
one joule-second
one Planck unit
```

without additional physics.

This explains why the existing physicalization programme repeatedly arrives at **one common scale**, not six independent directional scales.

The gate checks the integer additivity result on several rational scale choices as a finite regression control; the theorem itself is immediate from repeated addition.

---

# 11. A possible global arithmetic picture: finite places + one infinite place

A useful but still interpretive summary is

\[
\boxed{
\text{microscopic finite arithmetic}
\to\text{congruence data at finite places}
\to\mathbb Q
\to
\begin{cases}
\mathbb Q_p & \text{finite-place completions},\\
\mathbb R & \text{infinite/Archimedean completion}
\end{cases}
}
\]

and then, if an oriented two-plane/quarter-turn is selected,

\[
\boxed{\mathbb R\to\mathbb C.}
\]

This is stronger than saying "modular complex numbers can be converted into normal numbers".

The more precise statement is:

> Modular arithmetic supplies exact finite local arithmetic. Winding/CRT can lift bounded compatible data to ordinary integers and rationals. Ordinary real arithmetic requires the Archimedean topology and completion. Complex arithmetic can then be generated by an oriented real quarter-turn operator `J` with `J^2=-I`.

---

# 12. What this could mean for the information-graph theory

A new candidate chain is now visible:

\[
\boxed{
\begin{aligned}
&\text{binary q=2 distinctions}\\
&\to \mathbb Z_2^2\text{ Walsh geometry}\\
&\to C_4\text{ oriented transport}\\
&\to J^2=-I\\
&\to\text{finite phase spectrum}\\
&\leftrightarrow\text{modular integer arithmetic}\\
&\to\text{history/winding}\\
&\to\mathbb Z\\
&\to\mathbb Q\\
&\to\text{Archimedean completion}\\
&\to\mathbb R\\
&\xrightarrow{J}\mathbb C.
\end{aligned}}
\]

There are now three theory-specific killer questions.

### A. Why the Archimedean place?

The mathematics classifies the possible completions, but the microscopic dynamics must explain why macroscopic rods/clocks use the Archimedean one rather than a p-adic metric.

### B. Is J the same orientation degree already present in geometry?

Need a joint q=2 operator proof connecting logical geometric orientation `Y` to cyclic transport orientation `J`.

### C. Does phase/history refinement produce compatible winding?

Need an exact recursive refinement law, not an imposed unwrapping convention.

Until these are passed, this is an arithmetic architecture, not yet a physical derivation of quantum complex amplitudes.

---

# 13. New status ledger

| Statement | Status |
|---|---|
| modular congruence refinement is non-Archimedean in character | **STANDARD MATH** |
| nontrivial absolute values on Q are Archimedean or p-adic up to equivalence | **OSTROWSKI THEOREM** |
| Archimedean completion of Q is R | **STANDARD MATH** |
| exact nested rational brackets converge to sqrt(2) | **EXACT FINITE CONTROL** |
| same rational sequence can diverge real-wise and converge p-adically | **EXACT** |
| rational product formula | **STANDARD THEOREM + EXACT FINITE CONTROL** |
| `Q[J] ~= Q(i)` for `J^2=-I` | **EXACT** |
| `R[J] ~= C` | **EXACT STANDARD ALGEBRA** |
| integer additive physical map has form `Q(N)=sN` | **EXACT** |
| microscopic dynamics selects Archimedean physical place | **OPEN PHYSICAL** |
| geometry orientation `Y` is the same microscopic selector as `J` | **OPEN** |
| physical value of common scale `s` | **OPEN PHYSICAL** |

---

# 14. Reproduction

```bash
python scripts/archimedean_completion_complex_emergence_gate.py
```

or

```bash
python scripts/archimedean_completion_complex_emergence_gate.py \
  --json verification_results/archimedean_completion_complex_emergence_gate.json
```

The executable uses exact integer/Fraction arithmetic only.
