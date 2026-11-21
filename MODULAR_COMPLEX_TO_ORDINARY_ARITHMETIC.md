# Modular-complex numbers -> ordinary arithmetic

**Status:** exact arithmetic/representation bridge; physical interpretation beyond the finite phase layer is conditional.

This note explores a new direction for the information-graph programme: can a finite modular phase/complex arithmetic be rewritten in terms of ordinary arithmetic, and can that help explain why complex quantum amplitudes and ordinary real-valued observables coexist?

The answer is surprisingly strong at the algebraic level:

```text
finite cyclic/oriented combinatorics
-> modular phase
-> complex spectral representation
-> ordinary integer matrices / scalar residue channels
-> CRT bounded lift to ordinary integers
-> scale + continuum limit required for ordered real observables.
```

The important boundary is equally strong: **a single finite modular ring does not by itself contain the ordinary Archimedean order of the real numbers.** The finite-to-real bridge therefore needs an additional lift/scale/limit, not a relabelling trick.

---

## 1. What does "complex modulo N" mean?

Writing `C mod N` is ambiguous because the ordinary complex field is not naturally reduced modulo an integer in the way `Z` is.

The clean finite object used here is the Gaussian modular ring

\[
\boxed{
R_N=\mathbb Z[i]/N\mathbb Z[i]
\cong (\mathbb Z/N\mathbb Z)[x]/(x^2+1).
}
\]

Every element is represented by a pair

\[
z=a+bi,\qquad a,b\in\mathbb Z/N\mathbb Z,
\]

with

\[
i^2=-1.
\]

Addition and multiplication are

\[
(a,b)+(c,d)=(a+c,b+d),
\]

\[
(a,b)(c,d)=(ac-bd,ad+bc)
\]

with every component reduced modulo `N`.

---

# Theorem A. Complex modular arithmetic is ordinary matrix arithmetic

Define

\[
\boxed{
\Phi_N(a+bi)=
\begin{pmatrix}
a&-b\\
b&a
\end{pmatrix}
\pmod N.
}
\]

Then

\[
\Phi_N(z+w)=\Phi_N(z)+\Phi_N(w),
\]

and

\[
\boxed{\Phi_N(zw)=\Phi_N(z)\Phi_N(w)}.
\]

Therefore

\[
\boxed{
\mathbb Z[i]/N\mathbb Z[i]
\hookrightarrow
M_2(\mathbb Z/N\mathbb Z)
}
\]

as the subring of matrices of the form

\[
\begin{pmatrix}a&-b\\b&a\end{pmatrix}.
\]

Nothing complex is needed to execute the arithmetic: all entries are ordinary integers modulo `N`.

The familiar complex structures become ordinary matrix invariants:

\[
\boxed{\det\Phi_N(z)=a^2+b^2=N(z)\pmod N},
\]

\[
\boxed{\operatorname{tr}\Phi_N(z)=2a\pmod N},
\]

and conjugation is transpose:

\[
\boxed{\Phi_N(\bar z)=\Phi_N(z)^T}.
\]

This is the first exact decomplexification theorem.

---

## 2. The imaginary unit is a quarter-turn operator

Let

\[
J=
\begin{pmatrix}
0&-1\\
1&0
\end{pmatrix}.
\]

Then

\[
\boxed{J^2=-I}.
\]

Hence

\[
a+bi
\quad\leftrightarrow\quad
aI+bJ.
\]

The symbol `i` can therefore be interpreted as the spectral shorthand for an ordinary two-dimensional quarter-turn operator.

This is not merely analogy: composition of the ordinary matrices reproduces complex multiplication exactly.

---

# Theorem B. The q=2 square already contains the real block J

The frozen q=2 route adjacency is the square

\[
Q_2=C_4.
\]

Important distinction:

\[
\boxed{\mathbb Z_2^2\not\cong\mathbb Z_4}.
\]

The four microscopic labels `00,01,10,11` form the group `Z2^2`; we must not silently replace that group by `Z4`.

However, the **oriented cycle automorphism** of the square has order four.

Let `S` be the integer permutation operator

\[
S|k\rangle=|k+1\pmod4\rangle.
\]

In the ordinary four-state basis,

\[
S=
\begin{pmatrix}
0&0&0&1\\
1&0&0&0\\
0&1&0&0\\
0&0&1&0
\end{pmatrix},
\qquad S^4=I.
\]

Choose the real vectors

\[
u_0=(1,1,1,1),
\]

\[
u_2=(1,-1,1,-1),
\]

\[
u_c=(1,0,-1,0),
\]

\[
u_s=(0,1,0,-1).
\]

Then

\[
Su_0=u_0,
\]

\[
Su_2=-u_2,
\]

\[
Su_c=u_s,
\]

\[
Su_s=-u_c.
\]

Therefore in this real basis

\[
\boxed{
S\sim
1\oplus(-1)\oplus
\begin{pmatrix}0&-1\\1&0\end{pmatrix}.
}
\]

The two-dimensional block is exactly `J`.

After complexification its two eigenvalues are

\[
\boxed{+i,-i}.
\]

So the complex unit appears as the spectral coordinate of an ordinary integer permutation dynamics on an oriented four-cycle.

### Orientation reversal becomes complex conjugation

Reversing the cycle gives

\[
S\to S^{-1}.
\]

On the two-dimensional block,

\[
J\to J^{-1}=-J.
\]

Thus

\[
\boxed{i\leftrightarrow-i}
\]

under orientation reversal.

This dovetails naturally with the existing geometry-qubit orientation pseudoscalar, but the stronger physical identification is **not yet claimed**: the theorem is representation-theoretic.

---

# Theorem C. Roots of unity turn complex multiplication into ordinary modular addition

Let

\[
\omega_M=e^{2\pi i/M}.
\]

The finite phase set

\[
\mu_M=\{1,\omega_M,\ldots,\omega_M^{M-1}\}
\]

obeys

\[
\omega_M^a\omega_M^b
=\omega_M^{a+b\;({\rm mod}\;M)}.
\]

Hence

\[
\boxed{\mu_M\cong\mathbb Z_M}
\]

where **complex multiplication becomes ordinary addition modulo `M`**.

For `M=4`,

```text
0 <->  1
1 <->  i
2 <-> -1
3 <-> -i
```

and

```text
phase multiplication <-> integer addition mod 4.
```

This is the finite version of phase unwrapping.

---

# Theorem D. For p = 1 mod 4 a modular complex number splits into two ordinary scalar channels

Let `p` be an odd prime with

\[
p\equiv1\pmod4.
\]

Then `-1` has two square roots modulo `p`:

\[
r^2\equiv-1\pmod p.
\]

Therefore

\[
x^2+1=(x-r)(x+r)
\]

over `F_p`, and the Gaussian modular ring splits:

\[
\boxed{
\mathbb F_p[i]
\cong
\mathbb F_p\times\mathbb F_p.
}
\]

The explicit map is

\[
\boxed{
\Psi_r(a+bi)=
(a+rb,\;a-rb)\pmod p.
}
\]

The inverse is

\[
\boxed{
a=\frac{u+v}{2},
\qquad
b=\frac{u-v}{2r}
\pmod p.
}
\]

Thus, for these primes, the **entire** modular complex arithmetic is exactly two copies of ordinary scalar modular arithmetic.

### Smallest example: p=5

Since

\[
2^2=4\equiv-1\pmod5,
\]

we may use

\[
i\mapsto2
\]

in one scalar channel and

\[
i\mapsto-2\equiv3
\]

in the conjugate channel.

Thus

\[
\boxed{
a+bi\mapsto(a+2b,\;a-2b)\pmod5}.
\]

The four unit phases map multiplicatively as

```text
1   -> 1
i   -> 2
-1  -> 4
-i  -> 3
```

because

\[
2^4\equiv1\pmod5.
\]

This is an exact finite example of complex phase arithmetic represented only by ordinary integers.

---

# Theorem E. For p = 3 mod 4 the complex sector does not split

For an odd prime

\[
p\equiv3\pmod4,
\]

`-1` is not a quadratic residue. Therefore `x^2+1` is irreducible over `F_p`, and

\[
\boxed{
\mathbb F_p[i]\cong\mathbb F_{p^2}.
}
\]

In this case there is no decomposition into two scalar `F_p` channels.

Nevertheless Theorem A still applies: the same arithmetic is represented by ordinary `2x2` modular matrices.

Every nonzero `a+bi` has inverse

\[
\boxed{
(a+bi)^{-1}
=\frac{a-bi}{a^2+b^2}
}
\]

because `a^2+b^2` cannot vanish modulo `p` unless `a=b=0`.

So there are two arithmetic universality classes:

```text
p = 1 mod 4 -> split pair of scalar residue channels
p = 3 mod 4 -> genuine quadratic extension
```

This classification is mathematical. The project does **not** yet claim that nature selects a physical prime modulus.

---

## 3. The p=2 exception is ramified

Modulo two,

\[
x^2+1=(x+1)^2.
\]

Let

\[
\epsilon=1+i.
\]

Then

\[
\boxed{\epsilon\neq0,\qquad\epsilon^2=0\pmod2}.
\]

Thus the binary prime is special: instead of a split pair or a quadratic field, it contains a nilpotent direction.

This is potentially interesting for a binary microscopic theory, but no physical use is asserted yet. It is registered as an arithmetic structural fact to be investigated.

---

# Theorem F. Many modular channels can reconstruct ordinary integers exactly

Suppose we compute an integer `x` modulo pairwise-coprime moduli

\[
m_1,m_2,\ldots,m_k.
\]

Let

\[
M=\prod_jm_j.
\]

The Chinese remainder theorem reconstructs a unique residue

\[
x\pmod M.
\]

If physics/computation supplies a bound

\[
|x|<M/2,
\]

then the centered representative is the unique ordinary integer `x`.

For Gaussian integers the construction is componentwise:

\[
\boxed{
(a\bmod m_j,\;b\bmod m_j)_{j=1}^k
\xrightarrow{CRT}
a+bi
}
\]

whenever the true components lie inside the centered window.

The executable gate uses

```text
moduli = 5,7,11
M = 385
z = 7 + 11 i
w = -5 + 6 i
```

and reconstructs exactly

\[
z+w=2+17i,
\]

\[
zw=-101-13i.
\]

This is an exact route from several finite modular computations to ordinary integer arithmetic.

---

# No-go 1. One finite modulus cannot create the ordinary number line

A nontrivial finite cyclic additive group cannot carry a translation-invariant total order compatible with addition.

Assume `1>0` in `Z_N`. Then repeated addition gives

\[
1>0,\quad2>1,\quad\ldots,\quad N> N-1.
\]

But

\[
N\equiv0\pmod N,
\]

which would imply

\[
0> N-1>\cdots>0,
\]

a contradiction. The `1<0` case fails similarly.

Therefore

\[
\boxed{
\text{finite modular arithmetic}\not\Rightarrow
\text{ordinary ordered arithmetic}
}
\]

without extra structure.

This is important for the physics programme: metric lengths, energies and probabilities eventually require an ordered/Archimedean interpretation.

---

# No-go 2. A residue has no unique ordinary lift without a window or winding register

From

\[
x\equiv r\pmod N
\]

one only knows

\[
x=r+kN,
\qquad k\in\mathbb Z.
\]

A unique ordinary integer requires either:

1. an a priori bounded window, or
2. additional winding/history data telling us the integer `k`.

The same statement appears for phase:

\[
e^{i\theta}
\]

determines

\[
\theta\pmod{2\pi},
\]

not a unique real `theta`.

A global ordinary phase requires an **unwrapped history** plus integer winding.

This suggests a natural architecture:

```text
compact phase
+ history/winding register
-> ordinary additive phase on the universal cover.
```

That is a potentially useful connection to the repository's existing history/projector layer.

---

# 4. A new arithmetic emergence ladder

The strongest mathematically legal chain suggested by the calculations is

\[
\boxed{
\text{finite relation/order}
\to C_M\text{ shift}
\to\mu_M\text{ phase spectrum}
\leftrightarrow\mathbb Z_M
\to\text{multi-modular channels}
\to\text{CRT bounded integer lift}
\to\text{scaled rational observables}
\to\text{continuum completion}.
}
\]

A possible complex branch is

\[
\boxed{
\mathbb Z[i]/N
\to aI+bJ
\to\text{ordinary matrix arithmetic}
\to\text{CRT lift}
\to\mathbb Z[i].
}
\]

Only after scale normalization and a controlled limit may one discuss

\[
\mathbb Q,\quad\mathbb R,\quad\mathbb C
\]

as physical continuum number systems.

This is a much safer statement than claiming `the reals are literally one finite modular ring`.

---

# 5. How this may connect to q=2 quantum geometry

The new bridge creates several precise research questions.

## 5.1 Oriented C4 -> complex phase

The q=2 route graph is `C4`; its oriented shift has the exact `J^2=-I` block.

Question:

> Is the orientation needed to select `S` versus `S^{-1}` dynamically the same orientation channel already measured by the logical geometry-qubit pseudoscalar `Y`?

If yes, then

```text
orientation sign
-> complex conjugation sign
```

would become a derived bridge rather than a notation choice.

This is **OPEN**.

## 5.2 C4 -> U(1) under refinement

A single oriented four-cycle has only the discrete phase group `mu_4`.

A continuum U(1) requires a refinement sequence in which the available phase angles become dense while compatible holonomies converge.

A legal target is therefore

\[
\boxed{
\mu_{M_g}\to U(1),\qquad M_g\to\infty,
}
\]

with a refinement-compatible embedding of finite phase groups.

The existing Hopf/Pancharatnam U(1) bridge supplies the continuum target; this new work suggests a finite arithmetic precursor.

The dynamical/refinement derivation remains **OPEN**.

## 5.3 Modular phase -> physical real angle

A compact phase only gives an angle modulo `2pi`. The physical history may carry the winding data needed to unwrap it.

Candidate target:

\[
\boxed{
\text{compact link phase}
+\text{relational history}
\to\text{unwrapped additive phase}
\to\text{frequency/action observable}.
}
\]

This could become relevant to the existing constraint-to-history-to-physical-frequency problem, but no identification is frozen yet.

---

# 6. A surprising computational analogy: FFT -> NTT

The ordinary discrete Fourier transform diagonalizes cyclic shifts using complex roots of unity.

A number-theoretic transform uses roots of unity inside a finite modular ring instead, replacing floating complex arithmetic by exact integer modular arithmetic.

This is not proof about nature, but it demonstrates a concrete fact useful for the theory:

\[
\boxed{
\text{the same cyclic spectral algebra can have both a complex representation and an exact modular-integer representation.}
}
\]

Therefore the appearance of complex phases in an effective description does not logically require complex numbers to be primitive microscopic data.

---

# 7. Current theorem ledger

| Statement | Status |
|---|---|
| `a+bi mod N -> [[a,-b],[b,a]] mod N` ring embedding | **EXACT** |
| determinant = modular norm, transpose = conjugation | **EXACT** |
| oriented q=2 `C4` shift contains `J^2=-I` real block | **EXACT** |
| orientation reversal sends `J -> -J` | **EXACT** |
| `mu_M` multiplication = `Z_M` addition | **EXACT** |
| `p=1 mod4` Gaussian ring splits to `F_p x F_p` | **EXACT** |
| `p=3 mod4` Gaussian ring is `F_{p^2}` | **EXACT** |
| `p=2`: `(1+i)^2=0` | **EXACT** |
| bounded multi-modular data -> ordinary integer via CRT | **EXACT** |
| one finite modulus -> Archimedean ordered reals | **NO-GO** |
| q=2 microscopic orientation dynamically selects the phase orientation | **OPEN** |
| finite q=2 phase groups refine to physical U(1) | **OPEN** |
| modular winding/history derives physical real frequency/action | **OPEN** |
| physical constants derived from this arithmetic layer | **NOT CLAIMED** |

---

# 8. Executable gate

Run

```bash
python scripts/modular_complex_arithmetic_bridge_gate.py
```

or

```bash
python scripts/modular_complex_arithmetic_bridge_gate.py \
  --json verification_results/modular_complex_arithmetic_bridge_gate.json
```

The gate uses exact integer arithmetic only and checks:

- exhaustive matrix homomorphism tests for moduli `2,3,5,7,8,13`;
- exact split and inverse maps for primes `5,13,17`;
- nonsplit quadratic-field behavior for `3,7,11,19`;
- ramification at `p=2`;
- the `C4 -> J` block;
- the `mu4 -> mod5` phase table;
- componentwise CRT recovery of ordinary Gaussian integer sum/product.

---

# 9. Strongest current interpretation

The result does **not** prove that nature fundamentally computes modulo a prime.

What it does prove is more useful and more conservative:

> Complex phase arithmetic and ordinary arithmetic are not disjoint mathematical worlds. Finite cyclic dynamics can be represented by complex phases, by modular integers, or by ordinary real/integer matrices, with exact equivalences in declared sectors. Ordinary unbounded ordered arithmetic requires additional lift/scale/history information.

For this project the new candidate route is therefore

\[
\boxed{
\text{binary relations}
\to\text{oriented cyclic dynamics}
\to\text{finite phase arithmetic}
\to\text{complex spectral language}
\leftrightarrow\text{ordinary modular arithmetic}
\to\text{history/CRT/scale lift}
\to\text{ordinary physical observables}.
}
\]

The next nontrivial task is to determine whether the q=2 microscopic dynamics actually selects the oriented cyclic/phase lift and whether that lift composes consistently under recursive refinement.
