# Finite phase -> ordinary integer/rational arithmetic -> one physical scale

**Status:** exact arithmetic bridge plus an open physical scale/refinement interpretation.

This note continues `MODULAR_COMPLEX_TO_ORDINARY_ARITHMETIC.md` one step further.

The previous result showed that finite complex/modular arithmetic can be rewritten exactly as ordinary matrix arithmetic, scalar residue channels, or cyclic integer addition. The remaining question is more physical:

> How can a wrapped modular phase become the ordinary unbounded arithmetic used for action, frequency, length and other measured quantities?

The answer naturally separates into three layers:

```text
wrapped phase
-> integer residue + winding/history
-> ordinary integer/rational arithmetic
-> one common dimensional scale
-> real-valued continuum observable.
```

The first two arrows are exact. The last scale/continuum arrow is not fixed by modular arithmetic alone.

---

# 1. A phase is already ordinary modular addition in disguise

For

\[
\omega_M=e^{2\pi i/M},
\]

we have

\[
\omega_M^a\omega_M^b
=\omega_M^{a+b\pmod M}.
\]

Therefore

\[
\boxed{\mu_M\cong\mathbb Z_M}.
\]

A multiplicative complex phase can be stored as one ordinary residue `k mod M`.

The complex representation is useful because cyclic translations diagonalize in that basis. It is not the only exact arithmetic representation.

---

# 2. The missing information is winding

Knowing

\[
k\pmod M
\]

does not tell us the unique integer `k`.

All integers

\[
k+wM,\qquad w\in\mathbb Z
\]

have the same residue.

The extra integer `w` is a **winding/carry/history register**.

For ordinary nonnegative representatives `a,b in {0,...,M-1}`,

\[
a+b=r+wM,
\]

where

\[
r=(a+b)\bmod M
\]

and

\[
w=\left\lfloor\frac{a+b}{M}\right\rfloor.
\]

Thus

\[
\boxed{
\text{modular residue} + \text{winding}
= \text{ordinary integer}.
}
\]

For a whole history,

\[
K_{\rm total}=r_{\rm final}+M W_{\rm total}.
\]

The finite phase remembers only `r_final`; the history remembers how many times the arithmetic crossed the branch cut.

This suggests a potentially important physical interpretation:

```text
compact phase = local observable state
winding register = information carried by history
unwrapped integer = additive history count.
```

That interpretation is a candidate, not yet a derived law of the microscopic Hamiltonian.

---

# 3. Universal-cover version

The continuum analogue is the circle group

\[
U(1)=\mathbb R/(2\pi\mathbb Z).
\]

The exponential map

\[
\theta\mapsto e^{i\theta}
\]

forgets integer multiples of `2pi`.

Therefore a compact phase determines only

\[
\theta\pmod{2\pi}.
\]

An ordinary additive phase requires a lift to the universal cover plus winding:

\[
\boxed{
\text{phase on }U(1)
+\text{winding/history}
\to\theta\in\mathbb R.
}
\]

This is structurally similar to the repository's broader lesson that a constraint spectral label is not automatically a physical time variable: history carries information that a local compact/spectral object alone does not contain.

---

# 4. Why one common scale survives

Suppose the arithmetic/history layer has produced an ordinary integer count

\[
N\in\mathbb Z.
\]

Now ask for a dimensionful physical quantity `Q` that composes additively with history:

\[
Q(N_1+N_2)=Q(N_1)+Q(N_2).
\]

Because the domain is the integers, this immediately forces

\[
\boxed{Q(N)=sN}
\]

with

\[
s=Q(1).
\]

No continuity assumption is needed.

This gives a clean mathematical interpretation of the **one common scale** that repeatedly survives the information-graph physicalization programme:

```text
arithmetic fixes the integer composition law
history fixes the winding/lift
one scalar s converts a dimensionless count into physical units.
```

The arithmetic cannot determine whether `s` is one metre, one Planck length, one unit of action, or something else.

That is a scale-setting problem, not an arithmetic problem.

This is exactly why several unrelated physical observables must not be calibrated independently after the fact.

---

# 5. Modular arithmetic can reconstruct exact rational arithmetic

Integers are not enough for geometry; ratios matter.

Suppose a reduced rational number

\[
q=\frac ab
\]

is represented modulo `M` by

\[
x\equiv a b^{-1}\pmod M,
\]

with `gcd(b,M)=1`.

Assume bounds

\[
|a|\le A,
\qquad
0<b\le B,
\]

and

\[
\boxed{2AB<M}.
\]

Then there can be at most one reduced fraction inside those bounds with the same residue.

### Proof of uniqueness

If

\[
a_1b_1^{-1}\equiv a_2b_2^{-1}\pmod M,
\]

then

\[
M\mid(a_1b_2-a_2b_1).
\]

But

\[
|a_1b_2-a_2b_1|\le2AB<M.
\]

The only multiple of `M` with magnitude smaller than `M` is zero. Hence

\[
a_1b_2=a_2b_1,
\]

so the reduced fractions coincide.

Therefore

\[
\boxed{
\text{sufficiently large modular arithmetic}
+\text{bounds}
\to\text{exact ordinary rational arithmetic}.
}
\]

The executable gate reconstructs uniquely, with

```text
M = 5*7*11*13 = 5005
A = 20
B = 30
```

the ordinary fractions

\[
17/29
\]

and

\[
-19/23.
\]

---

# 6. Why compatible residues alone are still not the ordinary integers

It is tempting to say:

```text
all modular residues together = Z.
```

That is too strong.

A compatible family of residues over all moduli naturally lives in the **profinite completion** of the integers. The ordinary integers embed into that object, but the profinite completion contains additional elements.

Therefore a physical theory still needs a selection principle such as:

- a finite bound;
- a growth condition;
- a history/winding construction;
- an Archimedean scale/continuum condition.

This is another no-go against treating a formal inverse limit as automatically identical to the measured number line.

---

# 7. General finite Fourier phases can be modular integers

The q=2 case starts with the four-cycle `C4`, but the principle is general.

Let `S_M` be the integer cyclic shift on `M` states.

Over the complex numbers it is diagonalized by the discrete Fourier basis with eigenvalues

\[
1,\omega_M,\omega_M^2,\ldots,\omega_M^{M-1}.
\]

Now choose a prime `p` for which

\[
\boxed{M\mid(p-1)}.
\]

Then `F_p^*` contains an element `r` of exact order `M`.

The same integer shift is diagonalized over `F_p` by a modular Fourier/Vandermonde basis, with eigenvalues

\[
1,r,r^2,\ldots,r^{M-1}.
\]

Thus

\[
\boxed{
\text{complex DFT phase algebra}
\leftrightarrow
\text{exact modular-integer Fourier phase algebra}
}
\]

for the same cyclic operator.

The executable gate checks exact modular diagonalization for

```text
M=4,  p=5
M=8,  p=17
M=12, p=13
M=16, p=17.
```

It also exhaustively checks small prime/order pairs and finds an exact-order `M` element exactly when `M` divides `p-1` in the tested range.

This is the algebra behind the number-theoretic transform: complex roots of unity are not required to perform the cyclic spectral arithmetic.

---

# 8. A possible finite precursor of U(1)

One discrete four-cycle gives only

\[
\mu_4=\{1,i,-1,-i\}.
\]

That is not the full circle `U(1)`.

But consider a nested dyadic family

\[
\mu_4\subset\mu_8\subset\mu_{16}\subset\cdots.
\]

The embedding is

\[
\omega_{M}^k\mapsto\omega_{2M}^{2k}.
\]

The union of dyadic roots is dense in `U(1)`. For any target angle `theta`, choosing the nearest dyadic phase at order `M=2^g` gives

\[
\boxed{
|\theta-2\pi k/M|\le\pi/M.
}
\]

So there is a mathematically clean continuum target:

\[
\boxed{
\mu_{2^g}\xrightarrow[g\to\infty]{\rm dense}U(1).
}
\]

### Important physical boundary

The information-graph model currently has an exact q=2 `C4` structure, but it has **not yet derived** the law

\[
M_g=2^{g+2}
\]

for a physical phase refinement.

Therefore

```text
C4 -> mu4                         exact representation result
mu4 -> mu8 -> mu16 -> ... -> U1  mathematically valid candidate refinement
microscopic dynamics selects it  OPEN.
```

---

# 9. Why this is interesting for q=2

We now have two finite structures that should not be conflated but can interact:

1. the microscopic route-label group

\[
\mathbb Z_2^2,
\]

2. the oriented automorphism cycle

\[
C_4.
\]

The first gives the Walsh/tetrahedral geometry bridge.

The second gives the order-four shift whose real two-dimensional block satisfies

\[
J^2=-I
\]

and whose spectral phases are

\[
\pm i.
\]

This suggests a new division of labour:

```text
Z2^2 label algebra -> tetrahedral spatial geometry
oriented C4 transport -> finite phase/complex structure.
```

That is much cleaner than forcing one algebra to do both jobs.

---

# 10. Candidate connection to the geometry-qubit orientation Y

The existing logical geometry qubit contains an orientation pseudoscalar channel `Y`.

The new cyclic theorem says

```text
cycle orientation reversal
S -> S^-1
J -> -J
i -> -i.
```

Therefore there is a precise open question:

> Does the microscopic sign/orientation sector that selects the logical `Y` orientation also select the sign of the finite complex structure `J`?

If the answer is yes, the theory would gain a derived relation

\[
\boxed{
\text{geometric orientation reversal}
\leftrightarrow
\text{complex conjugation}.
}
\]

This would be a real structural result, not numerology.

It remains **OPEN** until one operator acts on both declared sectors and the covariance is checked.

---

# 11. Candidate connection to action and physical phase

Quantum amplitudes have the form

\[
e^{iS/\hbar}.
\]

A finite phase version could encode only an action index modulo `M`:

\[
A_h=\omega_M^{n_h}.
\]

History composition gives

\[
A_{h_1}A_{h_2}=\omega_M^{n_1+n_2\pmod M}.
\]

A winding/history lift supplies an integer

\[
N_h=n_h+MW_h.
\]

Then the most general additive physical action assignment on integer counts is

\[
\boxed{S_h=s_A N_h}.
\]

This does **not** derive `hbar` or the numerical action unit. But it gives a concrete interpretation of why a single conversion scale survives after the exact finite composition law is fixed.

Potential future target:

```text
microscopic oriented path
-> modular phase index
-> history winding
-> integer action count
-> one common action scale
-> continuum phase.
```

---

# 12. New arithmetic truth table

| Arrow | Status |
|---|---|
| finite root phase multiplication -> `Z_M` addition | **EXACT** |
| residue + winding -> ordinary integer | **EXACT** |
| bounded modular residue -> rational reconstruction | **EXACT** |
| `M | p-1` -> modular Fourier root in `F_p` | **EXACT theorem; finite gate checks cases** |
| integer cyclic shift -> complex or modular Fourier spectrum | **EXACT** |
| nested dyadic root groups dense in U(1) | **EXACT mathematical limit** |
| q=2 dynamics selects dyadic phase refinement | **OPEN** |
| logical orientation `Y` selects sign of complex structure `J` | **OPEN** |
| winding register is physically the relational history variable | **OPEN** |
| integer count -> physical dimensionful quantity `sN` | **EXACT if additivity is assumed** |
| numerical value of the common scale `s` | **OPEN PHYSICAL** |
| modular arithmetic alone -> ordered real numbers | **NO-GO** |

---

# 13. Strongest new conceptual chain

The combined result of the two arithmetic notes is now

\[
\boxed{
\begin{aligned}
&\text{binary q=2 relations}\\
&\to \mathbb Z_2^2\text{ geometry labels}\\
&\to C_4\text{ oriented transport}\\
&\to J^2=-I\\
&\to \mu_4\text{ finite complex phase}\\
&\leftrightarrow \mathbb Z_4\text{ modular addition}\\
&\leftrightarrow \text{ordinary modular matrices/scalars}\\
&\to \text{winding/CRT lift}\\
&\to \mathbb Z,\mathbb Q\\
&\to \text{one common physical scale}\\
&\to \mathbb R/\mathbb C\text{ continuum observables in a controlled limit}.
\end{aligned}}
\]

The first several arrows are now exact arithmetic/representation results. The physical selection of the refinement, winding/history interpretation and scale remains the next frontier.

---

# Reproduction

```bash
python scripts/finite_phase_integer_rational_lift_gate.py
```

or

```bash
python scripts/finite_phase_integer_rational_lift_gate.py \
  --json verification_results/finite_phase_integer_rational_lift_gate.json
```
