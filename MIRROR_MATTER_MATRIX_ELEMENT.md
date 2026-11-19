# Mirror matter matrix element: static-source gate

Status: **exact mirror-covariant mass-derivative theorem + exact free-Dirac source-selection result + finite positive-mass control; numerical microscopic coupling still open**.

The mirror branch has reached

```text
Y_L/Q -> Sigma -> sigma(x) -> Z_sigma -> beta_m -> alpha -> F_mirror.
```

The remaining matter coefficient is no longer assigned by hand: it is a rest-energy matrix element.

---

## 1. The old axial bridge is not a universal static mass source

The mirror/chirality construction naturally allows

```text
sigma * J5^0,
```

because both factors are parity odd.

That remains a legitimate chirality/spin-sensitive coupling. But for a free on-shell massive Dirac state,

```text
J5^0/J^0 = h |p|/E,
```

with helicity `h=+/-1`.

Therefore

```text
p=0 -> J5^0=0,
```

and for an unpolarized ensemble

```text
<J5^0> = 0
```

at every momentum. The diagonal pseudoscalar bilinear also vanishes:

```text
ubar i gamma5 u = 0.
```

`scripts/mirror_matter_matrix_element_gate.py` verifies these statements with explicit Dirac matrices and positive-energy spinors.

Thus the old axial bridge survives for chirality/spin/velocity-dependent physics, but not as a universal static charge of cold unpolarized massive matter.

---

## 2. The static source is an energy derivative

For a normalized stationary matter state, Hellmann--Feynman gives

```text
dE/dsigma = <dH_m/dsigma>.
```

For a rest state with `m=E_rest`, define

```text
beta_m
 = (1/(chi*m))
   <dH_m/dsigma>_rest.
```

Equivalently,

```text
Q_sigma = chi beta_m m.
```

This is the operational definition of the missing matter coefficient.

If the derivative vanishes for every admissible physical state,

```text
beta_m=0
```

and the static mirror-force branch fails regardless of how robust the geometric mirror order is.

---

## 3. General mirror-mass derivative theorem

Introduce a matter mirror label

```text
q=+/-1
```

with mirror transformation

```text
(sigma,q) -> (-sigma,-q).
```

Let the positive rest-mass function satisfy only mirror covariance:

```text
m_q(sigma)=m_-q(-sigma).
```

Take the two aligned mirror vacua

```text
(q,sigma)=(+1,+v)
```

and

```text
(q,sigma)=(-1,-v).
```

Then mirror covariance immediately gives

```text
m_+(v)=m_-(-v).
```

So exact mirror partners have equal rest mass.

Differentiate the covariance identity:

```text
m_+'(sigma)=-m_-'(-sigma).
```

At the aligned vacua,

```text
m_+'(v)=-m_-'(-v).
```

Thus, whenever the derivative is nonzero, mirror symmetry itself gives **equal positive mass and opposite static sigma charge**.

The charge-to-mass magnitude is

```text
beta_m
 = m_+'(v)/m_+(v)
 = -m_-'(-v)/m_-(-v).
```

This conclusion does **not** require the exponential mass law.

It is the general microscopic quantity to compute.

---

## 4. Constant-beta special case

If one further requires a constant logarithmic response,

```text
d ln m_q/dsigma = q beta,
```

the solution at fixed normalization is

```text
m_q(sigma)=m_* exp(q beta sigma).
```

It is positive for every real `sigma`.

For aligned mirror states `(sigma,q)=(chi v,chi)`,

```text
m_phys=m_* exp(beta v)
```

is identical in the two branches, while

```text
dm/dsigma=chi beta m_phys.
```

Hence

```text
beta_m=beta.
```

The exponential model is therefore a globally positive constant-`beta_m` control, not an assumption required by the general sign theorem.

---

## 5. Dirac Hellmann--Feynman control

For

```text
H_D = alpha.p + beta_D m_q(sigma),
```

we have

```text
dH_D/dsigma=beta_D dm_q/dsigma.
```

For a positive-energy Dirac state,

```text
<beta_D>=m/E,
```

so

```text
dE/dsigma=(dm/dsigma) m/E.
```

In the constant-beta control this becomes

```text
dE/dsigma=q beta m^2/E.
```

At rest,

```text
dE_rest/dsigma=q beta m.
```

The executable gate verifies these identities to machine precision.

---

## 6. Direct extraction from a microscopic spectrum

The completely general local extraction is

```text
beta_m
 = (1/m)
   d m_+(sigma)/d sigma |_(sigma=v).
```

For the constant-beta exponential control one may equivalently use a finite mirror ratio:

```text
beta
 = [ln m_+(sigma)-ln m_-(sigma)]/(2 sigma).
```

Thus the next actual matter calculation is

```text
construct H_m(sigma)
 -> diagonalize mirror-resolved rest spectrum
 -> differentiate the physical mass eigenvalue
 -> obtain beta_m
 -> compute alpha.
```

No fifth-force strength needs to be fitted.

---

## 7. Force strength after the matter theorem

The continuum normalization gives

```text
alpha=beta_m^2/(4*pi*G*Z_sigma).
```

Using the regular-seed Hodge stiffness

```text
Z_sigma=(2*sqrt(2)/3)J/ell,
```

one gets

```text
alpha
 = 3 beta_m^2 ell
   /(8*sqrt(2)*pi*G*J).
```

The repulsion criterion is

```text
alpha > exp(m_sigma r)/(1+m_sigma r).
```

Equivalently, the microscopic mass derivative must satisfy

```text
|d ln m/dsigma|_vac
 > sqrt(
     4*pi*G*Z_sigma
     * exp(m_sigma r)/(1+m_sigma r)
   ).
```

This is the sharp matter-side threshold.

---

## 8. Remaining selection problem

The theorem says what happens **if** low-energy matter has a mirror-covariant `sigma`-dependent rest spectrum. The repository has not yet derived:

- the realistic microscopic matter Hilbert space;
- its gauge/chiral representations;
- the actual `sigma`-dependent mass matrix;
- why the low-energy branch aligns its matter mirror label with the background mirror order;
- the numerical derivative `beta_m`;
- the fate of the mirror partner and all local/global anomalies.

So the form of the source is now fixed, but its physical magnitude remains an open matter calculation.

---

## Reproduction

```bash
python scripts/mirror_matter_matrix_element_gate.py \
  --output verification_results/MIRROR_MATTER_MATRIX_ELEMENT.json
```
