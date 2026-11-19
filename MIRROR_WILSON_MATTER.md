# Mirror-covariant Wilson-Dirac matter carrier

Status: **finite concrete H_m(sigma) carrier; realistic PL/gauge/chiral matter and numerical beta_m remain open**.

The matter matrix-element theorem defines

```text
beta_m=(chi m)^-1 <dH_m/dsigma>_rest.
```

This note supplies the first explicit lattice Hamiltonian on which that quantity can be evaluated while retaining a standard fermion-doubler control.

---

## 1. Hamiltonian

On a three-dimensional spatial lattice use

```text
H_q(k,sigma)
 = sum_i alpha_i sin(k_i)
 + beta_D M_q(k,sigma),
```

with

```text
M_q(k,sigma)
 = m_q(sigma)
 + r_W sum_i [1-cos(k_i)].
```

For the constant-beta positive-mass control,

```text
m_q(sigma)=m_* exp(q beta sigma).
```

The Wilson coefficient satisfies `r_W>0`.

---

## 2. Mirror covariance

Mirror acts on the matter doublet as

```text
(sigma,q)->(-sigma,-q).
```

The mass law obeys

```text
m_q(sigma)=m_-q(-sigma).
```

The Wilson term is independent of `q` and `sigma`, so

```text
M_q(k,sigma)=M_-q(k,-sigma).
```

Hence the complete one-particle spectra are identical:

```text
Spec H_q(k,sigma)
 = Spec H_-q(k,-sigma).
```

For the aligned mirror vacua `(sigma,q)=(chi v,chi)`, both branches therefore have the same positive rest mass and the same positive-energy dispersion.

---

## 3. Static mirror source

Only the physical mass term depends on `sigma`, so

```text
dH/dsigma
 = beta_D dm_q/dsigma.
```

For the exponential control,

```text
dm_q/dsigma=q beta m_q.
```

At `k=0` and positive mass,

```text
<beta_D>=1,
```

therefore

```text
<dH/dsigma>_rest
 = q beta m_q.
```

For aligned mirror states `q=chi`,

```text
Q_sigma=chi beta m_phys,
```

so

```text
beta_m=beta.
```

At finite momentum the Wilson effective mass is

```text
M=m_q+r_W sum_i(1-cos k_i)
```

and the Hellmann--Feynman control becomes

```text
dE/dsigma
 = (dm_q/dsigma) M/E.
```

The executable gate verifies this matrix element directly.

---

## 4. Wilson doubler control

Set the physical mass to zero only for the corner test.

At the eight three-dimensional Brillouin corners `k_i=0 or pi`, the naive kinetic sine terms vanish. Without the Wilson term all eight corners would be zero modes.

With `r_W=1`,

```text
M_W(k)=sum_i(1-cos k_i)
```

is zero only at

```text
k=(0,0,0).
```

The other seven corners acquire nonzero Wilson mass.

`scripts/mirror_wilson_matter_gate.py` verifies exactly one massless corner in this negative/continuum control.

---

## 5. Deterministic finite control

Using

```text
m_*=0.4
beta=0.37
v=0.8
r_W=1,
```

the two aligned mirror branches have

```text
m_phys=0.5377880627328211
```

and opposite rest sources

```text
Q_sigma=+0.1989815832111438
Q_sigma=-0.1989815832111438.
```

Their tested momentum spectra agree to machine precision.

Thus a nonzero static mirror matrix element is compatible, at finite free-fermion level, with

```text
positive energy
mirror spectral degeneracy
opposite sigma charge
Wilson corner-doubler removal.
```

---

## 6. What this does NOT solve

This is a carrier, not the final matter theory.

Still open:

- deriving `beta` rather than entering it;
- placing fermions on the actual irregular recursive PL/Peter-Weyl geometry;
- local gauge representations and realistic gauge group;
- a genuinely chiral low-energy spectrum;
- the Wilson finite-spacing chiral-symmetry tradeoff;
- generations, Higgs/Yukawa structure and anomaly cancellation;
- derivation of the alignment between matter mirror label `q` and background `chi`;
- full gravity x route x mirror-matter HDA.

The result nevertheless closes an important logical question: the required rest-energy derivative does not force negative energy or reintroduce naive corner doublers.

---

## Reproduction

```bash
python scripts/mirror_wilson_matter_gate.py \
  --output verification_results/MIRROR_WILSON_MATTER.json
```
