# Physical boundary source dressing through the enlarged master projector

Status: **exact finite operator protocol + deterministic positive control target; no BQG physical zero mode or cosmological observable is claimed here.**

## 1. Why this layer is now mandatory

The certified full `32 x 32` Euclidean logical Gram/master baseline is positive definite on the bare all-`j=1/2` logical carrier. Therefore the physical zero sector cannot be searched for by diagonalizing only a master operator compressed to that carrier.

The correct object is the zero projector of the **enlarged parity-complete Peter-Weyl / graph-changing habitat** and its overlap with the q=2 boundary carrier.

Let

\[
B=(|b_1\rangle,\ldots,|b_r\rangle)
\]

be the frozen boundary block and let the full finite-regulator master constraint be

\[
\mathbb M\ge0,
\qquad
P_0=\mathbf 1_{\{0\}}(\mathbb M).
\]

The first physical boundary observable is

\[
\boxed{G_0=B^\dagger P_0 B.}
\]

A positive-definite compressed master

\[
B^\dagger\mathbb M B>0
\]

does **not** imply `G_0=0`.

## 2. Heat-kernel approximation

Define

\[
K_\tau=e^{-\tau\mathbb M}.
\]

If the finite master gap is `Delta_M>0`, then

\[
\|K_\tau-P_0\|=e^{-\tau\Delta_M}.
\]

For an isometric boundary block, `B^dagger B=I`,

\[
\boxed{
G_\tau=B^\dagger K_\tau B
\longrightarrow
G_0
}
\]

with

\[
\|G_\tau-G_0\|\le e^{-\tau\Delta_M}.
\]

The heat parameter `tau` is a **projection depth / rigging regulator**. It is not physical time and must never be renamed `omega^{-1}` or a propagation time.

## 3. Project an observable before building the source functional

Let `O_a` be a microscopic/coarse operator such as a shape, volume or relational metric insertion.

The exact physical operator seen by the boundary-projected zero sector is represented by

\[
A_{a,0}=B^\dagger P_0 O_a P_0 B.
\]

When `G_0` has support rank `r_0`, whiten only on that support:

\[
\boxed{
\bar O_a
=G_0^{-1/2}A_{a,0}G_0^{-1/2}.
}
\]

This is equivalent to using the orthonormal physical boundary images

\[
Q_0=P_0 B G_0^{-1/2},
\qquad Q_0^\dagger Q_0=I_{r_0},
\]

because

\[
\bar O_a=Q_0^\dagger O_a Q_0.
\]

Thus the source is attached **after physical projection**, not to the constraint resolvent.

## 4. Heat-kernel estimator computable by block Krylov

The matrices required above have direct finite-`tau` estimators:

\[
\boxed{
G_\tau=B^\dagger K_\tau B,
}
\]

\[
\boxed{
A_{a,\tau}
=B^\dagger K_{\tau/2}O_aK_{\tau/2}B,
}
\]

and

\[
\boxed{
\bar O_{a,\tau}
=G_\tau^{-1/2}A_{a,\tau}G_\tau^{-1/2}
}
\]

on the numerically stable support of `G_tau`.

As `tau -> infinity`,

\[
A_{a,\tau}\to B^\dagger P_0O_aP_0B,
\]

so the whitened operator converges to the exact physical boundary operator.

Crucially, `G_tau` and `A_a,tau` can be evaluated from a block-Lanczos/Krylov representation started by `B`; a dense diagonalization of the ambient Peter-Weyl habitat is not required.

## 5. Physical source functional on the projected support

After the physical boundary density matrix `rho_phys` has been frozen, define

\[
\boxed{
Z[J]
=\operatorname{Tr}_{\rm supp G_0}
\left[
\rho_{phys}
\exp\left(\sum_a J_a\bar O_a\right)
\right].
}
\]

Then

\[
W[J]=\log Z[J]
\]

generates connected physical boundary cumulants.

For a maximally mixed finite positive control,

\[
\rho_{phys}=I/r_0.
\]

The production BQG calculation must instead freeze the boundary/semi-classical state independently; the maximally mixed choice is not a cosmological vacuum prediction.

The important ordering is

\[
\boxed{
\mathbb M
\to P_0
\to G_0,A_{a,0}
\to \bar O_a
\to Z[J]
\to W[J]
\to \Gamma.
}
\]

## 6. Why powers of the projected operator matter

One must distinguish

\[
P_0e^{JO}P_0
\]

from

\[
e^{J(P_0OP_0)}
\]

when `[P_0,O] != 0`.

The second expression is the correct exponential of the operator acting **inside the physical Hilbert space**. Its quadratic term contains

\[
P_0OP_0OP_0,
\]

not merely `P_0 O^2 P_0`.

This distinction is essential for connected two-point functions and is one reason the source functional is constructed only after the physical projection/whitening step.

## 7. Deterministic four-dimensional positive control

The new gate uses a full Hilbert space with orthonormal states

\[
p_1=(e_1+e_3)/\sqrt2,
\quad
p_2=(e_2+e_4)/\sqrt2,
\]

\[
q_1=(e_1-e_3)/\sqrt2,
\quad
q_2=(e_2-e_4)/\sqrt2.
\]

Set

\[
\mathbb M=2|q_1\rangle\langle q_1|+5|q_2\rangle\langle q_2|,
\]

so

\[
P_0=|p_1\rangle\langle p_1|+|p_2\rangle\langle p_2|.
\]

Use the boundary block

\[
B=(e_1,e_2).
\]

Then

\[
\boxed{
B^\dagger\mathbb M B
=\operatorname{diag}(1,5/2)
}
\]

is strictly positive and has no zero vector, while

\[
\boxed{
G_0=B^\dagger P_0B=\frac12I_2\ne0.
}
\]

This is the minimal exact analogue of the logical-carrier obstruction discovered in the BQG Peter-Weyl calculation.

Define physical Pauli operators

\[
O_X=|p_1\rangle\langle p_2|+|p_2\rangle\langle p_1|,
\]

\[
O_Z=|p_1\rangle\langle p_1|-|p_2\rangle\langle p_2|.
\]

After whitening,

\[
\bar O_X=X,
\qquad
\bar O_Z=Z.
\]

For the maximally mixed projected boundary state,

\[
\boxed{
Z(j_X,j_Z)=\cosh\sqrt{j_X^2+j_Z^2},
}
\]

and

\[
\boxed{
\partial_a\partial_b W|_0=\delta_{ab}.
}
\]

Thus the same source algebra as the earlier q=2 relational positive control can survive even when the physical zero modes live outside the original boundary carrier.

## 8. Production BQG target

The next genuine calculation replaces this four-dimensional control by

```text
B = 32 q=2 logical boundary columns
O = collective volume + transported metric/shape insertions
M = full regulated Euclidean + Lorentzian master on the parity-complete
    constraint-generated Peter-Weyl habitat
```

and computes, without fitting cosmological functions,

\[
G_\tau,
\quad
A_{V,\tau},
\quad
A_{shape,\tau},
\quad
\bar O_{a,\tau}
\]

along a frozen `tau` / Krylov-depth / cutoff staircase.

Only if these quantities stabilize and `G_tau` retains nonzero support may the programme proceed to connected two-block source response and then to

\[
\Gamma_{scalar}^{(2)}(\omega,k).
\]

## 9. Falsifiers

The physicalization branch must stop or change architecture if, under regulator/refinement control,

1. `B^dag P_0 B -> 0` for every declared q=2 boundary family;
2. the apparent support is created only by an arbitrary spectral window and disappears as the window closes;
3. source matrices fail to converge under the same projector sequence;
4. different positive master metrics `G^{AB}` change the limiting zero-sector response rather than only convergence/conditioning;
5. a claimed cosmological pole exists only before physical projection or only before constraint/gauge reduction.

This protocol closes no physical cosmology gate by itself. It freezes the operator order required to ask the question correctly.
