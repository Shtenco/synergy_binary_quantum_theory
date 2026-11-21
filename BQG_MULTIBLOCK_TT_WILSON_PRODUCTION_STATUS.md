# BQG five-block bridge — production status

Status: **the exact five-block 30-column microscopic P carrier is now produced from the frozen L1 amplitudes and passes its full-global Gram/rank guards. The target-independent E/S image contraction and the connected physical 1PI bridge remain open.**

No new microscopic dynamics is introduced.

## Frozen operator

\[
\boxed{G_{frozen}=-\frac23H_E^{sine}-\frac{32}{9}S},
\qquad
S=-\frac{i}{2}(L_{raw}-L_{raw}^\dagger).
\]

`R_op` is not part of this gravitational TT spatial precursor.

## Exact retained carrier now produced

The frozen six q directions of one block are the already-derived background-orthogonal strict-interior q4 coarse-edge microscopic tangents, ordered

```text
(01),(02),(03),(12),(13),(23).
```

They are the same local q frame used by the measured map `M_hq`. The static maximal-symmetric j=3 block is rank one and is not a production metric carrier.

The production script

```text
scripts/bqg_fiveblock_p_carrier_producer.py
```

uses the frozen upstream commit

```text
2b2c9f623544f5d38f7ffd7f37617f91f4dae306
```

and constructs the actual centered patch

```text
parent 0 + face neighbours 1,2,4,8
```

inside the full L1 384-node / 768-link Gauss basis. The five coarse parent tetrahedra are

```text
0 : [0,2,4,6]
1 : [0,2,4,7]
2 : [0,2,5,6]
4 : [0,3,4,6]
8 : [1,2,4,6]
```

Each parent contributes 24 exact strict q4 source columns, combined by the already frozen S4 parity/coset rule into six normalized coarse-edge columns. Thus

\[
P\in\mathbb C^{960\times 30}.
\]

Every normalized P column has microscopic support 32. The measured production result is

\[
\operatorname{rank}(P)=30,
\]

\[
\frac{\|P^\dagger P-I_{30}\|}{\|P^\dagger P\|}
=4.5324665183683945\times10^{-17},
\]

with Hermiticity defect 0, eigenvalues in `[0.9999999999999999, 1.0]`, cross-block maximum absolute overlap exactly 0, and cross-block global-basis collision count exactly 0.

The common unnormalized six-edge norm square in every parent is

\[
\nu=0.6982713481539777,
\]

and the relative common-norm-I6 defect is

\[
1.2981963330876717\times10^{-16}.
\]

The serialized NPZ reproduces the direct sparse Gram with relative defect

\[
1.7786683001560746\times10^{-16}.
\]

### Why shared-face recoupling is not part of P

For this strict-interior carrier every changed microscopic spin belongs to an edge internal to exactly one parent, all parent-boundary spins remain at the background value, and all exterior labels remain exactly background. Therefore P columns from distinct parents have disjoint full-global Gauss-basis support and are exactly orthogonal.

Consequently complete six-link shared-face recoupling must **not** be inserted artificially before P. It becomes necessary only after `H_E^sine` or `S` acts and produces boundary-touching/crossing histories.

The independent frozen complete face-recoupling prerequisite was re-run in the same production CI: all 61 q4 face patterns pass, dimensions range 16..144, total doubled-spin support is `[0,2,4,6,8]`, worst unitarity defect is `5.1388140815939885e-15`, and the historical maximal-J3 Dicke subblock defect is `1.1102230246251565e-16`.

## Production basis without dense Hilbert materialization

The downstream assembler continues to accept a general raw P Gram for fail-closed reuse, but for the actual strict five-block carrier the measured result is numerically `K_P=I30` to machine precision.

The next heavy producer must:

1. load the actual 30 P columns;
2. apply the already frozen `H_E^sine` and Hermitian `S` to them without target pruning;
3. retain the resulting boundary-touching histories in the complete shared-face SU(2) representation;
4. project out `span(P)` and whiten only the residual amplitude-level Q span;
5. form

\[
W=[P,Q_Q]
\]

for this measured orthonormal P carrier;
6. measure directly

\[
E_W=W^\dagger H_E^{sine}W,
\qquad
S_W=W^\dagger S W,
\]

plus compression leakage

\[
\eta_O=\|(1-WW^\dagger)OW\|/\|OW\|,
\qquad O\in\{H_E^{sine},S\}.
\]

The production seam requires the declared basis to be closed at the frozen tolerance and refuses target fitting.

## Signed compressed operator

`scripts/bqg_signed_5block_operator_assembler.py` consumes `E_basis`, `S_basis`, raw `P_gram`, labels/geometry and provenance and forms only

\[
G_W=-\frac23E_W-\frac{32}{9}S_W.
\]

It never asks for a full microscopic `N x N` Hilbert matrix.

## Schur in orthonormal P/Q space

With

\[
G_W=\begin{pmatrix}A&B\\B^\dagger&D\end{pmatrix}
\]

the extractor uses

\[
C_{eff}^{orth}=A-BD^+B^\dagger.
\]

A zero/gapless Q mode coupled to P causes the hard stop

```text
GAPLESS_COUPLED_Q_MODE_REQUIRES_PROMOTION
```

rather than an `i eta` or fitted denominator.

The effective quadratic form is returned to the original labelled q frame:

\[
C_{eff}^{q}=K_P^{1/2}C_{eff}^{orth}K_P^{1/2},
\]

and the normalized-state Hessian is

\[
\boxed{K_q=2\Re C_{eff}^{q}-2C_{00}\Re K_P}.
\]

Then

\[
K_h=M_{hq}^{-T}K_qM_{hq}^{-1}.
\]

The spatial center+four-neighbour Fourier symbol is TT-projected and tested for masslessness, positive/isotropic leading k2, reciprocity and exact closure in the frozen six-dimensional quartic basis.

A passing run may emit only

\[
\boxed{\mathbf c_{micro,spatial}^{BQG}}.
\]

`c_BQG_IR` remains null.

## Existing reusable upstream results

The recovered `research/bcqg-core-candidate-v1` branch already supplies:

- L1 closed 16-cell habitat: 384 fine tetrahedra / 768 dual links;
- 16 canonical parent blocks, four face-neighbours each;
- exact physical-sine Euclidean action, including generic `H_sine_state(state,v,Jmax2)` on a full Gauss state;
- exact complete six-link SU(2) face recoupling;
- exact strict q4 boundary contraction and six-edge tangent rank;
- exact background removal and canonical orthonormal six-edge frame;
- q4 crossing/environment sectors restoring full source rank 24;
- one-E face support through j<=4;
- conservative one-S face support through j<=6;
- depth-2 target-independent wall through j<=9;
- Hermitian Lorentzian primitives;
- measured full-rank `M_hq`, with condition number sqrt(2);
- sparse Gram/whitening and Schur-gap algorithms.

These are reused as operator primitives, not redefined.

## Remaining heavy amplitude calculation

```text
load actual 960 x 30 strict P carrier
 -> apply frozen H_E^sine and Hermitian S to P without target pruning
 -> complete shared-face recoupling only for the boundary-touching E/S images
 -> target-independent residual Q histories through the frozen closure wall
 -> Q residual whitening
 -> direct E_W and S_W matrix elements
 -> leakage / support-wall closure
 -> signed G_W
 -> Q-gap / Schur
 -> K_q -> measured M_hq with explicit neighbour-frame transport -> K_h(k) -> TT
 -> c_micro_spatial_BQG
```

No GR ratio, TT result, observed dispersion or Wilson target may alter the basis.

## Physical bridge still required

The microscopic spatial precursor is not a physical graviton 1PI function. The physical vector requires the separate existing-constraint chain

\[
\boxed{
\{C_A\}
\to P_{phys}/\eta
\to Z_{phys}[J_g]
\to W_{phys}[J_g]
\to\Gamma_{phys}[g]
\to\Gamma_{TT}^{(2)}(\omega,\mathbf k)
\to\mathbf c_{BQG}^{IR}.
}
\]

Only the final connected physical pole data may be called `c_BQG_IR`.
