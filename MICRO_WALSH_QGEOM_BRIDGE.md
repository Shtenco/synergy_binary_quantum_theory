# Exact q=2 Walsh-character → tetrahedral geometry-qubit bridge

## Scope

This result addresses the local carrier part of `MICRO_TO_QGEOM` without inserting a target metric, tetrad, connection or random B-field.

The frozen q=2 microscopic route labels are the four elements of

\[
G=\mathbb Z_2^2=\{00,01,10,11\}.
\]

The three nontrivial real characters of \(G\) define

\[
\Phi(g)=\frac{1}{\sqrt 3}
\big(\chi_{01}(g),\chi_{10}(g),\chi_{11}(g)\big)\in\mathbb R^3.
\]

Character orthogonality gives exactly

\[
\sum_g\Phi(g)=0,
\qquad
\Phi(g)\cdot\Phi(h)=
\begin{cases}
1,&g=h,\\
-1/3,&g\ne h.
\end{cases}
\]

Therefore the four frozen q=2 route labels themselves form the four unit face normals of a regular tetrahedron.

More generally, all \(2^q-1\) nontrivial real characters of \(\mathbb Z_2^q\) embed its \(2^q\) labels as a regular simplex in \(\mathbb R^{2^q-1}\). The q=2 case is distinguished because the character space has exactly three components.

## Quantum carrier

![Walsh-кодирование и тетраэдральный носитель](docs/figures/02-walsh-tetrahedron.svg)


For each derived unit vector \(n_f\), define the face-qubit density matrix

\[
\rho_f=\frac12(I+n_f^i\sigma_i).
\]

No continuous direction is fitted: the Bloch vector is the Walsh vector above. The four face qubits are then reduced with the same exact SU(2) Gauss-singlet projector used by the canonical spatial geometry gate.

The CI gate `scripts/micro_walsh_qgeom_gate.py` reports

```text
flux closure norm                    = 0
regular tetrahedron Gram error       < 1e-14
Gauss-singlet weight                 = 0.22222222222222215 ≈ 2/9
logical geometry Bloch vector        ≈ (0,+1,0)
logical oriented volume              = sqrt(3)/4
reconstructed edge relative spread   = 0
orientation-reversed logical Bloch   ≈ (0,-1,0)
```

Thus the binary-derived coherent face state has nonzero singlet support and projects to an exact logical volume eigenstate in the declared four-spin geometry-qubit representation.

## Global gluing on the selected q=2 PL completion

![Локальная q=2-структура и глобальная 16-cell-комплектация](docs/figures/03-q2-to-16cell.svg)


The boundary of the 4D cross-polytope has 16 tetrahedral cells and 32 triangular faces. Each tetrahedron contains one signed vertex from each of four coordinate axes. Its four faces are canonically labelled by the omitted axis; those four colours are identified with the four q=2 labels.

`scripts/q2_global_face_qubit_gluing_gate.py` checks exactly that:

- every triangular face is incident on exactly two tetrahedra;
- both incident cells assign the same q=2 carrier label to the shared face;
- the dual graph of the 16 tetrahedra is exactly \(Q_4\);
- neighboring tetrahedra differ by one sign bit;
- the parity orientation alternates across every dual edge;
- therefore outward Walsh fluxes cancel pairwise on every shared face.

This produces a globally compatible face-qubit/flux carrier on the selected PL completion without introducing a random B-field.

## What is proved and what is not

**Exact inside the declared construction:**

\[
q=2\ \text{route labels}
\to
\text{regular tetrahedral flux frame}
\to
\text{pure face qubits}
\to
\text{nonzero Gauss-singlet geometry qubit},
\]

plus exact kinematic gluing on the chosen 16-cell completion.

**Still open:** the stronger dynamical claim that the microscopic graph-changing Hamiltonian uniquely generates/selects this quantum-geometric phase and its semiclassical measure from a generic microscopic state. The three Walsh characters are commuting functions on the classical route-label set; the Pauli/Bloch lift is the declared canonical quantum representation, not a derivation of the full noncommutative SU(2)/Peter-Weyl algebra from character multiplication alone.
