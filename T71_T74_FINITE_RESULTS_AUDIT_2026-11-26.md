# T71–T74 finite-result audit — 2026-11-26

This note records a source/artifact audit of four historical OpenCode entries. It distinguishes exact statements inside reduced or finite models from the physical BQG claims they do not establish.

## T71: reduced free-fermion propagation

`BQG_MANY_NODE_FREE_FERMION_MAPPING.md` maps the declared nearest-neighbour logical chain by a local basis rotation and Jordan–Wigner transformation to an open XX/free-fermion chain with one-particle dispersion

\[
\varepsilon(k)=-3\cos k,
\qquad
|v(k)|=|\partial_k\varepsilon|=3|\sin k|\le3.
\]

This is an exact propagation scale for that one-dimensional reduced Hamiltonian. It does not derive a four-dimensional Lorentzian metric or calibrate the number 3 to the physical speed of light. The separate finite graph constraint-move bound is a combinatorial statement and is not that calibration either.

## T72: finite PL-to-causal Galerkin spectrum

The current artifact `verification_results/CAUSAL_PL_CORRESPONDENCE_LOCAL_20261123.json` is marked `passed: true`. Its stored finite Galerkin operators have dimensions 1, 8, and 64 at generations 0, 1, and 2. The recorded effective spectral-dimension maxima are 0, 1.1830050092, and 2.0275436507; the corresponding finite-scale minima include values near zero at the first refinement and 0.2014082467 at generation 2.

The generation-1 and generation-2 partitions are explicitly marked non-equitable. The artifact says these are finite mass-normalized Galerkin spectra; it does not establish a plateau across generations or a controlled continuum limit. Its graph-conductance semantics are PL adjacency multiplicities, not physical QFI conductances. Therefore the rising finite peak is not evidence that the physical BQG spacetime has limiting dimension three.

## T73: exact generation-3 preflight

Fresh run of `python scripts/q2_pl_generation3_preflight.py` passed. From the declared generation-2 f-vector \((1696,10912,18432,9216)\), the barycentric flag-count formula gives

\[
g_0=f_0+f_1+f_2+f_3=40256,\quad
g_1=2f_1+6f_2+14f_3=261440,
\]
\[
g_2=6f_2+36f_3=442368,\quad
g_3=24f_3=221184.
\]

The exact Euler check is zero, and each parent tetrahedron contributes 24 child tetrahedra. This is a combinatorial preflight only; the script does not build the generation-3 quotient spectrum.

## T74: generation-3 graph and causal fibres

Fresh run of `python scripts/q2_pl_generation3_stream.py` passed. It built 221,184 child tetrahedra and 442,368 unique dual edges, with four route classes of 110,592 edges each, 512 nonempty causal fibres, total pushforward mass 16, and one connected component. The direct construction exactly matches the incidence construction at generations 1 and 2.

These checks support the finite graph/fibre construction and mass bookkeeping. They do not show that PL adjacency equals a physical QFI metric, establish an interacting propagator, or prove a spacetime continuum limit. The historical description of this finite combinatorial result as a closed physical C1/C3 layer is stronger than its artifact warrants.

## Relation to the GR/QM bridge

The four entries are not a single implication chain. T73 supports the combinatorial input used by T74; T72 analyzes a finite Galerkin operator built from such graph data; T71 is an independent one-dimensional reduced Hamiltonian result. None supplies the missing common physical-history measure, Lorentzian constraint algebra, physical metric correlator, or calibrated TT pole. They can be premises for later work only after those objects are derived and their maps between sectors are proved.