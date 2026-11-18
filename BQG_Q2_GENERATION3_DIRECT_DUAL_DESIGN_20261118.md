# Q2 generation-3 direct dual adjacency construction

Date: 2026-11-18

## Purpose

Generation 3 contains 221,184 barycentric child tetrahedra. Building a Python
incidence dictionary with 884,736 face entries and 442,368 edge dictionaries
is unnecessary for the dual graph needed by the correspondence and matrix-free
Feshbach gates. The dual adjacency can be generated directly from the
generation-2 tetrahedra and their already verified dual-face pairings.

This note specifies an exact combinatorial construction. It is a design and
proof artifact; it does not claim that the generation-3 exporter or numerical
Feshbach run has already been implemented or executed.

## Child indexing

Let a generation-2 tetrahedron
\[
T=(v_0,v_1,v_2,v_3)
\]
and let \(p=(p_0,p_1,p_2,p_3)\) range over the 24 permutations of its
vertices in a fixed lexicographic order. Its barycentric child flag is
\[
F(T,p)=
(\{p_0\},\{p_0,p_1\},\{p_0,p_1,p_2\},T).
\]
Assign the child a deterministic compact index
\[
i(T,p)=24\,i(T)+r(p),\qquad 0\le r(p)<24.
\]
The child tetrahedron in the generation-3 primal complex is the sorted tuple
of the four face IDs in this flag. The ordered flag must be retained separately:
its order carries the four route-colour slots.

There are exactly
\[
9,216\times24=221,184
\]
child tetrahedra.

## Exhaustive classification of dual neighbours

Two maximal flags in a barycentric subdivision share a codimension-one face
if and only if they differ in exactly one rank of the flag.

### 1. Three within-parent neighbours

For each child \((T,p)\) and each \(k\in\{0,1,2\}\), swap the adjacent
entries \(p_k,p_{k+1}\). The two flags then agree at every rank except
rank \(k\), so they share the other three members of the chain. This gives
one neighbour of route colour \(k\).

Each undirected within-parent adjacency is generated from both endpoints.
Therefore the number of distinct such edges is
\[
E_{\rm in}=\frac{3\cdot221,184}{2}=331,776.
\]

### 2. One across-parent neighbour

Let generation-2 tetrahedra \(T,T'\) share a triangular face \(S\).
For each of the six orderings \(p=(p_0,p_1,p_2)\) of the three vertices
of \(S\), append the unique vertex of \(T\setminus S\) to form a flag
in \(T\), and append the unique vertex of \(T'\setminus S\) to form a
flag in \(T'\).

These two flags have identical first three chain members and differ only
in the final tetrahedron. They are therefore dual neighbours, in route slot
3. Every cross-parent codimension-one flag face is obtained this way, exactly
once, because every parent triangular face has exactly two incident parent
tetrahedra.

The generation-2 dual graph has 18,432 edges, hence
\[
E_{\rm out}=18,432\times6=110,592.
\]

### 3. Total and regularity

The two classes are disjoint and exhaustive:
\[
E_3=E_{\rm in}+E_{\rm out}
=331,776+110,592=442,368.
\]
Every child has exactly three within-parent neighbours (one for each adjacent
transposition) and one across-parent neighbour (the unique parent face
opposite its final permutation entry). Thus the generation-3 dual graph is
4-regular, with
\[
2E_3=4\times221,184.
\]

This direct construction is equivalent to matching each codimension-one
barycentric flag face by a two-entry incidence table, provided the parent
dual-face pairing is itself complete and two-sided.

## Route-colour consistency

For an ordered flag, deleting chain rank \(k\) gives the shared dual face
whose route slot is \(k\).

- Swapping \(p_k,p_{k+1}\), for \(k=0,1,2\), preserves all other chain
  members and pairs two flags with the same deleted rank \(k\).
- Across a parent face, the two appended flags share the first three chain
  members; the deleted rank is 3 on both sides.

Consequently, the route colour agrees across every generated dual edge.
The construction does not infer route colours from geometric coordinates.

## Causal-fibre labels

Once the generation-2 correspondence data are available, the generation-3
child associated with parent \(T\) and permutation \(p\) receives:

1. route label \(c=\mathrm{faceColour}_T(p_0)\), from the parent face
   opposite the first permutation entry;
2. side \(s=\mathrm{faceSide}(T,T\setminus\{p_3\})\), from the
   already audited oriented parent-face rule;
3. causal child edge \((T_{\rm causal},c,s)\), where
   \(T_{\rm causal}=\Phi_2(T)\).

This is the same rule used by the generation-0-to-2 correspondence exporter.
It should be stored as a compact integer fibre-label array of length 221,184,
not as per-child Python dictionaries.

## Compact implementation plan

Use arrays rather than per-face or per-edge Python objects:

- child ordered flags / sorted tetrahedra: uint32[221184, 4];
- parent ID and permutation rank: uint32[221184] each, or recoverable
  from the deterministic child index;
- causal fibre label: uint16[221184] (512 fibres);
- dual edge endpoints: uint32[442368, 2];
- dual route colour: uint8[442368].

Generate within-parent edges by adjacent permutation swaps, retaining only
the canonical endpoint order. Generate across-parent edges from each
generation-2 dual edge and the six permutations of its shared face. Validate:

1. every parent dual face is paired exactly twice and route-colour-consistently;
2. every generated child has degree exactly four;
3. the exact edge split is 331,776 within-parent plus 110,592 across-parent;
4. all 442,368 edge route colours agree at both endpoints;
5. the graph is connected;
6. every child has exactly one causal fibre label and every expected causal
   child fibre is nonempty;
7. pushforward parent mass is conserved under the eight-child causal split.

For Feshbach, construct CSR directly from the endpoint array. Do not create a
JSON payload containing 221,184 vertex dictionaries and 442,368 edge
dictionaries. If a stable digest is required, hash canonical fixed-width
integer records in deterministic order rather than serializing Python objects.

## Scope boundary

The formulas above prove the counts and adjacency construction conditional on
the generation-2 parent complex and its two-sided dual pairing. They do not
establish that the causal correspondence remains surjective at generation 3,
that the pushforward measure is physically selected, that PL adjacency is a
QFI conductance, or that a continuum/physical spacetime limit exists.

The next implementation gate is to realize this compact construction, compare
its generation-0/1/2 restriction against the existing incidence exporter,
then run the generation-3 causal-fibre and matrix-free Feshbach checks.
