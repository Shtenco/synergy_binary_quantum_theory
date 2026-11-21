# Peter–Weyl Lorentzian middle-prefix result

Status: **MEASURED COMPLETE 12-PREFIX PARTITION FOR THE PREREGISTERED FIRST LOGICAL INPUT. Six ordered middle-prefix classes are pathwise zero and six are nonzero. Outer continuation remains required for the six live classes.**

The canonical complete middle-prefix evidence is `lorentzian-all-middle-prefixes` run `33971440632`, head `d1907fc723517ac405133b23e6add0e816e270a6`, at source node `v=0`, logical input index `0` (all local `K=0`), sine-Hermitian Euclidean ordering and `Jmax=7/2`.

The earlier one-term run `33962645194` remains independent evidence for the six pathwise-zero classes, while middle-prefix run `33964666509` independently measured the six live classes. The batched all-prefix run then recomputed all twelve prefixes in one shared-cache session and compared the six live sparse middle states against those independent artifacts.

The frozen ordered Lorentzian term is

\[
T_{abc}=\operatorname{Tr}_{aux}[C_a(K)C_b(K)C_c(V)],
\qquad K=[V,H_E^{sine}].
\]

For fixed ordered pair `(b,c)`, define

\[
\boxed{\Xi_{bc}^{ijk}=C_b(K)_{jk}\,C_c(V)_{ki}|\psi_0\rangle}.
\]

This state is independent of the outer edge `a`. Therefore a pathwise-zero prefix kills both full triples that share `(b,c)` before the expensive outer `C_a(K)` action.

## 1. Complete measured partition

Using the frozen neighbor-order pair enumeration, the all-prefix packet reports

```text
zero_prefix_indices    = [2,5,8,9,10,11]
nonzero_prefix_indices = [0,1,3,4,6,7]
zero_prefix_count      = 6
nonzero_prefix_count   = 6
```

Equivalently,

\[
\boxed{
\text{zero prefixes}
=\{(3,4),(2,4),(1,4),(4,3),(4,2),(4,1)\}
}
\]

and

\[
\boxed{
\text{live prefixes}
=\{(1,2),(1,3),(2,1),(2,3),(3,1),(3,2)\}.
}
\]

The six zero classes imply exactly twelve of the twenty-four epsilon-ordered triples vanish before the outer `C(K)` hit. The six live classes imply the other twelve triples cannot be discarded at the middle stage.

Therefore the previous possibility

```text
all 24 ordered triples die before outer C(K)
```

is **falsified for this first logical input**.

This does **not** yet imply

\[
H_L^{raw}|\psi_0\rangle\neq0,
\]

because the twelve live outer continuations may still cancel after the final `C_a(K)` actions and epsilon-signed sum.

## 2. Pathwise-zero controls

The earlier one-term artifacts resolve the six zero classes more finely:

| middle prefix `(b,c)` | two full triples `(a,b,c)` | `C(V)` support over 8 paths | middle `C(K)` support over 8 paths | result |
|---|---|---|---|---|
| `(3,4)` | `(2,3,4)`, `(1,3,4)` | `0,0,0,0,0,0,0,0` | `0,0,0,0,0,0,0,0` | prefix zero |
| `(2,4)` | `(3,2,4)`, `(1,2,4)` | `0,0,0,0,0,0,0,0` | `0,0,0,0,0,0,0,0` | prefix zero |
| `(1,4)` | `(3,1,4)`, `(2,1,4)` | `0,0,0,0,0,0,0,0` | `0,0,0,0,0,0,0,0` | prefix zero |
| `(4,3)` | `(2,4,3)`, `(1,4,3)` | `2,1,2,1,1,2,1,2` | `0,0,0,0,0,0,0,0` | middle `C_4(K)` kills nonzero `C(V)` |
| `(4,2)` | `(3,4,2)`, `(1,4,2)` | `5,3,5,3,3,5,3,5` | `0,0,0,0,0,0,0,0` | middle `C_4(K)` kills nonzero `C(V)` |
| `(4,1)` | `(3,4,1)`, `(2,4,1)` | `5,3,5,3,3,5,3,5` | `0,0,0,0,0,0,0,0` | middle `C_4(K)` kills nonzero `C(V)` |

The last three rows are nontrivial selection-rule controls because the right-most volume commutator is nonzero before the middle `C_4(K)` annihilates it.

## 3. Batched-vs-independent regression

Run `33971440632` computed all twelve prefixes in one shared-cache session with

```text
unique C(V) states                 = 16
naive separate-prefix C(V) requests = 96
explicit unique C(V) evaluations    = 16
```

and then compared the six live sparse middle states path-by-path with the independently produced artifacts from run `33964666509`.

The measured worst discrepancies were

\[
\boxed{\max |\Delta a|\approx4.30\times10^{-16}}
\]

and

\[
\boxed{\max \|\Delta\|_2/\|\Xi\|_2\approx3.31\times10^{-15}}.
\]

Thus the cache-sharing optimization reproduces the independent middle-prefix results at floating-point roundoff and does not alter the frozen operator.

## 4. Exact outer-continuation reduction

Only the twelve live triples require the final outer action. Linearity permits increasingly strong exact regrouping without changing the operator:

\[
\sum_k C_a(K)_{ij}\Xi_{bc}^{ijk}
=C_a(K)_{ij}\sum_k\Xi_{bc}^{ijk},
\]

and, more strongly,

\[
\boxed{
\sum_{b,c,k}\epsilon_{abc}\,C_a(K)_{ij}\Xi_{bc}^{ijk}
=
C_a(K)_{ij}
\left[\sum_{b,c,k}\epsilon_{abc}\Xi_{bc}^{ijk}\right].
}
\]

The latter groups all live prefixes by frozen outer edge `a` and auxiliary indices `(i,j)`. It requires at most sixteen final `C(K)` calls for the entire first logical column, rather than ninety-six ungrouped calls or forty-eight prefix-grouped calls. This is pure operator linearity: no tolerance pruning, symmetry reconstruction, beta fitting or physical assumption is introduced.

The outer result must still be reported only as one of

```text
FULL_RAW_HL_COLUMN_ZERO
FULL_RAW_HL_COLUMN_NONZERO_WITH_LOGICAL_RETURN
FULL_RAW_HL_COLUMN_NONZERO_LOGICAL_RETURN_ZERO
```

for input index `0`.

## 5. Claim boundary

This is a finite Peter–Weyl microscopic result for one source and one boundary input. It closes the middle-prefix classification only. It does not prove a global property of `H_L`, does not choose the Hermitian physical Lorentzian convention, does not certify the quantum HDA residual, and does not emit the enlarged physical projector, `W_phys`, a scalar pole, dark matter or dark energy.
