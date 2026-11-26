# HDA targeted regression audit — 2026-11-26

This note supplements the HDA status in `BQG_GR_QM_GRAVITON_INFOTON_INTEGRATION_AUDIT.md`.

## Fresh targeted run

Command:

```powershell
python -c "import torch; torch.set_num_threads(1); from scripts.dewitt_hda_uniqueness_gate import one; rows=[one(c,n=7) for c in (0.5,0.4)]; print([{k:r[k] for k in ('c','analytic_formula_relative_error','HDA_defect_against_GR_rhs','extra_integral_I')} for r in rows])"
```

Observed output:

| `c` | analytic formula relative error | HDA defect vs. GR RHS | `I` |
|---:|---:|---:|---:|
| 0.5 | `4.649216241189818e-09` | `4.649216241189818e-09` | `-0.0016898362327741768` |
| 0.4 | `5.055917240485696e-09` | `0.028459416419483966` | `-0.0016898362327741768` |

This reproduces the analytic bracket formula numerically on the declared `L=7` smooth state for two values of `c`, and shows a nonzero HDA defect away from `c=1/2`. It is a finite spectral regression, not the proof of the continuum identity and not a test of the microscopic BQG constraint algebra.

## Resolution sensitivity

The same targeted calculation at `L=4` gave a formula residual/HDA defect `5.84610240590843e-4` at `c=0.5` and formula residual `4.53494658845964e-4` at `c=0.4`. The underresolved grid does not reproduce the `L=7` precision. This is a warning against treating the small-grid number as a closure failure or as continuum evidence.

The earlier full nine-value `dewitt_hda_uniqueness_gate.py` attempt timed out. The two-row `L=7` run above is the only newly completed HDA regression in this audit; it does not establish the full parameter sweep reported in `DEWITT_HDA_UNIQUENESS.md`.

## Mathematical boundary

Within the documented local ADM ansatz, the source states

\[
\{H_c[N],H_c[M]\}=D[\beta]+4(c-1/2)I[N,M].
\]

For generic canonical data and arbitrary lapses, this selects `c=1/2` for standard HDA closure. Its application to BQG remains conditional on deriving the same local constraint ansatz, canonical structure, and HDA from the microscopic/refined theory.