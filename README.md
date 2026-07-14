<!-- README.md is written by hand; edit here directly. -->

# SurInsulot <img src="man/figures/logo.png" align="right" height="130" alt="SurInsulot hex logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/ZanettaT/SurInsulot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ZanettaT/SurInsulot/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Compute **Sur**rogates for **Insul**in sensitivity, insulin secretion, or beta-cell function, from fasting bloods and the oral
glucose tolerance test (OGTT).

Every index is its own function, but there's also a single `compute_surrogates()` wrapper which can derive the
whole panel from a data frame.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ZanettaT/SurInsulot")
```

SurInsulot has no hard dependencies and bundles everything it needs, including
the HOMA2 and OGIS calculations.

## Quick start

Call an index directly and pay mind to the units in each function's help page:

```r

library(SurInsulot)

homa_ir(insulin = 12, glucose = 5.56) # fasting insulin uU/mL, glucose mmol/L
#> [1] 2.963

quicki(insulin = 12, glucose = 100)  # glucose mg/dL
#> [1] 0.3248

igi(insulin_0 = 12, insulin_t = 90, glucose_0 = 5.56, glucose_t = 8.89)
#> [1] 23.4                         # fasting insulin uU/mL, glucose mmol/L

```

Or derive everything estimable from a data frame in one call. Columns are
computed only when their inputs are present, so the same call works on an 
OGTT dataset or a fasting-only one:

```r

df <- data.frame(
  g000 = c(95, 110), g030 = c(160, 190), g120 = c(140, 170),  # glucose mg/dL
  i000 = c(10, 18),  i030 = c(80, 60),                        # insulin uU/mL
  trig = c(150, 200), chdl = c(45, 38), bmi = c(29, 34), sex = c(1, 2)
)

out <- compute_surrogates(df)

#> SurInsulot: computed 22 surrogates
out[, c("homa_ir", "quicki", "tyg", "matsuda_30", "igi_30")]

```

`compute_surrogates()` expects the short column names
(`g000`…`g120`, `i000`…`i120`, `cp000`…`cp120`, `glcg000/glcg030`,
`trig`, `chdl`, `bmi`, `weight_kg`, `waist_cm`, `hip_cm`, `height_cm`, `whr`,
`whtr`, `hba1c_perc`, `adipon`, `leptin`, `pin`, `sex`, `age`, `hypertension`).
See `?compute_surrogates` for the full schema and units.

## Units

By default `compute_surrogates()` expects **glucose in mg/dL** and **insulin in
uU/mL** (= mU/L), and converts internally where an index needs mmol/L or pmol/L.
If your data are in other units, just declare them — inputs are normalised
before anything is computed:

```r
compute_surrogates(df,
  glucose_unit  = "mmol/L",   # or "mg/dL" (default)
  insulin_unit  = "pmol/L",   # or "uU/mL" (default)
  cpeptide_unit = "ng/mL",    # or "nmol/L" (default), "pmol/L"
  trig_unit     = "mmol/L",   # or "mg/dL" (default)
  hdl_unit      = "mmol/L",   # or "mg/dL" (default)
  glucagon_unit = "pg/mL")    # or "pmol/L" (default)
```

When calling the individual index functions directly, follow the units in each
help page. The conversion helpers are also exported:

```r
mgdl_to_mmol_glucose(90)   # 5.0
insulin_uU_to_pmol(10)     # 69.45   (1 uU/mL = 6.945 pmol/L)
```

## Index reference

Notation: `I0`,`G0` = fasting insulin/glucose; `It`,`Gt` = value at *t* min;
`Imean`,`Gmean` = trapezoidal mean across the OGTT; `ln` = natural log,
`log10` = base-10; `wt` = weight (kg).

### Fasting sensitivity / resistance

| Function | Index | Formula | Units |
|---|---|---|---|
| `raynaud()` | Raynaud | `40 / I0` | I uU/mL |
| `homa_ir()` | HOMA-IR | `I0 * G0 / 22.5` | I uU/mL, G mmol/L |
| `firi()` | FIRI | `I0 * G0 / 25` | I uU/mL, G mmol/L |
| `quicki()` | QUICKI | `1 / (log10 I0 + log10 G0)` | I uU/mL, G mg/dL |
| `insulin_glucose_ratio()` | Insulin:glucose | `I0 / G0` | I uU/mL, G mmol/L |
| `glucose_insulin_ratio()` | Glucose:insulin | `G0 / I0` | I uU/mL, G mmol/L |
| `bennett_isi()` | Bennett ISI | `1 / (ln I0 * ln G0)` | I uU/mL, G mmol/L |
| `mcauley()` | McAuley | `exp(2.63 - 0.28 ln I0 - 0.31 ln TG)` | I uU/mL, TG mmol/L |
| `ohkura()` | Ohkura | `20 / (CP0 * G0)` | CP nmol/L, G mmol/L |
| `homa_ad()` | HOMA-AD | `I0 * G0 / (22.5 * adiponectin)` | I uU/mL, G mmol/L, adipo ug/mL |
| `belfiore_basal()` | Belfiore basal | `2 / ((I0/Iref)(G0/Gref) + 1)` | I uU/mL, G mmol/L |
| `egdr_bmi()` | eGDR (BMI) | `19.02 - 0.22 BMI - 3.26 HTN - 0.61 HbA1c` | mg/kg/min |
| `egdr_waist()` | eGDR (waist) | `21.158 - 0.09 WC - 3.407 HTN - 0.551 HbA1c` | mg/kg/min |
| `egdr_whr()` | eGDR (WHR) | `24.31 - 12.22 WHR - 3.29 HTN - 0.57 HbA1c` | mg/kg/min |

`1 / homa_ir()` gives the reciprocal-HOMA sensitivity index.

### Lipid / anthropometric

| Function | Index | Formula | Units |
|---|---|---|---|
| `tg_hdl()` | TG:HDL | `TG / HDL` | mg/dL |
| `tyg()` | TyG | `ln(TG * G0 / 2)` | mg/dL |
| `tyg_adjusted()` | TyG-BMI/WC/WHtR/WHR | `TyG * adiposity` | — |
| `mets_ir()` | METS-IR | `ln(2 G0 + TG) * BMI / ln(HDL)` | mg/dL, BMI kg/m² |
| `mets_vf()` | METS-VF | `4.466 + 0.011 ln(METS-IR)³ + 3.239 ln(WHtR)³ + 0.319 male + 0.594 ln(age)` | — |
| `spise()` | SPISE | `600 HDL^0.185 / (TG^0.2 BMI^1.338)` | mg/dL, BMI kg/m² |
| `vai()` | VAI | sex-specific (see `?vai`) | TG, HDL mmol/L |
| `lap()` | LAP | `(WC - k) * TG`, k = 65 (M) / 58 (F) | TG mmol/L |
| `leptin_adiponectin_ratio()` | LAR | `leptin / adiponectin` | ng/mL : ug/mL |

### OGTT sensitivity

| Function | Index | Formula | Units |
|---|---|---|---|
| `matsuda()` | Matsuda | `10000 / sqrt(G0 * I0 * Gmean * Imean)` | I uU/mL, G mg/dL |
| `matsuda_from_curve()` | Modified Matsuda | Matsuda with trapezoidal means over 0–*t* | I uU/mL, G mg/dL |
| `gutt_isi()` | Gutt ISI(0,120) | `(75000 + (G0 - G120) 0.19 wt) / (120 Gmean log10 Imean)` | I uU/mL, G mg/dL |
| `cederholm()` | Cederholm | `(75000 + (G0 - G120) 1.15·180·0.19 wt) / (120 Gmean log10 Imean)` | I uU/mL, G mmol/L |
| `stumvoll_isi()` | Stumvoll ISI | `0.226 - 0.0032 BMI - 6.45e-5 I120 - 0.00375 G90` | I pmol/L, G mmol/L |
| `stumvoll_isi_demo()` | Stumvoll (demographic) | `0.222 - 0.00333 BMI - 7.79e-5 I120 - 4.22e-4 age` | I pmol/L |
| `stumvoll_isi_bmi_independent()` | Stumvoll (BMI-free) | `0.156 - 4.59e-5 I120 - 3.21e-4 I0 - 0.00541 G120` | I pmol/L, G mmol/L |
| `stumvoll_mcr()` | Stumvoll MCR | `18.8 - 0.271 BMI - 0.0052 I120 - 0.27 G90` | I pmol/L, G mmol/L |
| `avignon_sib()` | Avignon Sib | `1e8 / (I0 G0 · 150 wt)` | I uU/mL, G mmol/L |
| `avignon_si2h()` | Avignon Si2h | `1e8 / (I120 G120 · 150 wt)` | I uU/mL, G mmol/L |
| `avignon_sim()` | Avignon SiM | `(0.137 Sib + Si2h) / 2` | — |
| `ogis()` | OGIS (2 h model) | closed-form Mari model (see `?ogis`) | I pmol/L, G mmol/L |
| `belfiore_ogtt()` | Belfiore OGTT | `2 / ((AI/AIref)(AG/AGref) + 1)` | 0–1–2 h areas |

### Secretion / beta-cell

| Function | Index | Formula | Units |
|---|---|---|---|
| `homa_b()` | HOMA-%B | `20 I0 / (G0 - 3.5)` | I uU/mL, G mmol/L |
| `igi()` | Insulinogenic index | `(It - I0) / (Gt - G0)` | I uU/mL, G mmol/L |
| `cir()` | Corrected insulin response | `It / (Gt (Gt - 70))`, else NA | I uU/mL, G mg/dL |
| `kadowaki()` | Kadowaki | `(I30 - I0) / G30` | I uU/mL, G mg/dL |
| `stumvoll_first_phase()` | Stumvoll 1st phase | `1283 + 1.829 I30 - 138.7 G30 + 3.772 I0` | I pmol/L, G mmol/L |
| `stumvoll_second_phase()` | Stumvoll 2nd phase | `287 + 0.416 I30 - 26.07 G30 + 0.926 I0` | I pmol/L, G mmol/L |
| `bigtt_air_30_120()` / `bigtt_air_60_120()` | BIGTT-AIR | see `?bigtt_air` | I uU/mL, G mmol/L |

### HOMA2 (Oxford DTU model, built in)

HOMA2 has no closed form, so these interpolate over reference tables from the
DTU HOMA2 Calculator v2.2.4 (bundled). Each returns a data frame of `homa2_b`
(%B), `homa2_s` (%S), and `homa2_ir` (= 100/%S); out-of-range inputs give `NA`.

| Function | Inputs | Units |
|---|---|---|
| `homa2_insulin()` | glucose + non-specific insulin | mmol/L, pmol/L |
| `homa2_specific_insulin()` | glucose + proinsulin-free insulin | mmol/L, pmol/L |
| `homa2_cpeptide()` | glucose + C-peptide | mmol/L, nmol/L |

### C-peptide

| Function | Index | Formula | Units |
|---|---|---|---|
| `cpeptide_glucose_ratio()` | C-peptide:glucose | `CPt / Gt` | CP nmol/L, G mmol/L |
| `cpi()` | C-peptide index | `(CPt - CP0) / (Gt - G0)` | CP nmol/L, G mmol/L |
| `insulin_cpeptide_ratio()` | Insulin:C-peptide (molar) | `I0 / (CP0 * 1000)` | I pmol/L, CP nmol/L |
| `cpeptide_insulin_auc_ratio()` | C-peptide:insulin (AUC, molar) | `(AUC_CP * 1000) / (AUC_I * 6.945)` | — |

### Disposition & glucagon

| Function | Index | Formula |
|---|---|---|
| `disposition_index()` | Oral disposition index & variants (oDI, ISSI-2) | `secretion * sensitivity` |
| `glucagon_insulin_ratio()` | Glucagon:insulin (molar) | `GCG0 / I0` (both pmol/L) |
| `glucagon_suppression()` | Glucagon suppression | `GCG0 - GCG30` or its % |

### Transforms

Many indices are skewed. Pass `transform = TRUE` to `compute_surrogates()` to
append a Yeo-Johnson-normalised `<index>_yj` column for every index (λ fitted per
index by maximum likelihood), reproducing the `*_yj` outputs of the reports; the
fitted λ are stored in `attr(result, "surinsulot_yj_lambda")`. `yeo_johnson()`
and `yeo_johnson_lambda()` are also exported for standalone use.

### Helpers

`auc_trapezoid()`, `iauc_trapezoid()` (total / incremental trapezoidal AUC),
`guard_ratio()` (NA-safe division), `row_mean()`, and the `*_to_*` unit
converters.

## Notes on the formulas

- **HOMA2** (`homa2_ir`, `homa2_b`, `homa2_s`, and C-peptide variants) is built
  in via `homa2_insulin()` / `homa2_cpeptide()` — bilinear interpolation over
  lookup tables from the Oxford DTU HOMA2 Calculator v2.2.4, since the model has
  no closed form. `compute_surrogates()` converts fasting insulin to pmol/L at
  6.945 before the lookup. The tables are vendored from
  [homa2calc](https://github.com/ZanettaT/homa2calc).
- **Stumvoll first/second-phase secretion** expect insulin in **pmol/L**,
  matching Stumvoll (2000) and the analysis code (the GRADE data dictionary's
  uU/mL label is a documentation slip).
- **`insulin_cpeptide_ratio()`** returns a *true molar ratio* (~0.1–0.2). The
  GRADE script computes `(I_pmol / CP_nmol) * 1000`, which is 10⁶× larger; that
  is rank-preserving (so it does not change the Spearman analyses in the report)
  but is not a molar ratio. Multiply by `1e6` to reproduce the script's raw
  numbers.
- **Gutt / Cederholm** use a base-10 log of mean insulin, as in the source
  scripts.

## Provenance

The formulas, unit conventions, and column names come from two analysis
pipelines over the GRADE and DPP (Diabetes Prevention Program) cohorts. This
package factors that logic into tested, reusable functions.

## License

MIT © Zanetta Toomata. See [LICENSE.md](LICENSE.md).
