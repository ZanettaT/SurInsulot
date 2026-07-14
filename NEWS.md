# SurInsulot 0.1.0

* Initial release.
* One vectorised function per surrogate index, grouped by construct domain:
  fasting sensitivity, lipid/anthropometric, OGTT sensitivity, OGIS,
  beta-cell secretion, C-peptide, disposition, and glucagon indices.
* Trapezoidal `auc_trapezoid()` / `iauc_trapezoid()` helpers, `guard_ratio()`,
  `row_mean()`, and unit conversions between mg/dL and mmol/L and between
  insulin uU/mL and pmol/L.
* `compute_surrogates()` derives every estimable index from a data frame of
  fasting and OGTT measurements, mirroring the GRADE and DPP analysis scripts.
* `compute_surrogates()` accepts input-unit declarations (`glucose_unit`,
  `insulin_unit`, `cpeptide_unit`, `trig_unit`, `hdl_unit`, `glucagon_unit`) and
  normalises to the canonical units before computing, so raw data can be in
  mg/dL or mmol/L, uU/mL or pmol/L, nmol/L or ng/mL, etc.
* Yeo-Johnson transform: exported `yeo_johnson()` and `yeo_johnson_lambda()`,
  and a `transform = TRUE` option on `compute_surrogates()` that appends a
  `<index>_yj` column per index (lambda fitted by maximum likelihood),
  reproducing the reports' `*_yj` output.
* HOMA2 (%B, %S, IR) is built in via `homa2_insulin()`,
  `homa2_specific_insulin()`, and `homa2_cpeptide()`, using bilinear
  interpolation over lookup tables from the Oxford DTU HOMA2 Calculator v2.2.4
  (vendored from [homa2calc](https://github.com/ZanettaT/homa2calc)). No
  external dependency is required.
