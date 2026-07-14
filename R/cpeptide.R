#' C-peptide secretion and hepatic-extraction indices
#'
#' C-peptide is co-secreted with insulin but escapes hepatic first-pass
#' extraction, so C-peptide indices track beta-cell output, and insulin:C-peptide
#' ratios track hepatic insulin extraction.
#'
#' \tabular{ll}{
#'   `cpeptide_glucose_ratio()` \tab \eqn{CP_t / G_t} (CP nmol/L, G mmol/L) \cr
#'   `cpi()`                     \tab \eqn{(CP_t - CP_0) / (G_t - G_0)} \cr
#'   `insulin_cpeptide_ratio()`  \tab \eqn{I_0 / (CP_0 \times 1000)} (molar) \cr
#'   `cpeptide_insulin_auc_ratio()` \tab \eqn{(AUC_{CP} \times 1000) / (AUC_I \times 6.945)}
#' }
#'
#' `cpeptide_glucose_ratio()` and `cpi()` (the C-peptide analogue of the
#' insulinogenic index) are timepoint-agnostic. `insulin_cpeptide_ratio()`
#' expects **fasting insulin in pmol/L** and **C-peptide in nmol/L** and returns
#' the dimensionless molar ratio: C-peptide is scaled by 1000 to pmol/L before
#' dividing, giving values around 0.1-0.2 in health (insulin and C-peptide are
#' co-secreted, but C-peptide is cleared more slowly so its molar concentration
#' is higher). `cpeptide_insulin_auc_ratio()` converts the C-peptide area
#' (nmol/L*min) and insulin area (uU/mL*min) to a molar output ratio.
#'
#' @note The GRADE analysis script computes this ratio as
#'   `(insulin_pmol / cpeptide_nmol) * 1000`, which is 10^6 times larger than a
#'   true molar ratio (the 1000 factor is on the wrong side). Because it is
#'   rank-preserving this does not change the Spearman analyses in that report,
#'   but SurInsulot returns the corrected molar ratio. To reproduce the script's
#'   raw values exactly, multiply this result by `1e6`.
#'
#' @param cpeptide_0,cpeptide_t Fasting and post-load C-peptide, nmol/L.
#' @param glucose_0,glucose_t Fasting and post-load glucose, mmol/L.
#' @param insulin_0 Fasting insulin, pmol/L (for `insulin_cpeptide_ratio()`).
#' @param auc_cpeptide Trapezoidal C-peptide area, nmol/L*min.
#' @param auc_insulin Trapezoidal insulin area, uU/mL*min.
#'
#' @return A numeric vector the length of the recycled inputs.
#'
#' @references
#' Polonsky KS, Rubenstein AH (1984) *Diabetes* 33:486-494 (C-peptide kinetics).
#'
#' @examples
#' cpeptide_glucose_ratio(cpeptide_t = 1.2, glucose_t = 8.5)
#' cpi(cpeptide_0 = 0.6, cpeptide_t = 1.8, glucose_0 = 5.3, glucose_t = 8.9)
#' insulin_cpeptide_ratio(insulin_0 = 70, cpeptide_0 = 0.6)
#'
#' @name cpeptide
NULL

#' @rdname cpeptide
#' @export
cpeptide_glucose_ratio <- function(cpeptide_t, glucose_t) {
  cpeptide_t / glucose_t
}

#' @rdname cpeptide
#' @export
cpi <- function(cpeptide_0, cpeptide_t, glucose_0, glucose_t) {
  (cpeptide_t - cpeptide_0) / (glucose_t - glucose_0)
}

#' @rdname cpeptide
#' @export
insulin_cpeptide_ratio <- function(insulin_0, cpeptide_0) {
  insulin_0 / (cpeptide_0 * 1000)
}

#' @rdname cpeptide
#' @export
cpeptide_insulin_auc_ratio <- function(auc_cpeptide, auc_insulin) {
  (auc_cpeptide * 1000) / (auc_insulin * 6.945)
}
