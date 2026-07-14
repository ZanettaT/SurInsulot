# HOMA2 calculator (bilinear interpolation over the Oxford DTU reference grids).
#
# Vendored from the homa2calc package (https://github.com/ZanettaT/homa2calc),
# by the same author, so SurInsulot computes HOMA2 with no external dependency.
# The lookup grids live in R/homa2-tables.R. They derive from the HOMA2
# Calculator v2.2.4 validation dataset, (c) Diabetes Trials Unit, University of
# Oxford (https://www.dtu.ox.ac.uk/homacalculator/). Unlike the other indices,
# HOMA2 has no published closed form; the DTU model's constants were never
# released, so a table lookup is the only faithful implementation.

# ---- internal bilinear interpolator ---------------------------------------

.homa2_interp_single <- function(glucose, hormone, gluc_axis, horm_axis,
                                 b_mat, s_mat) {
  n_g <- length(gluc_axis)
  n_h <- length(horm_axis)

  # out-of-range -> NA
  if (is.na(glucose) || is.na(hormone) ||
      glucose < gluc_axis[1] || glucose > gluc_axis[n_g] ||
      hormone < horm_axis[1] || hormone > horm_axis[n_h]) {
    return(c(pct_b = NA_real_, pct_s = NA_real_, ir = NA_real_))
  }

  # bounding indices
  g1 <- max(which(gluc_axis <= glucose))
  g2 <- min(which(gluc_axis >= glucose))
  h1 <- max(which(horm_axis <= hormone))
  h2 <- min(which(horm_axis >= hormone))

  # fractional positions
  tx <- if (g2 == g1) 0 else (glucose - gluc_axis[g1]) / (gluc_axis[g2] - gluc_axis[g1])
  ty <- if (h2 == h1) 0 else (hormone - horm_axis[h1]) / (horm_axis[h2] - horm_axis[h1])

  bilerp <- function(mat) {
    (1 - tx) * (1 - ty) * mat[g1, h1] +
      tx       * (1 - ty) * mat[g2, h1] +
      (1 - tx) *      ty  * mat[g1, h2] +
      tx       *      ty  * mat[g2, h2]
  }

  pct_b <- bilerp(b_mat)
  pct_s <- bilerp(s_mat)
  ir    <- 100 / pct_s

  c(pct_b = round(pct_b, 1), pct_s = round(pct_s, 1), ir = ir)
}

.homa2_vectorised <- function(glucose, hormone, gluc_axis, horm_axis,
                              b_mat, s_mat) {
  n <- length(glucose)
  out_b <- numeric(n)
  out_s <- numeric(n)
  out_ir <- numeric(n)
  for (i in seq_len(n)) {
    res <- .homa2_interp_single(glucose[i], hormone[i],
                                gluc_axis, horm_axis, b_mat, s_mat)
    out_b[i]  <- res["pct_b"]
    out_s[i]  <- res["pct_s"]
    out_ir[i] <- res["ir"]
  }
  data.frame(homa2_b = out_b, homa2_s = out_s, homa2_ir = out_ir)
}


#' HOMA2 beta-cell function, sensitivity, and insulin resistance
#'
#' The updated Homeostasis Model Assessment (HOMA2) estimates steady-state
#' beta-cell function (`%B`), insulin sensitivity (`%S`), and insulin resistance
#' (`IR = 100 / %S`) from a single fasting glucose plus fasting insulin or
#' C-peptide. Unlike HOMA-IR, HOMA2 has **no published closed-form equation** -
#' the Oxford DTU model's constants were never released - so these functions use
#' bilinear interpolation over reference grids extracted from the HOMA2
#' Calculator v2.2.4 validation dataset. Grid-point values reproduce the official
#' DTU calculator exactly; interpolated values are accurate to rounding.
#'
#' * `homa2_insulin()` - non-specific (conventional RIA) insulin, which
#'   cross-reacts with proinsulin.
#' * `homa2_specific_insulin()` - proinsulin-free (specific) insulin.
#' * `homa2_cpeptide()` - fasting C-peptide, reflecting secretion without
#'   hepatic insulin extraction.
#'
#' Inputs outside the valid range return `NA`.
#'
#' @param glucose Fasting plasma glucose in **mmol/L** (valid range 3.0-25.0).
#' @param insulin Fasting plasma insulin in **pmol/L** (valid range 20-400).
#'   From uU/mL, multiply by 6.945 (or 6.0 for the older DTU convention).
#' @param specific_insulin Fasting specific (proinsulin-free) insulin in
#'   **pmol/L** (valid range 20-300).
#' @param cpeptide Fasting C-peptide in **nmol/L** (valid range 0.2-3.5). From
#'   ng/mL divide by 3.02; from pmol/L divide by 1000.
#'
#' @return A data frame with one row per input and columns `homa2_b` (%B),
#'   `homa2_s` (%S) and `homa2_ir` (= 100 / %S). Out-of-range rows are `NA`.
#'
#' @references
#' Levy JC, Matthews DR, Hermans MP (1998) Correct homeostasis model assessment
#' (HOMA) evaluation uses the computer program. *Diabetes Care* 21:2191-2192.
#' Wallace TM, Levy JC, Matthews DR (2004) Use and abuse of HOMA modeling.
#' *Diabetes Care* 27:1487-1495.
#' HOMA2 Calculator v2.2.4, Diabetes Trials Unit, University of Oxford,
#' <https://www.dtu.ox.ac.uk/homacalculator/>.
#'
#' @examples
#' homa2_insulin(glucose = 5.0, insulin = 60)
#' homa2_cpeptide(glucose = 5.0, cpeptide = 0.7)
#' # vectorised; convert insulin uU/mL -> pmol/L first
#' homa2_insulin(c(5.0, 6.0), insulin_uU_to_pmol(c(8, 14)))
#'
#' @name homa2
NULL

#' @rdname homa2
#' @export
homa2_insulin <- function(glucose, insulin) {
  glucose <- as.numeric(glucose)
  insulin <- as.numeric(insulin)
  stopifnot(length(glucose) == length(insulin))
  .homa2_vectorised(glucose, insulin,
                    homa2_ins_gluc_axis, homa2_ins_horm_axis,
                    homa2_ins_b_mat, homa2_ins_s_mat)
}

#' @rdname homa2
#' @export
homa2_specific_insulin <- function(glucose, specific_insulin) {
  glucose <- as.numeric(glucose)
  specific_insulin <- as.numeric(specific_insulin)
  stopifnot(length(glucose) == length(specific_insulin))
  .homa2_vectorised(glucose, specific_insulin,
                    homa2_spec_gluc_axis, homa2_spec_horm_axis,
                    homa2_spec_b_mat, homa2_spec_s_mat)
}

#' @rdname homa2
#' @export
homa2_cpeptide <- function(glucose, cpeptide) {
  glucose <- as.numeric(glucose)
  cpeptide <- as.numeric(cpeptide)
  stopifnot(length(glucose) == length(cpeptide))
  .homa2_vectorised(glucose, cpeptide,
                    homa2_cpep_gluc_axis, homa2_cpep_horm_axis,
                    homa2_cpep_b_mat, homa2_cpep_s_mat)
}
