#' Matsuda whole-body insulin-sensitivity index
#'
#' The Matsuda index summarises whole-body insulin sensitivity from an OGTT.
#' `matsuda()` is the faithful formula given the fasting values and the mean
#' glucose and insulin across the curve:
#' \deqn{10000 / \sqrt{G_0 \times I_0 \times \bar{G} \times \bar{I}}.}
#' `matsuda_from_curve()` is a convenience wrapper that computes the means for
#' you as trapezoidal time-averages over the sampled window, reproducing the
#' "modified Matsuda" variants (0-15, 0-30, 0-60, 0-90) used when only part of
#' the curve is available.
#'
#' Glucose is in **mg/dL** and insulin in **uU/mL**. The mean is the
#' trapezoidal area under the curve divided by the covered time span, matching
#' the reference implementation (this is the "modified" Matsuda mean rather than
#' the simple arithmetic mean of the original paper).
#'
#' @param glucose_0 Fasting (0 min) glucose, mg/dL.
#' @param insulin_0 Fasting (0 min) insulin, uU/mL.
#' @param glucose_mean,insulin_mean Mean glucose (mg/dL) and insulin (uU/mL)
#'   across the OGTT window.
#' @param times Numeric vector of OGTT sampling times in minutes, starting at 0.
#' @param glucose,insulin Numeric matrices (subjects in rows, timepoints in
#'   columns) or single numeric vectors, aligned to `times`. Glucose mg/dL,
#'   insulin uU/mL.
#' @param upto Optional upper time bound (minutes). Timepoints at or before
#'   `upto` are used for the mean; `NULL` (default) uses the whole curve.
#'
#' @return A numeric vector of Matsuda index values.
#'
#' @references
#' Matsuda M, DeFronzo RA (1999) *Diabetes Care* 22:1462-1470.
#'
#' @examples
#' matsuda(glucose_0 = 95, insulin_0 = 10,
#'         glucose_mean = 140, insulin_mean = 60)
#'
#' g <- c(95, 160, 150, 130, 110)
#' i <- c(10, 80, 90, 70, 45)
#' matsuda_from_curve(c(0, 30, 60, 90, 120), g, i)          # full
#' matsuda_from_curve(c(0, 30, 60, 90, 120), g, i, upto = 30) # modified 0-30
#'
#' @name matsuda
NULL

#' @rdname matsuda
#' @export
matsuda <- function(glucose_0, insulin_0, glucose_mean, insulin_mean) {
  10000 / sqrt(glucose_0 * insulin_0 * glucose_mean * insulin_mean)
}

#' @rdname matsuda
#' @export
matsuda_from_curve <- function(times, glucose, insulin, upto = NULL) {
  times <- as.numeric(times)
  if (is.null(dim(glucose))) glucose <- matrix(glucose, nrow = 1L)
  if (is.null(dim(insulin))) insulin <- matrix(insulin, nrow = 1L)
  glucose <- as.matrix(glucose)
  insulin <- as.matrix(insulin)
  keep <- if (is.null(upto)) rep(TRUE, length(times)) else times <= upto
  if (sum(keep) < 2L) {
    stop("Need at least two timepoints at or before `upto`.", call. = FALSE)
  }
  t_sub <- times[keep]
  span <- max(t_sub) - min(t_sub)
  g_mean <- auc_trapezoid(t_sub, glucose[, keep, drop = FALSE]) / span
  i_mean <- auc_trapezoid(t_sub, insulin[, keep, drop = FALSE]) / span
  matsuda(glucose[, 1L], insulin[, 1L], g_mean, i_mean)
}


#' Clamp-analogue OGTT sensitivity indices (Gutt, Cederholm)
#'
#' Gutt's ISI(0,120) and the Cederholm index estimate glucose uptake per unit
#' insulin from the fasting-to-2 h glucose fall, the glucose load, body weight,
#' and the mean glucose and insulin across the curve.
#'
#' \deqn{Gutt = \frac{75000 + (G_0 - G_{120})\,0.19\,wt}{120\,\bar{G}\,\log_{10}\bar{I}}}
#' \deqn{Cederholm = \frac{75000 + (G_0 - G_{120})\,1.15\,180\,0.19\,wt}{120\,\bar{G}\,\log_{10}\bar{I}}}
#'
#' `gutt_isi()` takes glucose in **mg/dL**; `cederholm()` takes glucose in
#' **mmol/L**. Both take insulin in uU/mL and weight in kg, and use a base-10
#' logarithm of mean insulin, matching the reference implementation.
#'
#' @param glucose_0,glucose_120 Fasting and 2 h glucose (mg/dL for `gutt_isi()`,
#'   mmol/L for `cederholm()`).
#' @param glucose_mean Mean glucose across the OGTT (same unit as `glucose_0`).
#' @param insulin_mean Mean insulin across the OGTT, uU/mL.
#' @param weight Body weight, kg.
#'
#' @return A numeric vector of sensitivity values.
#'
#' @references
#' Gutt M et al. (2000) *Diabetes Res Clin Pract* 47:177-184.
#' Cederholm J, Wibell L (1990) *Diabetes Res Clin Pract* 10:167-175.
#'
#' @examples
#' gutt_isi(glucose_0 = 95, glucose_120 = 150, glucose_mean = 140,
#'          insulin_mean = 60, weight = 85)
#'
#' @name gutt_cederholm
NULL

#' @rdname gutt_cederholm
#' @export
gutt_isi <- function(glucose_0, glucose_120, glucose_mean, insulin_mean,
                     weight) {
  (75000 + (glucose_0 - glucose_120) * 0.19 * weight) /
    (120 * glucose_mean * log10(insulin_mean))
}

#' @rdname gutt_cederholm
#' @export
cederholm <- function(glucose_0, glucose_120, glucose_mean, insulin_mean,
                      weight) {
  (75000 + (glucose_0 - glucose_120) * 1.15 * 180 * 0.19 * weight) /
    (120 * glucose_mean * log10(insulin_mean))
}


#' Stumvoll OGTT sensitivity indices
#'
#' Regression-based sensitivity and metabolic-clearance estimates from Stumvoll
#' et al. (2000). All expect **insulin in pmol/L** and **glucose in mmol/L**, as
#' in the original derivation.
#'
#' \tabular{ll}{
#'   `stumvoll_isi()` \tab \eqn{0.226 - 0.0032\,BMI - 0.0000645\,I_{120} - 0.00375\,G_{90}} \cr
#'   `stumvoll_isi_demo()` \tab \eqn{0.222 - 0.00333\,BMI - 0.0000779\,I_{120} - 0.000422\,age} \cr
#'   `stumvoll_isi_bmi_independent()` \tab \eqn{0.156 - 0.0000459\,I_{120} - 0.000321\,I_0 - 0.00541\,G_{120}} \cr
#'   `stumvoll_mcr()` \tab \eqn{18.8 - 0.271\,BMI - 0.0052\,I_{120} - 0.27\,G_{90}}
#' }
#'
#' @param bmi Body mass index, kg/m^2.
#' @param insulin_0,insulin_120 Fasting and 2 h insulin, pmol/L.
#' @param glucose_90,glucose_120 90 min and 2 h glucose, mmol/L.
#' @param age Age, years.
#'
#' @return A numeric vector (sensitivity index, or metabolic clearance rate for
#'   `stumvoll_mcr()`).
#'
#' @references
#' Stumvoll M et al. (2000) *Diabetes Care* 23:295-301.
#'
#' @examples
#' stumvoll_isi(bmi = 30, insulin_120 = 400, glucose_90 = 8.5)
#'
#' @name stumvoll_isi
NULL

#' @rdname stumvoll_isi
#' @export
stumvoll_isi <- function(bmi, insulin_120, glucose_90) {
  0.226 - 0.0032 * bmi - 0.0000645 * insulin_120 - 0.00375 * glucose_90
}

#' @rdname stumvoll_isi
#' @export
stumvoll_isi_demo <- function(bmi, insulin_120, age) {
  0.222 - 0.00333 * bmi - 0.0000779 * insulin_120 - 0.000422 * age
}

#' @rdname stumvoll_isi
#' @export
stumvoll_isi_bmi_independent <- function(insulin_0, insulin_120, glucose_120) {
  0.156 - 0.0000459 * insulin_120 - 0.000321 * insulin_0 - 0.00541 * glucose_120
}

#' @rdname stumvoll_isi
#' @export
stumvoll_mcr <- function(bmi, insulin_120, glucose_90) {
  18.8 - 0.271 * bmi - 0.0052 * insulin_120 - 0.27 * glucose_90
}


#' Avignon insulin-sensitivity indices
#'
#' Avignon's Si indices scale the reciprocal of the insulin-glucose product by a
#' fixed glucose distribution volume (150 mL per kg body weight). `avignon_sib()`
#' uses the fasting (basal) state, `avignon_si2h()` the 2 h state, and
#' `avignon_sim()` combines them.
#'
#' \deqn{S_{i,b} = 10^8 / (I_0 \times G_0 \times 150\,wt)}
#' \deqn{S_{i,2h} = 10^8 / (I_{120} \times G_{120} \times 150\,wt)}
#' \deqn{S_{i,M} = (0.137 \times S_{i,b} + S_{i,2h}) / 2}
#'
#' Insulin is in uU/mL, glucose in mmol/L, weight in kg.
#'
#' @param insulin_0,insulin_120 Fasting and 2 h insulin, uU/mL.
#' @param glucose_0,glucose_120 Fasting and 2 h glucose, mmol/L.
#' @param weight Body weight, kg.
#' @param sib,si2h Basal and 2 h Avignon indices (as returned by
#'   `avignon_sib()` / `avignon_si2h()`), for `avignon_sim()`.
#'
#' @return A numeric vector of Avignon sensitivity values.
#'
#' @references
#' Avignon A et al. (1999) *Int J Obes Relat Metab Disord* 23:512-517.
#'
#' @examples
#' sib  <- avignon_sib(insulin_0 = 10, glucose_0 = 5.3, weight = 85)
#' si2h <- avignon_si2h(insulin_120 = 60, glucose_120 = 7.8, weight = 85)
#' avignon_sim(sib, si2h)
#'
#' @name avignon
NULL

#' @rdname avignon
#' @export
avignon_sib <- function(insulin_0, glucose_0, weight) {
  1e8 / (insulin_0 * glucose_0 * (150 * weight))
}

#' @rdname avignon
#' @export
avignon_si2h <- function(insulin_120, glucose_120, weight) {
  1e8 / (insulin_120 * glucose_120 * (150 * weight))
}

#' @rdname avignon
#' @export
avignon_sim <- function(sib, si2h) {
  ((0.137 * sib) + si2h) / 2
}


#' Belfiore OGTT insulin-sensitivity index
#'
#' Belfiore's index normalises the insulin-glucose product against reference
#' (normal) values and maps it to a 0-2 sensitivity scale. The OGTT form uses
#' the 0-1-2 h glucose and insulin areas; `belfiore_area()` computes that area
#' from the 0, 60 and 120 min values as \eqn{0.5 v_0 + v_{60} + 0.5 v_{120}}
#' (a trapezoidal 0-1-2 h area in mmol/L*h or uU/mL*h). See [belfiore_basal()]
#' for the fasting-only form.
#'
#' \deqn{2 / ((A_I / A_{I,ref})(A_G / A_{G,ref}) + 1)}
#'
#' @param v0,v60,v120 Analyte values at 0, 60 and 120 min (glucose mmol/L or
#'   insulin uU/mL) for `belfiore_area()`.
#' @param insulin_area,glucose_area 0-1-2 h insulin (uU/mL*h) and glucose
#'   (mmol/L*h) areas, as from `belfiore_area()`.
#' @param ref_insulin_area,ref_glucose_area Reference (normal) areas. Defaults
#'   are the Belfiore (1998) normal means, 91.87 uU/mL*h and 11.36 mmol/L*h.
#'
#' @return A numeric vector of Belfiore OGTT sensitivity values.
#'
#' @references
#' Belfiore F et al. (1998) *Mol Genet Metab* 63:134-141.
#'
#' @examples
#' gA <- belfiore_area(5.1, 9.0, 7.0)     # glucose area, mmol/L*h
#' iA <- belfiore_area(10, 80, 45)        # insulin area, uU/mL*h
#' belfiore_ogtt(insulin_area = iA, glucose_area = gA)
#'
#' @name belfiore_ogtt
NULL

#' @rdname belfiore_ogtt
#' @export
belfiore_area <- function(v0, v60, v120) {
  0.5 * v0 + v60 + 0.5 * v120
}

#' @rdname belfiore_ogtt
#' @export
belfiore_ogtt <- function(insulin_area, glucose_area,
                          ref_insulin_area = 91.87,
                          ref_glucose_area = 11.36) {
  2 / (((insulin_area / ref_insulin_area) *
          (glucose_area / ref_glucose_area)) + 1)
}
