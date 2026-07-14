#' Beta-cell secretion / insulinogenic indices
#'
#' Surrogates of insulin secretion and beta-cell function from fasting bloods
#' and the OGTT.
#'
#' \tabular{ll}{
#'   `homa_b()`   \tab \eqn{20 I_0 / (G_0 - 3.5)} (I uU/mL, G mmol/L) \cr
#'   `igi()`      \tab \eqn{(I_t - I_0) / (G_t - G_0)} (I uU/mL, G mmol/L) \cr
#'   `cir()`      \tab \eqn{I_t / (G_t (G_t - 70))} if \eqn{G_t > 70} (G mg/dL) \cr
#'   `kadowaki()` \tab \eqn{(I_{30} - I_0) / G_{30}} (I uU/mL, G mg/dL)
#' }
#'
#' `igi()` (the insulinogenic index) and `cir()` (corrected insulin response)
#' are timepoint-agnostic: pass the post-load insulin/glucose at whichever time
#' you want (30 min is the classic insulinogenic index). `cir()` returns `NA`
#' when \eqn{G_t \le 70} mg/dL to avoid a non-physiological denominator.
#'
#' @param insulin,insulin_0,insulin_t,insulin_30 Insulin, uU/mL. `insulin_0` is
#'   fasting; `insulin_t` is the post-load value at the chosen time.
#' @param glucose Fasting glucose, mmol/L (for `homa_b()`).
#' @param glucose_0 Fasting glucose, mmol/L (for `igi()`).
#' @param glucose_t Post-load glucose: mmol/L for `igi()`, mg/dL for `cir()`.
#' @param glucose_30 30 min glucose, mg/dL (for `kadowaki()`).
#'
#' @return A numeric vector the length of the recycled inputs.
#'
#' @references
#' Matthews DR et al. (1985) HOMA. *Diabetologia* 28:412-419.
#' Seltzer HS et al. (1967) insulinogenic index. *J Clin Invest* 46:323-335.
#' Sluiter WJ et al. (1976) corrected insulin response. *Diabetes* 25:245-249.
#'
#' @examples
#' homa_b(insulin = 10, glucose = 5.3)
#' igi(insulin_0 = 10, insulin_t = 80, glucose_0 = 5.3, glucose_t = 8.9)
#' cir(insulin_t = 80, glucose_t = 160)
#'
#' @name secretion
NULL

#' @rdname secretion
#' @export
homa_b <- function(insulin, glucose) {
  (20 * insulin) / (glucose - 3.5)
}

#' @rdname secretion
#' @export
igi <- function(insulin_0, insulin_t, glucose_0, glucose_t) {
  (insulin_t - insulin_0) / (glucose_t - glucose_0)
}

#' @rdname secretion
#' @export
cir <- function(insulin_t, glucose_t) {
  ifelse(glucose_t > 70,
         insulin_t / (glucose_t * (glucose_t - 70)),
         NA_real_)
}

#' @rdname secretion
#' @export
kadowaki <- function(insulin_0, insulin_30, glucose_30) {
  (insulin_30 - insulin_0) / glucose_30
}


#' Stumvoll first- and second-phase insulin secretion
#'
#' Regression estimates of first- and second-phase insulin secretion from the
#' 0 and 30 min OGTT samples (Stumvoll et al. 2000). Both expect **insulin in
#' pmol/L** and **glucose in mmol/L**, matching the original derivation. (Note
#' the GRADE report's data dictionary labels these as uU/mL, but the analysis
#' code and the source paper use pmol/L; SurInsulot follows pmol/L.)
#'
#' \deqn{1st = 1283 + 1.829 I_{30} - 138.7 G_{30} + 3.772 I_0}
#' \deqn{2nd = 287 + 0.416 I_{30} - 26.07 G_{30} + 0.926 I_0}
#'
#' @param insulin_0,insulin_30 Fasting and 30 min insulin, pmol/L.
#' @param glucose_30 30 min glucose, mmol/L.
#'
#' @return A numeric vector of estimated secretion (pmol/L).
#'
#' @references
#' Stumvoll M et al. (2000) *Diabetes Care* 23:295-301.
#'
#' @examples
#' stumvoll_first_phase(insulin_0 = 70, insulin_30 = 500, glucose_30 = 9)
#' stumvoll_second_phase(insulin_0 = 70, insulin_30 = 500, glucose_30 = 9)
#'
#' @name stumvoll_secretion
NULL

#' @rdname stumvoll_secretion
#' @export
stumvoll_first_phase <- function(insulin_0, insulin_30, glucose_30) {
  1283 + (1.829 * insulin_30) - (138.7 * glucose_30) + (3.772 * insulin_0)
}

#' @rdname stumvoll_secretion
#' @export
stumvoll_second_phase <- function(insulin_0, insulin_30, glucose_30) {
  287 + (0.416 * insulin_30) - (26.07 * glucose_30) + (0.926 * insulin_0)
}


#' BIGTT acute insulin response
#'
#' The BIGTT-AIR indices estimate the acute (first-phase) insulin response from
#' reduced OGTT sampling, sex, and BMI (Hansen et al. 2007). Two forms use
#' different intermediate timepoints (30 or 60 min). **Insulin is in uU/mL** and
#' **glucose in mmol/L**.
#'
#' \deqn{AIR_{0,30,120} = \exp(8.20 + 0.00178 I_0 + 0.00168 I_{30} - 0.000383 I_{120} - 0.314 G_0 - 0.109 G_{30} + 0.0781 G_{120} + 0.180\,male + 0.032\,BMI)}
#' \deqn{AIR_{0,60,120} = \exp(8.19 + 0.00339 I_0 + 0.00152 I_{60} - 0.000959 I_{120} - 0.389 G_0 - 0.142 G_{60} + 0.164 G_{120} + 0.256\,male + 0.038\,BMI)}
#'
#' @param insulin_0,insulin_30,insulin_60,insulin_120 Insulin at 0/30/60/120
#'   min, uU/mL.
#' @param glucose_0,glucose_30,glucose_60,glucose_120 Glucose at 0/30/60/120
#'   min, mmol/L.
#' @param bmi Body mass index, kg/m^2.
#' @param male Sex indicator: `1`/`TRUE` for men, `0`/`FALSE` for women.
#'
#' @return A numeric vector of estimated acute insulin response.
#'
#' @references
#' Hansen T et al. (2007) BIGTT. *Diabetes Care* 30:257-262.
#'
#' @examples
#' bigtt_air_30_120(insulin_0 = 10, insulin_30 = 80, insulin_120 = 60,
#'                  glucose_0 = 5.3, glucose_30 = 9, glucose_120 = 7.5,
#'                  bmi = 30, male = 1)
#'
#' @name bigtt_air
NULL

#' @rdname bigtt_air
#' @export
bigtt_air_30_120 <- function(insulin_0, insulin_30, insulin_120,
                             glucose_0, glucose_30, glucose_120, bmi, male) {
  m <- .male_indicator(male)
  exp(8.20 + (0.00178 * insulin_0) + (0.00168 * insulin_30) -
        (0.000383 * insulin_120) - (0.314 * glucose_0) -
        (0.109 * glucose_30) + (0.0781 * glucose_120) +
        (0.180 * m) + (0.032 * bmi))
}

#' @rdname bigtt_air
#' @export
bigtt_air_60_120 <- function(insulin_0, insulin_60, insulin_120,
                             glucose_0, glucose_60, glucose_120, bmi, male) {
  m <- .male_indicator(male)
  exp(8.19 + (0.00339 * insulin_0) + (0.00152 * insulin_60) -
        (0.000959 * insulin_120) - (0.389 * glucose_0) -
        (0.142 * glucose_60) + (0.164 * glucose_120) +
        (0.256 * m) + (0.038 * bmi))
}
