#' Lipid- and anthropometry-based insulin-resistance indices
#'
#' Insulin-independent surrogates that estimate insulin resistance (or its
#' adiposity substrate) from lipids and body measurements. They are attractive
#' because they need no insulin assay.
#'
#' \tabular{ll}{
#'   `tg_hdl()`      \tab \eqn{TG / HDL} (mg/dL) \cr
#'   `tyg()`         \tab \eqn{\ln(TG \times G_0 / 2)} (mg/dL) \cr
#'   `tyg_adjusted()`\tab \eqn{TyG \times A} for an adiposity measure \eqn{A} \cr
#'   `mets_ir()`     \tab \eqn{\ln(2 G_0 + TG)\,BMI / \ln(HDL)} (mg/dL) \cr
#'   `mets_vf()`     \tab visceral-fat score from METS-IR, WHtR, sex, age \cr
#'   `spise()`       \tab \eqn{600\,HDL^{0.185} / (TG^{0.2} BMI^{1.338})} (mg/dL) \cr
#'   `vai()`         \tab sex-specific visceral adiposity index (mmol/L) \cr
#'   `lap()`         \tab sex-specific lipid accumulation product (mmol/L) \cr
#'   `leptin_adiponectin_ratio()` \tab \eqn{leptin / adiponectin}
#' }
#'
#' `tg_hdl()`, `tyg()`, `mets_ir()` and `spise()` use triglycerides, glucose and
#' HDL in **mg/dL**. `vai()` and `lap()` use triglycerides and HDL in **mmol/L**.
#' For `tyg_adjusted()`, multiply the TyG index by BMI (kg/m^2), waist (cm),
#' waist-to-height ratio, or waist-to-hip ratio to obtain TyG-BMI, TyG-WC,
#' TyG-WHtR or TyG-WHR respectively. Except for `spise()` and
#' `leptin_adiponectin_ratio()` (which rise with sensitivity), higher values
#' indicate more resistance / adiposity.
#'
#' @param trig Triglycerides. mg/dL for `tg_hdl()`, `tyg()`; mmol/L for `vai()`,
#'   `lap()`.
#' @param hdl HDL cholesterol. mg/dL for `tg_hdl()`, `mets_ir()`, `spise()`;
#'   mmol/L for `vai()`.
#' @param glucose Fasting glucose, mg/dL.
#' @param bmi Body mass index, kg/m^2.
#' @param waist Waist circumference, cm.
#' @param tyg_value A TyG index value, as returned by `tyg()`.
#' @param adiposity An adiposity measure to scale TyG by (BMI, waist, WHtR or
#'   WHR).
#' @param mets_ir_value A METS-IR value, as returned by `mets_ir()`.
#' @param whtr Waist-to-height ratio (unitless).
#' @param age Age, years.
#' @param male Sex indicator: `1`/`TRUE` for men, `0`/`FALSE` for women. Pass an
#'   explicit indicator (e.g. `sex == 1`), not a raw study sex code.
#' @param leptin Leptin, ng/mL.
#' @param adiponectin Adiponectin, ug/mL.
#'
#' @return A numeric vector the length of the recycled inputs.
#'
#' @references
#' Simental-Mendia LE et al. (2008) TyG. *Metab Syndr Relat Disord* 6:299-304.
#' Bello-Chavolla OY et al. (2018) METS-IR and METS-VF.
#' *Eur J Endocrinol* 178:533-544 / *Clin Nutr* 39:1613-1621.
#' Paulmichl K et al. (2016) SPISE. *Clin Chem* 62:1211-1219.
#' Amato MC et al. (2010) VAI. *Diabetes Care* 33:920-922.
#' Kahn HS (2005) LAP. *BMC Cardiovasc Disord* 5:26.
#'
#' @examples
#' tyg(trig = 150, glucose = 95)
#' spise(hdl = 45, trig = 150, bmi = 30)
#' vai(waist = 100, bmi = 30, trig = 1.7, hdl = 1.1, male = 1)
#'
#' @name lipid_anthro
NULL

#' @rdname lipid_anthro
#' @export
tg_hdl <- function(trig, hdl) {
  trig / hdl
}

#' @rdname lipid_anthro
#' @export
tyg <- function(trig, glucose) {
  log(glucose * trig / 2)
}

#' @rdname lipid_anthro
#' @export
tyg_adjusted <- function(tyg_value, adiposity) {
  tyg_value * adiposity
}

#' @rdname lipid_anthro
#' @export
mets_ir <- function(glucose, trig, hdl, bmi) {
  (log((2 * glucose) + trig) * bmi) / log(hdl)
}

#' @rdname lipid_anthro
#' @export
mets_vf <- function(mets_ir_value, whtr, male, age) {
  m <- .male_indicator(male)
  4.466 + 0.011 * (log(mets_ir_value))^3 + 3.239 * (log(whtr))^3 +
    0.319 * m + 0.594 * log(age)
}

#' @rdname lipid_anthro
#' @export
spise <- function(hdl, trig, bmi) {
  (600 * hdl^0.185) / (trig^0.2 * bmi^1.338)
}

#' @rdname lipid_anthro
#' @export
vai <- function(waist, bmi, trig, hdl, male) {
  m <- .male_indicator(male)
  ifelse(
    m == 1,
    (waist / (39.68 + 1.88 * bmi)) * (trig / 1.03) * (1.31 / hdl),
    (waist / (36.58 + 1.89 * bmi)) * (trig / 0.81) * (1.52 / hdl)
  )
}

#' @rdname lipid_anthro
#' @export
lap <- function(waist, trig, male) {
  m <- .male_indicator(male)
  ifelse(m == 1, (waist - 65) * trig, (waist - 58) * trig)
}

#' @rdname lipid_anthro
#' @export
leptin_adiponectin_ratio <- function(leptin, adiponectin) {
  leptin / adiponectin
}
