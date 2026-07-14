#' Unit conversions
#'
#' Convert between the mass and molar units used across the insulin surrogate
#' literature. Conversions are element-wise and preserve `NA`.
#'
#' The molar mass factors are 18.0 for glucose, 88.57 for triglycerides, and
#' 38.67 for cholesterol fractions (HDL, LDL, total). Insulin is converted with
#' the IUPAC / consensus factor 1 uU/mL (= 1 mU/L) = 6.945 pmol/L.
#'
#' @param x Numeric vector of measurements in the source unit.
#'
#' @return A numeric vector in the target unit.
#'
#' @examples
#' mgdl_to_mmol_glucose(90)      # 5.0 mmol/L
#' mmol_to_mgdl_glucose(5)       # 90 mg/dL
#' insulin_uU_to_pmol(10)        # 69.45 pmol/L
#' insulin_pmol_to_uU(69.45)     # 10 uU/mL
#'
#' @name unit_conversions
NULL

#' @rdname unit_conversions
#' @export
mgdl_to_mmol_glucose <- function(x) x / 18.0

#' @rdname unit_conversions
#' @export
mmol_to_mgdl_glucose <- function(x) x * 18.0

#' @rdname unit_conversions
#' @export
mgdl_to_mmol_trig <- function(x) x / 88.57

#' @rdname unit_conversions
#' @export
mmol_to_mgdl_trig <- function(x) x * 88.57

#' @rdname unit_conversions
#' @export
mgdl_to_mmol_chol <- function(x) x / 38.67

#' @rdname unit_conversions
#' @export
mmol_to_mgdl_chol <- function(x) x * 38.67

#' @rdname unit_conversions
#' @export
insulin_uU_to_pmol <- function(x) x * 6.945

#' @rdname unit_conversions
#' @export
insulin_pmol_to_uU <- function(x) x / 6.945
