#' Glucagon indices
#'
#' Simple fasting and OGTT glucagon measures. `glucagon_insulin_ratio()` is the
#' molar fasting glucagon:insulin ratio (both in pmol/L). `glucagon_suppression()`
#' quantifies the fall in glucagon from fasting to 30 min, either as an absolute
#' change or as a percentage of the fasting value.
#'
#' \deqn{glucagon:insulin = GCG_0 / I_0}
#' \deqn{suppression_{abs} = GCG_0 - GCG_{30}}
#' \deqn{suppression_{pct} = 100 (GCG_0 - GCG_{30}) / GCG_0}
#'
#' @param glucagon_0,glucagon_30 Fasting and 30 min glucagon, pmol/L.
#' @param insulin_0 Fasting insulin, pmol/L.
#' @param percent If `TRUE`, `glucagon_suppression()` returns the percentage
#'   suppression; if `FALSE` (default), the absolute change.
#'
#' @return A numeric vector.
#'
#' @examples
#' glucagon_insulin_ratio(glucagon_0 = 15, insulin_0 = 70)
#' glucagon_suppression(glucagon_0 = 15, glucagon_30 = 9)              # absolute
#' glucagon_suppression(glucagon_0 = 15, glucagon_30 = 9, percent = TRUE)
#'
#' @name glucagon
NULL

#' @rdname glucagon
#' @export
glucagon_insulin_ratio <- function(glucagon_0, insulin_0) {
  glucagon_0 / insulin_0
}

#' @rdname glucagon
#' @export
glucagon_suppression <- function(glucagon_0, glucagon_30, percent = FALSE) {
  if (percent) {
    ((glucagon_0 - glucagon_30) / glucagon_0) * 100
  } else {
    glucagon_0 - glucagon_30
  }
}
