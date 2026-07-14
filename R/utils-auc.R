#' Trapezoidal area under the curve
#'
#' `auc_trapezoid()` returns the total area under a concentration-time curve by
#' the trapezoidal rule. `iauc_trapezoid()` returns the *incremental* area,
#' subtracting the baseline (first timepoint) from every timepoint before
#' integrating, so below-baseline excursions reduce the total (net iAUC).
#'
#' Both functions are vectorised over subjects: pass a numeric matrix with one
#' row per subject and one column per timepoint, and a length-*k* `times`
#' vector giving the sampling times. A single numeric vector is treated as one
#' subject's curve. Any `NA` among the timepoints used for a subject propagates
#' to that subject's area.
#'
#' @param times Numeric vector of sampling times (for example
#'   `c(0, 30, 60, 90, 120)`), length equal to the number of columns of
#'   `values`. Spacing need not be regular.
#' @param values Numeric matrix (subjects in rows, timepoints in columns) or a
#'   single numeric vector for one curve.
#'
#' @return A numeric vector of areas, one per subject (length 1 for a single
#'   curve). Units are the product of the measurement unit and the time unit
#'   (for example mg/dL x min).
#'
#' @examples
#' # one glucose curve at 0/30/60/90/120 min
#' g <- c(90, 150, 140, 120, 110)
#' auc_trapezoid(c(0, 30, 60, 90, 120), g)
#' iauc_trapezoid(c(0, 30, 60, 90, 120), g)
#'
#' # several subjects at once
#' m <- rbind(c(90, 150, 110), c(100, 180, 130))
#' auc_trapezoid(c(0, 30, 120), m)
#'
#' @seealso [guard_ratio()], [row_mean()]
#' @name auc
NULL

#' @rdname auc
#' @export
auc_trapezoid <- function(times, values) {
  times <- as.numeric(times)
  if (length(times) < 2L) {
    stop("`times` must have at least two timepoints.", call. = FALSE)
  }
  if (is.null(dim(values))) {
    values <- matrix(values, nrow = 1L)
  }
  values <- as.matrix(values)
  if (ncol(values) != length(times)) {
    stop(
      "`values` must have one column per element of `times` (",
      length(times), " expected, ", ncol(values), " supplied).",
      call. = FALSE
    )
  }
  w <- diff(times)
  seg <- vapply(
    seq_along(w),
    function(i) w[i] * (values[, i] + values[, i + 1L]) / 2,
    numeric(nrow(values))
  )
  if (is.null(dim(seg))) sum(seg) else rowSums(seg)
}

#' @rdname auc
#' @export
iauc_trapezoid <- function(times, values) {
  if (is.null(dim(values))) {
    values <- matrix(values, nrow = 1L)
  }
  values <- as.matrix(values)
  base <- values[, 1L]
  auc_trapezoid(times, values - base)
}

#' Guarded ratio
#'
#' Divide `num` by `den`, returning `NA` wherever the denominator is not
#' strictly greater than `min_den`. This protects ratio indices (for example
#' incremental AUC ratios) from dividing by zero or by a non-positive
#' increment.
#'
#' @param num,den Numeric vectors, the numerator and denominator.
#' @param min_den Numeric scalar; the ratio is computed only where
#'   `den > min_den`. Defaults to `0`.
#'
#' @return A numeric vector of the same length as `num` / `den`, with `NA`
#'   where the guard fails.
#'
#' @examples
#' guard_ratio(c(10, 20, 5), c(2, 0, -1))   # 5, NA, NA
#'
#' @export
guard_ratio <- function(num, den, min_den = 0) {
  ifelse(den > min_den, num / den, NA_real_)
}

#' Row-wise mean of repeated measures
#'
#' Average several vectors of repeated measurements element-wise, ignoring
#' `NA`. Rows for which every value is missing return `NA` (rather than `NaN`).
#'
#' @param ... Two or more numeric vectors of equal length (for example repeated
#'   blood-pressure or anthropometry readings).
#'
#' @return A numeric vector of row means.
#'
#' @examples
#' row_mean(c(140, NA, 130), c(138, NA, 132))
#'
#' @export
row_mean <- function(...) {
  m <- rowMeans(cbind(...), na.rm = TRUE)
  m[is.nan(m)] <- NA_real_
  m
}
