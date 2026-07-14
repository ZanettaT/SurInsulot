#' Yeo-Johnson power transform
#'
#' The Yeo-Johnson transform normalises a numeric vector while remaining defined
#' for zero and negative values (unlike a log or Box-Cox). `yeo_johnson()`
#' applies the transform for a given power parameter `lambda`;
#' `yeo_johnson_lambda()` estimates `lambda` by maximum likelihood - the value
#' that makes the transformed vector most nearly Gaussian. This is the transform
#' behind the `*_yj` columns in the GRADE/DPP reports and the `transform = TRUE`
#' option of [compute_surrogates()].
#'
#' For \eqn{x \ge 0}: \eqn{((x+1)^\lambda - 1)/\lambda}, or \eqn{\ln(x+1)} when
#' \eqn{\lambda = 0}. For \eqn{x < 0}:
#' \eqn{-((-x+1)^{2-\lambda} - 1)/(2-\lambda)}, or \eqn{-\ln(-x+1)} when
#' \eqn{\lambda = 2}. \eqn{\lambda = 1} is the identity. Non-finite values pass
#' through unchanged. `lambda` is fitted over all finite values; refit it on each
#' training fold under cross-validation to avoid leakage.
#'
#' @param x A numeric vector.
#' @param lambda The power parameter (for example as returned by
#'   `yeo_johnson_lambda()`).
#' @param interval Search interval for `lambda` (default `c(-5, 5)`).
#'
#' @return `yeo_johnson()` returns the transformed vector, the same length as
#'   `x`. `yeo_johnson_lambda()` returns the single fitted `lambda`, or `NA` when
#'   `x` has fewer than 10 finite values or zero variance.
#'
#' @references
#' Yeo IK, Johnson RA (2000) A new family of power transformations to improve
#' normality or symmetry. *Biometrika* 87:954-959.
#'
#' @examples
#' x <- c(1, 2, 2, 3, 5, 8, 13, 21, 34, 55)   # right-skewed
#' lam <- yeo_johnson_lambda(x)
#' yeo_johnson(x, lam)
#' yeo_johnson(x, 1)                           # lambda = 1 is the identity
#'
#' @name yeo_johnson
NULL

#' @rdname yeo_johnson
#' @export
yeo_johnson <- function(x, lambda) {
  eps <- 1e-8
  out <- x
  pos <- is.finite(x) & x >= 0
  neg <- is.finite(x) & x <  0
  if (abs(lambda) < eps) {
    out[pos] <- log1p(x[pos])
  } else {
    out[pos] <- ((x[pos] + 1)^lambda - 1) / lambda
  }
  if (abs(lambda - 2) < eps) {
    out[neg] <- -log1p(-x[neg])
  } else {
    out[neg] <- -(((-x[neg] + 1)^(2 - lambda) - 1) / (2 - lambda))
  }
  out
}

# Profile log-likelihood over lambda (normal likelihood + transform Jacobian).
.yj_loglik <- function(lambda, x) {
  z <- yeo_johnson(x, lambda)
  s2 <- mean((z - mean(z))^2)
  if (!is.finite(s2) || s2 <= 0) return(-Inf)
  -0.5 * length(x) * log(s2) + (lambda - 1) * sum(sign(x) * log1p(abs(x)))
}

#' @rdname yeo_johnson
#' @export
yeo_johnson_lambda <- function(x, interval = c(-5, 5)) {
  x <- x[is.finite(x)]
  if (length(x) < 10 || stats::sd(x) == 0) {
    return(NA_real_)
  }
  stats::optimize(.yj_loglik, interval = interval, x = x, maximum = TRUE)$maximum
}
