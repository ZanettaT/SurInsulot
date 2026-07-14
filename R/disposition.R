#' Disposition index
#'
#' The (oral) disposition index expresses insulin secretion relative to the
#' prevailing insulin sensitivity as their product: a given amount of secretion
#' means healthier beta cells if it occurs in an insulin-sensitive person.
#' `disposition_index()` is simply
#' \deqn{DI = secretion \times sensitivity}
#' and every named variant in the literature is this product with a particular
#' pair of surrogates plugged in.
#'
#' \tabular{ll}{
#'   oral DI (oDI)      \tab `igi()`\eqn{_{30}} x `1 / fasting insulin` \cr
#'   oDI (Matsuda)      \tab `igi()`\eqn{_{30}} x `matsuda()` \cr
#'   oDI (QUICKI)       \tab `igi()`\eqn{_{30}} x `quicki()` \cr
#'   oDI (1/HOMA-IR)    \tab `igi()`\eqn{_{30}} x `1 / homa_ir()` \cr
#'   oDI (C-peptide)    \tab `cpi()`\eqn{_{30}} x `matsuda()` \cr
#'   ISSI-2             \tab (AUC insulin:glucose, 0-120) x `matsuda()`
#' }
#'
#' Pass a secretion index (for example the insulinogenic index [igi()]) and a
#' *sensitivity* index (higher = more sensitive, for example [matsuda()],
#' [quicki()], or `1 / homa_ir()`). Do not pass a resistance index such as
#' HOMA-IR directly; invert it first.
#'
#' @param secretion A secretion surrogate (for example [igi()] or [cpi()]).
#' @param sensitivity A sensitivity surrogate where higher means more sensitive
#'   (for example [matsuda()], [quicki()], or `1 / homa_ir()`).
#'
#' @return A numeric vector of disposition-index values.
#'
#' @references
#' Utzschneider KM et al. (2009) oral disposition index.
#' *Diabetes Care* 32:335-341.
#' Retnakaran R et al. (2009) ISSI-2. *Obesity* 17:1232-1238.
#'
#' @examples
#' ig  <- igi(insulin_0 = 10, insulin_t = 80, glucose_0 = 5.3, glucose_t = 8.9)
#' mat <- matsuda(glucose_0 = 95, insulin_0 = 10,
#'                glucose_mean = 140, insulin_mean = 60)
#' disposition_index(secretion = ig, sensitivity = mat)   # oDI (Matsuda)
#'
#' @export
disposition_index <- function(secretion, sensitivity) {
  secretion * sensitivity
}
