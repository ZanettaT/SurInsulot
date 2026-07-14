#' OGIS oral glucose insulin sensitivity (2-hour model)
#'
#' The Oral Glucose Insulin Sensitivity index (OGIS) estimates a clamp-like
#' glucose clearance from a 2 h OGTT using glucose at 0/90/120 min, insulin at
#' 0/90 min, body size, and the glucose dose. It is a closed-form
#' reimplementation of the model published by Mari et al. and distributed as a
#' spreadsheet at <http://webmet.pd.cnr.it/ogis/>.
#'
#' The calculation proceeds through body-surface area (BSA), the dose per unit
#' area \eqn{D_0}, a clearance term, and a quadratic solution:
#' \deqn{BSA = 0.1640443958298 \times wt^{0.515} \times (height/100)^{0.422}}
#' \deqn{D_0 = 5.551 \times dose / BSA}
#' \deqn{CL = p_4 \left(\frac{p_1 D_0 - V (G_{120}-G_{90})/\tau}{G_{90}} + \frac{p_3}{G_0}\right) / (I_{90} - I_0 + p_2)}
#' \deqn{B = (p_5 (G_{90} - G_{cl}) + 1)\,CL}
#' \deqn{OGIS = \tfrac{1}{2}\left(B + \sqrt{B^2 + 4 p_5 p_6 (G_{90}-G_{cl})\,CL}\right)}
#'
#' Glucose is in **mmol/L**, insulin in **pmol/L**, weight in kg, height in cm.
#' Note the model constant \eqn{p_2} is on the pmol/L insulin scale; convert
#' uU/mL insulin with [insulin_uU_to_pmol()] first.
#'
#' @param glucose_0,glucose_90,glucose_120 Glucose at 0, 90 and 120 min, mmol/L.
#' @param insulin_0,insulin_90 Insulin at 0 and 90 min, pmol/L.
#' @param weight Body weight, kg.
#' @param height Height, cm.
#' @param dose Oral glucose load, grams (default 75).
#' @param p1,p2,p3,p4,p5,p6 Model constants (defaults from the published
#'   OGIS 2 h model: 6.5, 1951, 4514, 792, 0.0118, 173).
#' @param volume Glucose distribution volume, mL/m^2 (default 10000).
#' @param clamp_glucose Reference clamp glucose, mmol/L (default 5).
#' @param interval Interval between the 90 and 120 min samples, minutes
#'   (default 30).
#'
#' @return A numeric vector of OGIS values (mL/min/m^2).
#'
#' @references
#' Mari A et al. (2001) A model-based method for assessing insulin sensitivity
#' from the oral glucose tolerance test. *Diabetes Care* 24:539-548.
#'
#' @examples
#' ogis(glucose_0 = 5.0, glucose_90 = 8.0, glucose_120 = 7.0,
#'      insulin_0 = 60, insulin_90 = 400, weight = 85, height = 170)
#'
#' @export
ogis <- function(glucose_0, glucose_90, glucose_120, insulin_0, insulin_90,
                 weight, height, dose = 75,
                 p1 = 6.5, p2 = 1951, p3 = 4514, p4 = 792,
                 p5 = 0.0118, p6 = 173,
                 volume = 10000, clamp_glucose = 5, interval = 30) {
  bsa <- 0.1640443958298 * weight^0.515 * (0.01 * height)^0.422
  d0 <- 5.551 * dose / bsa
  cl <- p4 * ((p1 * d0 - volume * (glucose_120 - glucose_90) / interval) /
                glucose_90 + p3 / glucose_0) /
    (insulin_90 - insulin_0 + p2)
  b <- (p5 * (glucose_90 - clamp_glucose) + 1) * cl
  (b + sqrt(b^2 + 4 * p5 * p6 * (glucose_90 - clamp_glucose) * cl)) / 2
}
