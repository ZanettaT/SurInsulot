#' Fasting insulin-sensitivity / resistance indices
#'
#' Surrogate indices of insulin sensitivity or resistance computed from a single
#' fasting blood draw. Each function is a direct transcription of the published
#' formula and is vectorised over its inputs.
#'
#' \tabular{ll}{
#'   `raynaud()`                 \tab \eqn{40 / I_0} \cr
#'   `homa_ir()`                 \tab \eqn{(I_0 \times G_0) / 22.5} \cr
#'   `firi()`                    \tab \eqn{(I_0 \times G_0) / 25} \cr
#'   `quicki()`                  \tab \eqn{1 / (\log_{10} I_0 + \log_{10} G_0)} \cr
#'   `insulin_glucose_ratio()`   \tab \eqn{I_0 / G_0} \cr
#'   `glucose_insulin_ratio()`   \tab \eqn{G_0 / I_0} \cr
#'   `bennett_isi()`             \tab \eqn{1 / (\ln I_0 \times \ln G_0)} \cr
#'   `mcauley()`                 \tab \eqn{\exp(2.63 - 0.28 \ln I_0 - 0.31 \ln TG)} \cr
#'   `ohkura()`                  \tab \eqn{20 / (CP_0 \times G_0)} \cr
#'   `homa_ad()`                 \tab \eqn{(I_0 \times G_0) / (22.5 \times ADIPO)} \cr
#'   `belfiore_basal()`          \tab \eqn{2 / ((I_0/I_{ref})(G_0/G_{ref}) + 1)}
#' }
#'
#' `homa_ir()`, `firi()`, `insulin_glucose_ratio()`, `glucose_insulin_ratio()`,
#' `bennett_isi()`, `homa_ad()` and `belfiore_basal()` take fasting glucose in
#' **mmol/L**; `quicki()` takes fasting glucose in **mg/dL**. `mcauley()` takes
#' triglycerides in **mmol/L**. `ohkura()` takes fasting C-peptide in nmol/L and
#' glucose in mmol/L. Insulin is in uU/mL (= mU/L) throughout, adiponectin in
#' ug/mL. Higher `homa_ir()`, `firi()`, `insulin_glucose_ratio()`, `homa_ad()`
#' mean more resistance; the others increase with sensitivity. `1 / homa_ir()`
#' gives the reciprocal-HOMA sensitivity index.
#'
#' @param insulin Fasting insulin, uU/mL (= mU/L).
#' @param glucose Fasting glucose. mmol/L for all functions except `quicki()`,
#'   which expects mg/dL.
#' @param trig Fasting triglycerides, mmol/L.
#' @param cpeptide Fasting C-peptide, nmol/L.
#' @param adiponectin Fasting adiponectin, ug/mL.
#' @param ref_insulin,ref_glucose Reference (normal) fasting insulin (uU/mL) and
#'   glucose (mmol/L) for the Belfiore normalisation. Defaults are the
#'   Belfiore (1998) normal means, 9.46 uU/mL and 5.08 mmol/L.
#'
#' @return A numeric vector the length of the recycled inputs.
#'
#' @references
#' Matthews DR et al. (1985) Homeostasis model assessment. *Diabetologia*
#' 28:412-419.
#' Katz A et al. (2000) QUICKI. *J Clin Endocrinol Metab* 85:2402-2410.
#' Raynaud E et al. (1999) *Diabetes Care* 22:1003-1004.
#' McAuley KA et al. (2001) *Diabetes Care* 24:460-464.
#' Ohkura T et al. (2013) *Diabetol Metab Syndr* 5:43.
#' Belfiore F et al. (1998) *Mol Genet Metab* 63:134-141.
#'
#' @examples
#' homa_ir(insulin = 10, glucose = 5.5)
#' quicki(insulin = 10, glucose = 99)          # glucose in mg/dL
#' mcauley(insulin = 10, trig = 1.3)
#'
#' @name fasting_sensitivity
NULL

#' @rdname fasting_sensitivity
#' @export
raynaud <- function(insulin) {
  40 / insulin
}

#' @rdname fasting_sensitivity
#' @export
homa_ir <- function(insulin, glucose) {
  (insulin * glucose) / 22.5
}

#' @rdname fasting_sensitivity
#' @export
firi <- function(insulin, glucose) {
  (insulin * glucose) / 25
}

#' @rdname fasting_sensitivity
#' @export
quicki <- function(insulin, glucose) {
  1 / (log10(insulin) + log10(glucose))
}

#' @rdname fasting_sensitivity
#' @export
insulin_glucose_ratio <- function(insulin, glucose) {
  insulin / glucose
}

#' @rdname fasting_sensitivity
#' @export
glucose_insulin_ratio <- function(insulin, glucose) {
  glucose / insulin
}

#' @rdname fasting_sensitivity
#' @export
bennett_isi <- function(insulin, glucose) {
  1 / (log(insulin) * log(glucose))
}

#' @rdname fasting_sensitivity
#' @export
mcauley <- function(insulin, trig) {
  exp(2.63 - 0.28 * log(insulin) - 0.31 * log(trig))
}

#' @rdname fasting_sensitivity
#' @export
ohkura <- function(cpeptide, glucose) {
  20 / (cpeptide * glucose)
}

#' @rdname fasting_sensitivity
#' @export
homa_ad <- function(insulin, glucose, adiponectin) {
  (insulin * glucose) / (22.5 * adiponectin)
}

#' @rdname fasting_sensitivity
#' @export
belfiore_basal <- function(insulin, glucose, ref_insulin = 9.46,
                           ref_glucose = 5.08) {
  2 / (((insulin / ref_insulin) * (glucose / ref_glucose)) + 1)
}


#' Estimated glucose disposal rate (eGDR)
#'
#' The eGDR estimates insulin-stimulated glucose disposal (mg/kg/min) from
#' routine clinical variables; lower values indicate more insulin resistance.
#' Three published forms differ only in the adiposity term (BMI, waist
#' circumference, or waist-to-hip ratio).
#'
#' \tabular{ll}{
#'   `egdr_bmi()`   \tab \eqn{19.02 - 0.22\,BMI - 3.26\,HTN - 0.61\,HbA1c} \cr
#'   `egdr_waist()` \tab \eqn{21.158 - 0.09\,WC - 3.407\,HTN - 0.551\,HbA1c} \cr
#'   `egdr_whr()`   \tab \eqn{24.31 - 12.22\,WHR - 3.29\,HTN - 0.57\,HbA1c}
#' }
#'
#' @param hba1c HbA1c, percent (NGSP / DCCT units).
#' @param hypertension Hypertension indicator: `1` (or `TRUE`) if hypertensive,
#'   `0` (or `FALSE`) otherwise.
#' @param bmi Body mass index, kg/m^2.
#' @param waist Waist circumference, cm.
#' @param whr Waist-to-hip ratio (unitless).
#'
#' @return A numeric vector of eGDR values (mg/kg/min).
#'
#' @references
#' Williams KV et al. (2000) *Diabetes* 49:626-632 (BMI form).
#' Epanomeritakis / Garofolo et al. (waist and WHR forms) as used in the
#' EURODIAB-derived equations.
#'
#' @examples
#' egdr_bmi(hba1c = 7.2, hypertension = 0, bmi = 31)
#' egdr_waist(hba1c = 7.2, hypertension = 1, waist = 104)
#'
#' @name egdr
NULL

#' @rdname egdr
#' @export
egdr_bmi <- function(hba1c, hypertension, bmi) {
  19.02 - (0.22 * bmi) - (3.26 * hypertension) - (0.61 * hba1c)
}

#' @rdname egdr
#' @export
egdr_waist <- function(hba1c, hypertension, waist) {
  21.158 - (0.09 * waist) - (3.407 * hypertension) - (0.551 * hba1c)
}

#' @rdname egdr
#' @export
egdr_whr <- function(hba1c, hypertension, whr) {
  24.31 - (12.22 * whr) - (3.29 * hypertension) - (0.57 * hba1c)
}
