#' SurInsulot: Surrogate Indices of Insulin Sensitivity and Secretion
#'
#' SurInsulot computes surrogate indices of insulin resistance / sensitivity
#' and insulin deficiency / secretion from fasting bloods and the oral glucose
#' tolerance test (OGTT). Each index has its own small, vectorised function
#' documented with its formula, the units it expects, and a primary reference,
#' so the package doubles as a reproducible reference for how each measure is
#' defined. [compute_surrogates()] derives every estimable index from a data
#' frame in one call.
#'
#' @section Unit conventions:
#' Unless a function's documentation says otherwise, inputs are expected in:
#' * insulin: uU/mL (equivalently mU/L),
#' * C-peptide: nmol/L,
#' * glucagon and proinsulin: pmol/L,
#' * glucose in the mmol/L-based indices: mmol/L, and in the mg/dL-based
#'   indices (QUICKI, TyG, METS-IR, SPISE, CIR, Kadowaki, Matsuda, Gutt):
#'   mg/dL,
#' * triglycerides and HDL cholesterol: mg/dL unless an index is defined in
#'   mmol/L (VAI, LAP, McAuley),
#' * weight: kg, height: cm, waist/hip: cm, BMI: kg/m^2, age: years.
#'
#' Helpers [mgdl_to_mmol_glucose()], [insulin_uU_to_pmol()] and friends convert
#' between the two systems. A handful of indices were originally derived in
#' pmol/L (the Stumvoll indices); their functions expect insulin in pmol/L and
#' say so explicitly.
#'
#' @section Index families:
#' * Fasting sensitivity: [homa_ir()], [quicki()], [raynaud()], [firi()],
#'   [insulin_glucose_ratio()], [glucose_insulin_ratio()], [bennett_isi()],
#'   [mcauley()], [ohkura()], [homa_ad()], [belfiore_basal()], [egdr_bmi()].
#' * HOMA2 (Oxford DTU model, built in): [homa2_insulin()],
#'   [homa2_specific_insulin()], [homa2_cpeptide()].
#' * Lipid / anthropometric: [tg_hdl()], [tyg()], [tyg_adjusted()],
#'   [mets_ir()], [mets_vf()], [spise()], [vai()], [lap()],
#'   [leptin_adiponectin_ratio()].
#' * OGTT sensitivity: [matsuda()], [gutt_isi()], [stumvoll_isi()],
#'   [cederholm()], [avignon_sib()], [ogis()], [belfiore_ogtt()].
#' * Secretion: [homa_b()], [igi()], [cir()], [kadowaki()],
#'   [stumvoll_first_phase()], [stumvoll_second_phase()], [bigtt_air()].
#' * C-peptide: [cpeptide_glucose_ratio()], [cpi()],
#'   [insulin_cpeptide_ratio()], [cpeptide_insulin_auc_ratio()].
#' * Disposition: [disposition_index()].
#' * Glucagon: [glucagon_insulin_ratio()], [glucagon_suppression()].
#' * Helpers: [auc_trapezoid()], [iauc_trapezoid()], [guard_ratio()],
#'   [row_mean()].
#'
#' @keywords internal
"_PACKAGE"
