#' Compute every estimable surrogate index from a data frame
#'
#' `compute_surrogates()` derives the full panel of insulin sensitivity /
#' resistance and secretion surrogates from a data frame of fasting and OGTT
#' measurements, reproducing the GRADE and DPP analysis pipelines. It computes
#' each index only when its required inputs are present as columns, so the same
#' call works on a rich OGTT dataset or on a fasting-only dataset - indices that
#' cannot be computed are simply skipped.
#'
#' @section Expected columns:
#' Supply any subset of the following, using these names and units. Missing
#' columns just switch off the indices that need them.
#' \tabular{lll}{
#'   **Column** \tab **Meaning** \tab **Unit** \cr
#'   `g000,g015,g030,g060,g090,g120` \tab glucose at 0-120 min \tab mg/dL \cr
#'   `i000,i015,i030,i060,i090,i120` \tab insulin at 0-120 min \tab uU/mL \cr
#'   `cp000,cp015,cp030,cp060,cp090,cp120` \tab C-peptide at 0-120 min \tab nmol/L \cr
#'   `glcg000,glcg030` \tab glucagon at 0 / 30 min \tab pmol/L \cr
#'   `trig,chdl` \tab triglycerides, HDL cholesterol \tab mg/dL \cr
#'   `hba1c_perc` \tab HbA1c \tab percent \cr
#'   `adipon,leptin,pin` \tab adiponectin, leptin, fasting proinsulin \tab ug/mL, ng/mL, pmol/L \cr
#'   `bmi,weight_kg,waist_cm,hip_cm,height_cm,whr,whtr` \tab anthropometry \tab kg/m^2, kg, cm, cm, cm, ratio, ratio \cr
#'   `sex,age,hypertension` \tab covariates \tab see `sex_col`/`male_code`; years; 1/0 \cr
#' }
#' The units in the table are the defaults. If your data are in other units,
#' declare them with the `*_unit` arguments and they are converted internally
#' first - e.g. `glucose_unit = "mmol/L"`, `insulin_unit = "pmol/L"`,
#' `cpeptide_unit = "ng/mL"`. Everything is normalised to the canonical units
#' above before any index is computed, and glucose/insulin are then converted to
#' mmol/L / pmol/L for the indices that need them.
#'
#' @param data A data frame with one row per observation (for example one
#'   participant-visit) and the columns above.
#' @param glucose_unit Unit of the `g*` columns: `"mg/dL"` (default) or `"mmol/L"`.
#' @param insulin_unit Unit of the `i*` columns: `"uU/mL"` (default) or `"pmol/L"`.
#' @param cpeptide_unit Unit of the `cp*` columns: `"nmol/L"` (default),
#'   `"ng/mL"` (divided by 3.02), or `"pmol/L"` (divided by 1000).
#' @param trig_unit Unit of `trig`: `"mg/dL"` (default) or `"mmol/L"`.
#' @param hdl_unit Unit of `chdl`: `"mg/dL"` (default) or `"mmol/L"`.
#' @param glucagon_unit Unit of `glcg*`: `"pmol/L"` (default) or `"pg/mL"`
#'   (divided by 3.483, glucagon molar mass ~3483 g/mol).
#' @param sex_col Name of the column holding sex (default `"sex"`). Used only by
#'   the sex-specific indices (VAI, LAP, METS-VF, BIGTT).
#' @param male_code The value of `sex_col` that denotes men (default `1`). Men
#'   are encoded as 1 and everyone else as 0 for the sex-specific indices.
#' @param verbose If `TRUE` (default), message how many index columns were
#'   computed.
#' @param transform If `TRUE`, also append a Yeo-Johnson-normalised
#'   `<index>_yj` column for every computed index (lambda fitted per index by
#'   maximum likelihood over the supplied rows), reproducing the `*_yj` outputs
#'   of the GRADE/DPP reports. The fitted lambdas are stored in
#'   `attr(result, "surinsulot_yj_lambda")`. Default `FALSE`.
#'
#' @return `data` with one appended column per computed index (names match the
#'   GRADE/DPP scripts, for example `homa_ir`, `matsuda`, `igi_30`, `ogis`), and
#'   an extra `<index>_yj` column per index when `transform = TRUE`. The names of
#'   the computed indices are stored in `attr(result, "surinsulot_computed")`.
#'
#' @details
#' HOMA2 indices (`homa2_ir`, `homa2_b_ins`, `homa2_s`, and the C-peptide
#' variants `homa2_ir_cp`, `homa2_b_cp`, `homa2_s_cp`) are computed with the
#' built-in [homa2_insulin()] and [homa2_cpeptide()] functions (Oxford DTU
#' reference tables); fasting insulin is converted to pmol/L at 6.945 before the
#' table lookup.
#'
#' @examples
#' df <- data.frame(
#'   g000 = c(95, 110), g030 = c(160, 190), g120 = c(140, 170),
#'   i000 = c(10, 18),  i030 = c(80, 60),
#'   trig = c(150, 200), chdl = c(45, 38), bmi = c(29, 34),
#'   sex = c(1, 2)
#' )
#' out <- compute_surrogates(df)
#' out[, c("homa_ir", "quicki", "tyg", "igi_30")]
#'
#' @seealso The individual index functions, e.g. [homa_ir()], [matsuda()],
#'   [igi()], [ogis()].
#' @export
compute_surrogates <- function(data,
                               glucose_unit = c("mg/dL", "mmol/L"),
                               insulin_unit = c("uU/mL", "pmol/L"),
                               cpeptide_unit = c("nmol/L", "ng/mL", "pmol/L"),
                               trig_unit = c("mg/dL", "mmol/L"),
                               hdl_unit = c("mg/dL", "mmol/L"),
                               glucagon_unit = c("pmol/L", "pg/mL"),
                               sex_col = "sex", male_code = 1,
                               verbose = TRUE, transform = FALSE) {
  stopifnot(is.data.frame(data))
  glucose_unit  <- match.arg(glucose_unit)
  insulin_unit  <- match.arg(insulin_unit)
  cpeptide_unit <- match.arg(cpeptide_unit)
  trig_unit     <- match.arg(trig_unit)
  hdl_unit      <- match.arg(hdl_unit)
  glucagon_unit <- match.arg(glucagon_unit)
  get_col <- function(nm) if (nm %in% names(data)) data[[nm]] else NULL

  ## ---- raw inputs -----------------------------------------------------------
  trig <- get_col("trig"); chdl <- get_col("chdl")
  bmi <- get_col("bmi"); weight <- get_col("weight_kg")
  waist <- get_col("waist_cm"); hip <- get_col("hip_cm")
  height <- get_col("height_cm"); whr <- get_col("whr"); whtr <- get_col("whtr")
  hba1c <- get_col("hba1c_perc"); adipon <- get_col("adipon")
  leptin <- get_col("leptin"); pin <- get_col("pin")
  age <- get_col("age"); htn <- get_col("hypertension")
  glcg000 <- get_col("glcg000"); glcg030 <- get_col("glcg030")

  sexv <- get_col(sex_col)
  male <- if (!is.null(sexv)) as.numeric(sexv == male_code) else NULL

  ## timepoint-indexed lists (NULL where absent)
  tp <- c("0", "15", "30", "60", "90", "120")
  suffix <- c("000", "015", "030", "060", "090", "120")
  G <- lapply(paste0("g", suffix), get_col);  names(G) <- tp
  I <- lapply(paste0("i", suffix), get_col);  names(I) <- tp
  CP <- lapply(paste0("cp", suffix), get_col); names(CP) <- tp

  ## ---- normalise raw inputs to the package's canonical units ---------------
  ## canonical: glucose mg/dL, insulin uU/mL, C-peptide nmol/L, TG/HDL mg/dL,
  ## glucagon pmol/L. Everything downstream then runs unchanged.
  if (glucose_unit == "mmol/L")
    G <- lapply(G, function(x) if (!is.null(x)) mmol_to_mgdl_glucose(x) else NULL)
  if (insulin_unit == "pmol/L")
    I <- lapply(I, function(x) if (!is.null(x)) insulin_pmol_to_uU(x) else NULL)
  if (cpeptide_unit != "nmol/L") {
    cp_div <- if (cpeptide_unit == "ng/mL") 3.02 else 1000
    CP <- lapply(CP, function(x) if (!is.null(x)) x / cp_div else NULL)
  }
  if (!is.null(trig) && trig_unit == "mmol/L") trig <- mmol_to_mgdl_trig(trig)
  if (!is.null(chdl) && hdl_unit == "mmol/L")  chdl <- mmol_to_mgdl_chol(chdl)
  if (glucagon_unit == "pg/mL") {
    if (!is.null(glcg000)) glcg000 <- glcg000 / 3.483
    if (!is.null(glcg030)) glcg030 <- glcg030 / 3.483
  }

  Gm <- lapply(G, function(x) if (!is.null(x)) mgdl_to_mmol_glucose(x) else NULL)
  Ip <- lapply(I, function(x) if (!is.null(x)) insulin_uU_to_pmol(x) else NULL)
  trig_m <- if (!is.null(trig)) mgdl_to_mmol_trig(trig) else NULL
  chdl_m <- if (!is.null(chdl)) mgdl_to_mmol_chol(chdl) else NULL

  ## ---- output accumulator ---------------------------------------------------
  out <- list()
  add <- function(name, value) if (!is.null(value)) out[[name]] <<- value
  present <- function(...) !any(vapply(list(...), is.null, logical(1)))

  ## AUC over a window if all its columns exist, else NULL
  auc_win <- list("15" = c(0, 15), "30" = c(0, 30), "60" = c(0, 30, 60),
                  "90" = c(0, 30, 60, 90), "120" = c(0, 30, 60, 90, 120))
  auc_key <- list("15" = c("0", "15"), "30" = c("0", "30"),
                  "60" = c("0", "30", "60"), "90" = c("0", "30", "60", "90"),
                  "120" = c("0", "30", "60", "90", "120"))
  make_auc <- function(L, incremental = FALSE) {
    res <- list()
    for (w in names(auc_win)) {
      cols <- L[auc_key[[w]]]
      if (any(vapply(cols, is.null, logical(1)))) { res[[w]] <- NULL; next }
      mat <- do.call(cbind, cols)
      res[[w]] <- if (incremental) iauc_trapezoid(auc_win[[w]], mat) else
        auc_trapezoid(auc_win[[w]], mat)
    }
    res
  }
  aG <- make_auc(G);  aI <- make_auc(I);  aCP <- make_auc(CP)
  iG <- make_auc(G, TRUE); iI <- make_auc(I, TRUE); iCP <- make_auc(CP, TRUE)

  ## ---- fasting sensitivity --------------------------------------------------
  if (present(I[["0"]])) { add("fasting_ins", I[["0"]]); add("raynaud", raynaud(I[["0"]])) }
  if (present(G[["0"]])) add("fasting_gluc", G[["0"]])
  if (present(I[["0"]], Gm[["0"]])) {
    hi <- homa_ir(I[["0"]], Gm[["0"]])
    add("homa_ir", hi); add("isi", 1 / hi)
    add("firi", firi(I[["0"]], Gm[["0"]]))
    add("igr", insulin_glucose_ratio(I[["0"]], Gm[["0"]]))
    add("gir", glucose_insulin_ratio(I[["0"]], Gm[["0"]]))
    add("isi_bennett", bennett_isi(I[["0"]], Gm[["0"]]))
  }
  if (present(I[["0"]], G[["0"]])) add("quicki", quicki(I[["0"]], G[["0"]]))
  if (present(I[["0"]], trig_m)) add("mcauley", mcauley(I[["0"]], trig_m))
  if (present(CP[["0"]], Gm[["0"]])) add("ohkura", ohkura(CP[["0"]], Gm[["0"]]))
  if (present(I[["0"]], Gm[["0"]], adipon)) add("homa_ad", homa_ad(I[["0"]], Gm[["0"]], adipon))
  if (present(I[["0"]], Gm[["0"]])) add("belfiore_basal", belfiore_basal(I[["0"]], Gm[["0"]]))
  if (present(hba1c, htn, bmi)) add("egdr_bmi", egdr_bmi(hba1c, htn, bmi))
  if (present(hba1c, htn, waist)) add("egdr_wc", egdr_waist(hba1c, htn, waist))
  if (present(hba1c, htn, whr)) add("egdr_whr", egdr_whr(hba1c, htn, whr))

  ## ---- lipid / anthropometric ----------------------------------------------
  if (present(trig, chdl)) add("tg_hdl", tg_hdl(trig, chdl))
  tyg_v <- if (present(trig, G[["0"]])) tyg(trig, G[["0"]]) else NULL
  if (!is.null(tyg_v)) {
    add("tyg", tyg_v)
    if (present(bmi))   add("tyg_bmi",  tyg_adjusted(tyg_v, bmi))
    if (present(waist)) add("tyg_wc",   tyg_adjusted(tyg_v, waist))
    if (present(whtr))  add("tyg_whtr", tyg_adjusted(tyg_v, whtr))
    if (present(whr))   add("tyg_whr",  tyg_adjusted(tyg_v, whr))
  }
  mets <- if (present(G[["0"]], trig, chdl, bmi)) mets_ir(G[["0"]], trig, chdl, bmi) else NULL
  if (!is.null(mets)) add("mets_ir", mets)
  if (present(mets, whtr, male, age)) add("mets_vf", mets_vf(mets, whtr, male, age))
  if (present(chdl, trig, bmi)) add("spise", spise(chdl, trig, bmi))
  if (present(waist, bmi, trig_m, chdl_m, male)) add("vai", vai(waist, bmi, trig_m, chdl_m, male))
  if (present(waist, trig_m, male)) add("lap", lap(waist, trig_m, male))
  if (present(leptin, adipon)) add("lar", leptin_adiponectin_ratio(leptin, adipon))

  ## ---- OGTT sensitivity -----------------------------------------------------
  gmean <- if (!is.null(aG[["120"]])) aG[["120"]] / 120 else NULL
  imean <- if (!is.null(aI[["120"]])) aI[["120"]] / 120 else NULL
  gmean_mmol <- if (!is.null(gmean)) mgdl_to_mmol_glucose(gmean) else NULL
  if (present(G[["0"]], I[["0"]], gmean, imean))
    add("matsuda", matsuda(G[["0"]], I[["0"]], gmean, imean))
  for (w in c("15", "30", "60", "90")) {
    if (!is.null(aG[[w]]) && !is.null(aI[[w]]) && present(G[["0"]], I[["0"]])) {
      add(paste0("matsuda_", w),
          matsuda(G[["0"]], I[["0"]], aG[[w]] / as.numeric(w), aI[[w]] / as.numeric(w)))
    }
  }
  if (present(G[["0"]], G[["120"]], gmean, imean, weight))
    add("gutt_isi0120", gutt_isi(G[["0"]], G[["120"]], gmean, imean, weight))
  if (present(Gm[["0"]], Gm[["120"]], gmean_mmol, imean, weight))
    add("cederholm", cederholm(Gm[["0"]], Gm[["120"]], gmean_mmol, imean, weight))
  sib <- if (present(I[["0"]], Gm[["0"]], weight)) avignon_sib(I[["0"]], Gm[["0"]], weight) else NULL
  si2h <- if (present(I[["120"]], Gm[["120"]], weight)) avignon_si2h(I[["120"]], Gm[["120"]], weight) else NULL
  if (!is.null(sib)) add("avignon_sib", sib)
  if (!is.null(si2h)) add("avignon_si2h", si2h)
  if (present(sib, si2h)) add("avignon_sim", avignon_sim(sib, si2h))
  if (present(bmi, Ip[["120"]], Gm[["90"]])) add("stumvoll_isi", stumvoll_isi(bmi, Ip[["120"]], Gm[["90"]]))
  if (present(bmi, Ip[["120"]], age)) add("stumvoll_dem", stumvoll_isi_demo(bmi, Ip[["120"]], age))
  if (present(Ip[["0"]], Ip[["120"]], Gm[["120"]]))
    add("stumvoll_modi", stumvoll_isi_bmi_independent(Ip[["0"]], Ip[["120"]], Gm[["120"]]))
  if (present(bmi, Ip[["120"]], Gm[["90"]])) add("stumvoll_mcr", stumvoll_mcr(bmi, Ip[["120"]], Gm[["90"]]))
  if (present(Gm[["0"]], Gm[["90"]], Gm[["120"]], Ip[["0"]], Ip[["90"]], weight, height))
    add("ogis", ogis(Gm[["0"]], Gm[["90"]], Gm[["120"]], Ip[["0"]], Ip[["90"]], weight, height))
  if (present(I[["0"]], I[["60"]], I[["120"]], Gm[["0"]], Gm[["60"]], Gm[["120"]])) {
    iarea <- belfiore_area(I[["0"]], I[["60"]], I[["120"]])
    garea <- belfiore_area(Gm[["0"]], Gm[["60"]], Gm[["120"]])
    add("belfiore_ogtt", belfiore_ogtt(iarea, garea))
  }

  ## ---- secretion ------------------------------------------------------------
  if (present(I[["0"]], Gm[["0"]])) add("homa_b", homa_b(I[["0"]], Gm[["0"]]))
  for (t in c("15", "30", "60", "90", "120")) {
    if (present(I[["0"]], I[[t]], Gm[["0"]], Gm[[t]]))
      add(paste0("igi_", t), igi(I[["0"]], I[[t]], Gm[["0"]], Gm[[t]]))
    if (present(I[[t]], G[[t]])) add(paste0("cir", t), cir(I[[t]], G[[t]]))
  }
  if (present(I[["0"]], I[["30"]], G[["30"]])) add("kadowaki", kadowaki(I[["0"]], I[["30"]], G[["30"]]))
  if (present(Ip[["0"]], Ip[["30"]], Gm[["30"]])) {
    add("stumvoll_1st", stumvoll_first_phase(Ip[["0"]], Ip[["30"]], Gm[["30"]]))
    add("stumvoll_2nd", stumvoll_second_phase(Ip[["0"]], Ip[["30"]], Gm[["30"]]))
  }
  if (present(I[["0"]], I[["30"]], I[["120"]], Gm[["0"]], Gm[["30"]], Gm[["120"]], bmi, male))
    add("bigtt_30_120", bigtt_air_30_120(I[["0"]], I[["30"]], I[["120"]],
                                         Gm[["0"]], Gm[["30"]], Gm[["120"]], bmi, male))
  if (present(I[["0"]], I[["60"]], I[["120"]], Gm[["0"]], Gm[["60"]], Gm[["120"]], bmi, male))
    add("bigtt_60_120", bigtt_air_60_120(I[["0"]], I[["60"]], I[["120"]],
                                         Gm[["0"]], Gm[["60"]], Gm[["120"]], bmi, male))

  ## ---- C-peptide ------------------------------------------------------------
  for (t in c("0", "15", "30", "60", "90", "120")) {
    if (present(CP[[t]], Gm[[t]])) add(paste0("cpep_gluc_", t), cpeptide_glucose_ratio(CP[[t]], Gm[[t]]))
  }
  for (t in c("15", "30", "60", "90", "120")) {
    if (present(CP[["0"]], CP[[t]], Gm[["0"]], Gm[[t]]))
      add(paste0("cpi_", t), cpi(CP[["0"]], CP[[t]], Gm[["0"]], Gm[[t]]))
  }
  if (present(Ip[["0"]], CP[["0"]])) add("icpr", insulin_cpeptide_ratio(Ip[["0"]], CP[["0"]]))
  if (present(aCP[["120"]], aI[["120"]]))
    add("cpep_ins_molar_auc", cpeptide_insulin_auc_ratio(aCP[["120"]], aI[["120"]]))

  ## ---- AUCs and ratios ------------------------------------------------------
  for (w in names(auc_win)) {
    if (!is.null(aG[[w]])) add(paste0("auc_gluc_0_", w), aG[[w]])
    if (!is.null(aI[[w]])) add(paste0("auc_ins_0_", w), aI[[w]])
    if (!is.null(aCP[[w]])) add(paste0("auc_cp_0_", w), aCP[[w]])
    if (!is.null(aI[[w]]) && !is.null(aG[[w]])) add(paste0("auc_ins_gluc_0_", w), aI[[w]] / aG[[w]])
    if (!is.null(aCP[[w]]) && !is.null(aG[[w]])) add(paste0("auc_cp_gluc_0_", w), aCP[[w]] / aG[[w]])
    if (!is.null(iG[[w]])) add(paste0("iauc_gluc_0_", w), iG[[w]])
    if (!is.null(iI[[w]])) add(paste0("iauc_ins_0_", w), iI[[w]])
    if (!is.null(iCP[[w]])) add(paste0("iauc_cp_0_", w), iCP[[w]])
  }
  for (w in c("30", "60", "90", "120")) {
    if (!is.null(iI[[w]]) && !is.null(iG[[w]])) add(paste0("iauc_ins_gluc_0_", w), guard_ratio(iI[[w]], iG[[w]]))
    if (!is.null(iCP[[w]]) && !is.null(iG[[w]])) add(paste0("iauc_cp_gluc_0_", w), guard_ratio(iCP[[w]], iG[[w]]))
  }

  ## ---- disposition ----------------------------------------------------------
  igi30 <- out[["igi_30"]]
  if (!is.null(igi30) && present(I[["0"]])) add("odi", disposition_index(igi30, 1 / I[["0"]]))
  if (!is.null(igi30) && !is.null(out[["matsuda"]])) add("odi_matsuda", disposition_index(igi30, out[["matsuda"]]))
  if (!is.null(igi30) && !is.null(out[["matsuda_30"]])) add("odi_matsuda_30", disposition_index(igi30, out[["matsuda_30"]]))
  if (!is.null(igi30) && !is.null(out[["quicki"]])) add("odi_quicki", disposition_index(igi30, out[["quicki"]]))
  if (!is.null(igi30) && !is.null(out[["isi"]])) add("odi_isi", disposition_index(igi30, out[["isi"]]))
  if (!is.null(out[["cpi_30"]]) && !is.null(out[["matsuda"]])) add("odi_cp", disposition_index(out[["cpi_30"]], out[["matsuda"]]))
  if (!is.null(out[["auc_ins_gluc_0_120"]]) && !is.null(out[["matsuda"]]))
    add("issi2", disposition_index(out[["auc_ins_gluc_0_120"]], out[["matsuda"]]))

  ## ---- glucagon -------------------------------------------------------------
  if (present(glcg000)) add("fasting_glucagon", glcg000)
  if (present(glcg000, Ip[["0"]])) add("glucagon_ins_ratio", glucagon_insulin_ratio(glcg000, Ip[["0"]]))
  if (present(glcg000, glcg030)) {
    add("glucagon_suppr_abs", glucagon_suppression(glcg000, glcg030))
    add("glucagon_suppr_pct", glucagon_suppression(glcg000, glcg030, percent = TRUE))
  }

  ## ---- proinsulin -----------------------------------------------------------
  if (present(pin)) add("fasting_proinsulin", pin)
  if (present(pin, Ip[["0"]])) add("proinsulin_insulin_ratio", pin / Ip[["0"]])

  ## ---- HOMA2 (built-in Oxford DTU tables) -----------------------------------
  ## Insulin is converted uU/mL -> pmol/L at 6.945 (matching the DPP pipeline);
  ## see [homa2_insulin()] to use the 6.0 DTU convention instead.
  if (present(Gm[["0"]], Ip[["0"]])) {
    h <- homa2_insulin(Gm[["0"]], Ip[["0"]])
    add("homa2_ir", h$homa2_ir); add("homa2_b_ins", h$homa2_b); add("homa2_s", h$homa2_s)
  }
  if (present(Gm[["0"]], CP[["0"]])) {
    hc <- homa2_cpeptide(Gm[["0"]], CP[["0"]])
    add("homa2_ir_cp", hc$homa2_ir); add("homa2_b_cp", hc$homa2_b); add("homa2_s_cp", hc$homa2_s)
  }

  ## ---- assemble -------------------------------------------------------------
  for (nm in names(out)) data[[nm]] <- out[[nm]]
  attr(data, "surinsulot_computed") <- names(out)

  ## ---- optional Yeo-Johnson transform of every computed index ---------------
  if (transform) {
    lambdas <- numeric(0)
    for (nm in names(out)) {
      lam <- yeo_johnson_lambda(data[[nm]])
      lambdas[nm] <- lam
      data[[paste0(nm, "_yj")]] <-
        if (is.na(lam)) NA_real_ else yeo_johnson(data[[nm]], lam)
    }
    attr(data, "surinsulot_yj_lambda") <- lambdas
  }

  if (verbose) {
    message(sprintf("SurInsulot: computed %d surrogate index columns%s.",
                    length(out),
                    if (transform) sprintf(" (+ %d _yj columns)", length(out)) else ""))
  }
  data
}
