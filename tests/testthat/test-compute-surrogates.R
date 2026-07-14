full_df <- data.frame(
  g000 = subj$g000, g015 = subj$g015, g030 = subj$g030,
  g060 = subj$g060, g090 = subj$g090, g120 = subj$g120,
  i000 = subj$i000, i015 = subj$i015, i030 = subj$i030,
  i060 = subj$i060, i090 = subj$i090, i120 = subj$i120,
  cp000 = subj$cp000, cp015 = subj$cp015, cp030 = subj$cp030,
  cp060 = subj$cp060, cp090 = subj$cp090, cp120 = subj$cp120,
  glcg000 = subj$glcg000, glcg030 = subj$glcg030,
  trig = subj$trig, chdl = subj$chdl,
  bmi = subj$bmi, weight_kg = subj$weight, waist_cm = subj$waist,
  hip_cm = subj$hip, height_cm = subj$height, whr = subj$whr, whtr = subj$whtr,
  hba1c_perc = subj$hba1c, adipon = subj$adipon, leptin = subj$leptin,
  pin = subj$pin, sex = subj$male, age = subj$age, hypertension = subj$htn
)

test_that("compute_surrogates reproduces the individual index values", {
  out <- compute_surrogates(full_df, verbose = FALSE)
  checks <- c("homa_ir", "quicki", "tyg", "tyg_bmi", "spise", "vai",
              "matsuda", "matsuda_30", "gutt_isi0120", "cederholm", "ogis",
              "belfiore_ogtt", "homa_b", "igi_30", "cir30", "kadowaki",
              "stumvoll_1st", "bigtt_30_120", "icpr", "cpi_30",
              "odi_matsuda", "issi2", "glucagon_suppr_pct", "auc_gluc_0_120")
  ref <- c(homa_ir = EXP$homa_ir, quicki = EXP$quicki, tyg = EXP$tyg,
           tyg_bmi = EXP$tyg_bmi, spise = EXP$spise, vai = EXP$vai,
           matsuda = EXP$matsuda, matsuda_30 = EXP$matsuda_30,
           gutt_isi0120 = EXP$gutt_isi0120, cederholm = EXP$cederholm,
           ogis = EXP$ogis, belfiore_ogtt = EXP$belfiore_ogtt,
           homa_b = EXP$homa_b, igi_30 = EXP$igi_30, cir30 = EXP$cir30,
           kadowaki = EXP$kadowaki, stumvoll_1st = EXP$stumvoll_1st,
           bigtt_30_120 = EXP$bigtt_30_120, icpr = EXP$icpr, cpi_30 = EXP$cpi_30,
           odi_matsuda = EXP$odi_matsuda, issi2 = EXP$issi2,
           glucagon_suppr_pct = EXP$glucagon_suppr_pct,
           auc_gluc_0_120 = EXP$auc_gluc_0_120)
  for (nm in checks) {
    expect_true(nm %in% names(out), info = nm)
    expect_equal(out[[nm]][1], unname(ref[nm]), tolerance = TOL, info = nm)
  }
})

test_that("compute_surrogates preserves rows and records computed names", {
  out <- compute_surrogates(full_df, verbose = FALSE)
  expect_equal(nrow(out), nrow(full_df))
  expect_true(all(names(full_df) %in% names(out)))
  expect_true(length(attr(out, "surinsulot_computed")) > 40)
})

test_that("fasting-only input computes fasting indices and skips OGTT ones", {
  fast_df <- data.frame(
    g000 = c(95, 110), i000 = c(10, 18),
    trig = c(150, 200), chdl = c(45, 38), bmi = c(29, 34), sex = c(1, 2)
  )
  out <- compute_surrogates(fast_df, verbose = FALSE)
  expect_true(all(c("homa_ir", "quicki", "tyg", "tg_hdl", "spise") %in% names(out)))
  expect_false("matsuda" %in% names(out))
  expect_false("igi_30" %in% names(out))
  expect_equal(out$homa_ir, c(10 * 95 / 18 / 22.5, 18 * 110 / 18 / 22.5), tolerance = TOL)
})

test_that("declared input units convert and give identical results", {
  base <- compute_surrogates(full_df, verbose = FALSE)

  # re-express the same data in mmol/L glucose, pmol/L insulin, ng/mL C-peptide,
  # mmol/L lipids, pg/mL glucagon
  alt <- full_df
  gcols <- paste0("g", c("000", "015", "030", "060", "090", "120"))
  icols <- paste0("i", c("000", "015", "030", "060", "090", "120"))
  cpcols <- paste0("cp", c("000", "015", "030", "060", "090", "120"))
  alt[gcols]  <- lapply(alt[gcols],  function(x) x / 18)      # mg/dL -> mmol/L
  alt[icols]  <- lapply(alt[icols],  function(x) x * 6.945)   # uU/mL -> pmol/L
  alt[cpcols] <- lapply(alt[cpcols], function(x) x * 3.02)    # nmol/L -> ng/mL
  alt$trig <- alt$trig / 88.57                                # mg/dL -> mmol/L
  alt$chdl <- alt$chdl / 38.67                                # mg/dL -> mmol/L
  alt$glcg000 <- alt$glcg000 * 3.483                          # pmol/L -> pg/mL
  alt$glcg030 <- alt$glcg030 * 3.483

  out <- compute_surrogates(
    alt, verbose = FALSE,
    glucose_unit = "mmol/L", insulin_unit = "pmol/L", cpeptide_unit = "ng/mL",
    trig_unit = "mmol/L", hdl_unit = "mmol/L", glucagon_unit = "pg/mL"
  )

  for (nm in attr(base, "surinsulot_computed")) {
    expect_equal(out[[nm]], base[[nm]], tolerance = 1e-6, info = nm)
  }
})

test_that("invalid unit strings are rejected", {
  expect_error(compute_surrogates(full_df, glucose_unit = "mmHg", verbose = FALSE))
})

test_that("HOMA2 indices are computed from the built-in DTU tables", {
  out <- compute_surrogates(full_df, verbose = FALSE)
  expect_true(all(c("homa2_ir", "homa2_b_ins", "homa2_s",
                    "homa2_ir_cp", "homa2_b_cp", "homa2_s_cp") %in% names(out)))
  # subject: g000 100 mg/dL -> 5.556 mmol/L; i000 12 uU/mL -> 83.34 pmol/L
  expect_equal(out$homa2_b_ins[1], 106.8, tolerance = 1e-4)
  expect_equal(out$homa2_s[1], 63.8, tolerance = 1e-4)
  expect_equal(out$homa2_ir[1], 1.5676, tolerance = 1e-3)
  # C-peptide route (cp000 = 0.6 nmol/L)
  expect_equal(out$homa2_b_cp[1], 96.0, tolerance = 1e-4)
  expect_equal(out$homa2_s_cp[1], 74.6, tolerance = 1e-4)
})
