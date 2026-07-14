test_that("disposition index is the product of secretion and sensitivity", {
  expect_equal(disposition_index(EXP$igi_30, 1 / subj$i000), EXP$odi, tolerance = TOL)
  expect_equal(disposition_index(EXP$igi_30, EXP$matsuda), EXP$odi_matsuda, tolerance = TOL)
  expect_equal(disposition_index(EXP$igi_30, EXP$matsuda_30), EXP$odi_matsuda_30, tolerance = TOL)
  expect_equal(disposition_index(EXP$igi_30, EXP$quicki), EXP$odi_quicki, tolerance = TOL)
  expect_equal(disposition_index(EXP$igi_30, EXP$isi), EXP$odi_isi, tolerance = TOL)
  expect_equal(disposition_index(EXP$cpi_30, EXP$matsuda), EXP$odi_cp, tolerance = TOL)
  expect_equal(disposition_index(EXP$auc_ins_gluc_0_120, EXP$matsuda), EXP$issi2, tolerance = TOL)
})

test_that("glucagon indices match reference values", {
  expect_equal(glucagon_insulin_ratio(subj$glcg000, subj$i000p),
               EXP$glucagon_ins_ratio, tolerance = TOL)
  expect_equal(glucagon_suppression(subj$glcg000, subj$glcg030),
               EXP$glucagon_suppr_abs, tolerance = TOL)
  expect_equal(glucagon_suppression(subj$glcg000, subj$glcg030, percent = TRUE),
               EXP$glucagon_suppr_pct, tolerance = TOL)
})
