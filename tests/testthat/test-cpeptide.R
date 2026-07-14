test_that("C-peptide indices match reference values", {
  expect_equal(cpeptide_glucose_ratio(subj$cp000, subj$g000m), EXP$cpep_gluc_0, tolerance = TOL)
  expect_equal(cpeptide_glucose_ratio(subj$cp030, subj$g030m), EXP$cpep_gluc_30, tolerance = TOL)
  expect_equal(cpi(subj$cp000, subj$cp030, subj$g000m, subj$g030m), EXP$cpi_30, tolerance = TOL)
  expect_equal(insulin_cpeptide_ratio(subj$i000p, subj$cp000), EXP$icpr, tolerance = TOL)

  cpcurve <- c(subj$cp000, subj$cp030, subj$cp060, subj$cp090, subj$cp120)
  icurve  <- c(subj$i000, subj$i030, subj$i060, subj$i090, subj$i120)
  aCP <- auc_trapezoid(c(0, 30, 60, 90, 120), cpcurve)
  aI  <- auc_trapezoid(c(0, 30, 60, 90, 120), icurve)
  expect_equal(cpeptide_insulin_auc_ratio(aCP, aI), EXP$cpep_ins_molar_auc, tolerance = TOL)
})

test_that("insulin_cpeptide_ratio returns a plausible molar ratio", {
  # co-secreted 1:1 but C-peptide clears slower -> insulin:C-peptide well below 1
  expect_lt(insulin_cpeptide_ratio(subj$i000p, subj$cp000), 1)
})
