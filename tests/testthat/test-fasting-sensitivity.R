test_that("fasting sensitivity indices match reference values", {
  expect_equal(raynaud(subj$i000), EXP$raynaud, tolerance = TOL)
  expect_equal(homa_ir(subj$i000, subj$g000m), EXP$homa_ir, tolerance = TOL)
  expect_equal(1 / homa_ir(subj$i000, subj$g000m), EXP$isi, tolerance = TOL)
  expect_equal(firi(subj$i000, subj$g000m), EXP$firi, tolerance = TOL)
  expect_equal(quicki(subj$i000, subj$g000), EXP$quicki, tolerance = TOL)
  expect_equal(insulin_glucose_ratio(subj$i000, subj$g000m), EXP$igr, tolerance = TOL)
  expect_equal(glucose_insulin_ratio(subj$i000, subj$g000m), EXP$gir, tolerance = TOL)
  expect_equal(bennett_isi(subj$i000, subj$g000m), EXP$isi_bennett, tolerance = TOL)
  expect_equal(mcauley(subj$i000, subj$trig_m), EXP$mcauley, tolerance = TOL)
  expect_equal(ohkura(subj$cp000, subj$g000m), EXP$ohkura, tolerance = TOL)
  expect_equal(homa_ad(subj$i000, subj$g000m, subj$adipon), EXP$homa_ad, tolerance = TOL)
  expect_equal(belfiore_basal(subj$i000, subj$g000m), EXP$belfiore_basal, tolerance = TOL)
})

test_that("eGDR forms match reference values", {
  expect_equal(egdr_bmi(subj$hba1c, subj$htn, subj$bmi), EXP$egdr_bmi, tolerance = TOL)
  expect_equal(egdr_waist(subj$hba1c, subj$htn, subj$waist), EXP$egdr_wc, tolerance = TOL)
  expect_equal(egdr_whr(subj$hba1c, subj$htn, subj$whr), EXP$egdr_whr, tolerance = TOL)
})

test_that("fasting indices are vectorised", {
  expect_equal(homa_ir(c(10, 12), c(5, 5.5)),
               c(10 * 5 / 22.5, 12 * 5.5 / 22.5), tolerance = TOL)
})
