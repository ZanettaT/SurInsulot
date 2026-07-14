gt <- c(0, 30, 60, 90, 120)
gcurve <- c(subj$g000, subj$g030, subj$g060, subj$g090, subj$g120)
icurve <- c(subj$i000, subj$i030, subj$i060, subj$i090, subj$i120)
gmean <- auc_trapezoid(gt, gcurve) / 120
imean <- auc_trapezoid(gt, icurve) / 120

test_that("Matsuda (full and windowed) match reference values", {
  expect_equal(matsuda(subj$g000, subj$i000, gmean, imean), EXP$matsuda, tolerance = TOL)
  expect_equal(matsuda_from_curve(gt, gcurve, icurve), EXP$matsuda, tolerance = TOL)
  expect_equal(matsuda_from_curve(gt, gcurve, icurve, upto = 30), EXP$matsuda_30, tolerance = TOL)
  expect_equal(matsuda_from_curve(gt, gcurve, icurve, upto = 60), EXP$matsuda_60, tolerance = TOL)
})

test_that("Gutt and Cederholm match reference values", {
  expect_equal(gutt_isi(subj$g000, subj$g120, gmean, imean, subj$weight),
               EXP$gutt_isi0120, tolerance = TOL)
  expect_equal(cederholm(subj$g000m, subj$g120m, gmean / 18, imean, subj$weight),
               EXP$cederholm, tolerance = TOL)
})

test_that("Stumvoll ISI family matches reference values", {
  expect_equal(stumvoll_isi(subj$bmi, subj$i120p, subj$g090m), EXP$stumvoll_isi, tolerance = TOL)
  expect_equal(stumvoll_isi_demo(subj$bmi, subj$i120p, subj$age), EXP$stumvoll_dem, tolerance = TOL)
  expect_equal(stumvoll_isi_bmi_independent(subj$i000p, subj$i120p, subj$g120m),
               EXP$stumvoll_modi, tolerance = TOL)
  expect_equal(stumvoll_mcr(subj$bmi, subj$i120p, subj$g090m), EXP$stumvoll_mcr, tolerance = TOL)
})

test_that("Avignon and Belfiore OGTT match reference values", {
  sib  <- avignon_sib(subj$i000, subj$g000m, subj$weight)
  si2h <- avignon_si2h(subj$i120, subj$g120m, subj$weight)
  expect_equal(sib, EXP$avignon_sib, tolerance = TOL)
  expect_equal(si2h, EXP$avignon_si2h, tolerance = TOL)
  expect_equal(avignon_sim(sib, si2h), EXP$avignon_sim, tolerance = TOL)
  iarea <- belfiore_area(subj$i000, subj$i060, subj$i120)
  garea <- belfiore_area(subj$g000m, subj$g060m, subj$g120m)
  expect_equal(belfiore_ogtt(iarea, garea), EXP$belfiore_ogtt, tolerance = TOL)
})

test_that("OGIS matches reference value", {
  expect_equal(
    ogis(subj$g000m, subj$g090m, subj$g120m, subj$i000p, subj$i090p,
         subj$weight, subj$height),
    EXP$ogis, tolerance = TOL)
})
