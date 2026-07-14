test_that("secretion indices match reference values", {
  expect_equal(homa_b(subj$i000, subj$g000m), EXP$homa_b, tolerance = TOL)
  expect_equal(igi(subj$i000, subj$i030, subj$g000m, subj$g030m), EXP$igi_30, tolerance = TOL)
  expect_equal(igi(subj$i000, subj$i120, subj$g000m, subj$g120m), EXP$igi_120, tolerance = TOL)
  expect_equal(cir(subj$i030, subj$g030), EXP$cir30, tolerance = TOL)
  expect_equal(kadowaki(subj$i000, subj$i030, subj$g030), EXP$kadowaki, tolerance = TOL)
  expect_equal(stumvoll_first_phase(subj$i000p, subj$i030p, subj$g030m),
               EXP$stumvoll_1st, tolerance = TOL)
  expect_equal(stumvoll_second_phase(subj$i000p, subj$i030p, subj$g030m),
               EXP$stumvoll_2nd, tolerance = TOL)
})

test_that("cir returns NA when glucose <= 70 mg/dL", {
  expect_true(is.na(cir(50, 65)))
  expect_false(is.na(cir(50, 90)))
})

test_that("BIGTT-AIR forms match reference values", {
  expect_equal(
    bigtt_air_30_120(subj$i000, subj$i030, subj$i120,
                     subj$g000m, subj$g030m, subj$g120m, subj$bmi, subj$male),
    EXP$bigtt_30_120, tolerance = TOL)
  expect_equal(
    bigtt_air_60_120(subj$i000, subj$i060, subj$i120,
                     subj$g000m, subj$g060m, subj$g120m, subj$bmi, subj$male),
    EXP$bigtt_60_120, tolerance = TOL)
})
