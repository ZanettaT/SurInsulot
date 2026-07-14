test_that("lipid / anthropometric indices match reference values", {
  expect_equal(tg_hdl(subj$trig, subj$chdl), EXP$tg_hdl, tolerance = TOL)
  tyg_v <- tyg(subj$trig, subj$g000)
  expect_equal(tyg_v, EXP$tyg, tolerance = TOL)
  expect_equal(tyg_adjusted(tyg_v, subj$bmi), EXP$tyg_bmi, tolerance = TOL)
  mir <- mets_ir(subj$g000, subj$trig, subj$chdl, subj$bmi)
  expect_equal(mir, EXP$mets_ir, tolerance = TOL)
  expect_equal(mets_vf(mir, subj$whtr, subj$male, subj$age), EXP$mets_vf, tolerance = TOL)
  expect_equal(spise(subj$chdl, subj$trig, subj$bmi), EXP$spise, tolerance = TOL)
  expect_equal(vai(subj$waist, subj$bmi, subj$trig_m, subj$chdl_m, subj$male),
               EXP$vai, tolerance = TOL)
  expect_equal(lap(subj$waist, subj$trig_m, subj$male), EXP$lap, tolerance = TOL)
  expect_equal(leptin_adiponectin_ratio(subj$leptin, subj$adipon), EXP$lar, tolerance = TOL)
})

test_that("sex-specific indices differ by sex and NA on unknown codes", {
  expect_false(isTRUE(all.equal(
    vai(subj$waist, subj$bmi, subj$trig_m, subj$chdl_m, male = 1),
    vai(subj$waist, subj$bmi, subj$trig_m, subj$chdl_m, male = 0))))
  expect_true(is.na(lap(subj$waist, subj$trig_m, male = 2)))  # raw sex code -> NA
})
