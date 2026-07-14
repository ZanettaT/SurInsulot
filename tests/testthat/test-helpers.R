test_that("unit conversions round-trip and match known values", {
  expect_equal(mgdl_to_mmol_glucose(90), 5, tolerance = TOL)
  expect_equal(mmol_to_mgdl_glucose(5), 90, tolerance = TOL)
  expect_equal(insulin_uU_to_pmol(10), 69.45, tolerance = TOL)
  expect_equal(insulin_pmol_to_uU(69.45), 10, tolerance = TOL)
  expect_equal(mmol_to_mgdl_glucose(mgdl_to_mmol_glucose(123)), 123, tolerance = TOL)
  expect_equal(insulin_pmol_to_uU(insulin_uU_to_pmol(7.3)), 7.3, tolerance = TOL)
  expect_true(is.na(mgdl_to_mmol_glucose(NA_real_)))
})

test_that("auc_trapezoid matches hand calculation for one and many subjects", {
  expect_equal(auc_trapezoid(c(0, 30), c(100, 160)), 3900, tolerance = TOL)
  full <- auc_trapezoid(c(0, 30, 60, 90, 120), c(100, 160, 150, 140, 130))
  expect_equal(full, 16950, tolerance = TOL)
  m <- rbind(c(90, 150, 110), c(100, 180, 130))
  expect_equal(auc_trapezoid(c(0, 30, 120), m), c(15300, 18150), tolerance = TOL)
})

test_that("iauc_trapezoid subtracts baseline", {
  expect_equal(iauc_trapezoid(c(0, 30), c(100, 160)), 900, tolerance = TOL)
  # baseline 12 -> increments c(0,78,73,58,43); net iAUC over 0/30/60/90/120
  expect_equal(iauc_trapezoid(c(0, 30, 60, 90, 120), c(12, 90, 85, 70, 55)),
               6915, tolerance = TOL)
  expect_gt(iauc_trapezoid(c(0, 30), c(12, 90)), 0)
})

test_that("auc errors on mismatched lengths", {
  expect_error(auc_trapezoid(c(0, 30, 60), c(1, 2)))
  expect_error(auc_trapezoid(5, c(1)))
})

test_that("guard_ratio and row_mean behave", {
  expect_equal(guard_ratio(c(10, 20, 5), c(2, 0, -1)),
               c(5, NA, NA), tolerance = TOL)
  expect_equal(row_mean(c(140, NA, 130), c(138, NA, 132)),
               c(139, NA, 131), tolerance = TOL)
})
