# HOMA2 grid-point values reproduce the Oxford DTU HOMA2 Calculator v2.2.4
# exactly (these anchors are the homa2calc validation set). Because the tables
# are the implementation, these are exact, not tolerance-based, checks.

test_that("homa2_insulin matches DTU calculator at grid points", {
  r <- homa2_insulin(3.0, 20)
  expect_equal(r$homa2_b, 139.8, tolerance = 1e-6)
  expect_equal(r$homa2_s, 308.9, tolerance = 1e-6)

  r2 <- homa2_insulin(5.2, 58)
  expect_equal(r2$homa2_b, 93.0, tolerance = 1e-6)
  expect_equal(r2$homa2_s, 91.4, tolerance = 1e-6)

  r3 <- homa2_insulin(14.0, 96)
  expect_equal(r3$homa2_b, 21.7, tolerance = 1e-6)
  expect_equal(r3$homa2_s, 42.7, tolerance = 1e-6)
})

test_that("homa2_cpeptide and homa2_specific_insulin match at grid points", {
  r <- homa2_cpeptide(3.0, 0.2)
  expect_equal(r$homa2_b, 151.8, tolerance = 1e-6)
  expect_equal(r$homa2_s, 272.6, tolerance = 1e-6)

  r2 <- homa2_cpeptide(14.0, 3.5)
  expect_equal(r2$homa2_b, 97.2, tolerance = 1e-6)
  expect_equal(r2$homa2_s, 8.9, tolerance = 1e-6)

  r3 <- homa2_specific_insulin(3.0, 20)
  expect_equal(r3$homa2_b, 150.3, tolerance = 1e-6)
  expect_equal(r3$homa2_s, 276.9, tolerance = 1e-6)
})

test_that("bilinear interpolation between grid points is correct", {
  # midway example: verified independently, NOT the (incorrect) homa2calc
  # README value of 131.9/92.8
  r <- homa2_insulin(5.0, 60)
  expect_equal(r$homa2_b, 105.6, tolerance = 1e-6)
  expect_equal(r$homa2_s, 90.1, tolerance = 1e-6)
})

test_that("out-of-range inputs return NA", {
  expect_true(is.na(homa2_insulin(2.0, 50)$homa2_ir))    # glucose < 3.0
  expect_true(is.na(homa2_insulin(5.0, 500)$homa2_ir))   # insulin > 400
  expect_true(is.na(homa2_cpeptide(5.0, 5.0)$homa2_b))   # cpeptide > 3.5
})

test_that("HOMA2 output is a vectorised data frame with IR = 100 / %S", {
  res <- homa2_insulin(c(4.0, 6.0, 8.5), c(60, 80, 120))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
  expect_named(res, c("homa2_b", "homa2_s", "homa2_ir"))
  # homa2_ir = 100 / (unrounded %S); homa2_s is %S rounded to 0.1. Recover %S
  # from IR and confirm it matches the reported (rounded) value.
  expect_equal(round(100 / res$homa2_ir, 1), res$homa2_s)
})
