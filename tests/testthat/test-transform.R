test_that("yeo_johnson is identity at lambda = 1 and passes NA/Inf through", {
  x <- c(-3, -1, 0, 2, 5)
  expect_equal(yeo_johnson(x, 1), x)
  y <- yeo_johnson(c(1, NA, Inf, -2), 0.5)
  expect_true(is.na(y[2]))
  expect_true(is.infinite(y[3]))
})

test_that("yeo_johnson matches the closed-form values", {
  expect_equal(yeo_johnson(2, 0.5), (3^0.5 - 1) / 0.5, tolerance = 1e-8)
  expect_equal(yeo_johnson(-2, 0.5), -((3^1.5 - 1) / 1.5), tolerance = 1e-8)
  expect_equal(yeo_johnson(2, 0), log1p(2), tolerance = 1e-8)
  expect_equal(yeo_johnson(-2, 2), -log1p(2), tolerance = 1e-8)
})

test_that("yeo_johnson is monotonic increasing", {
  x <- seq(-6, 6, by = 0.25)
  expect_true(all(diff(yeo_johnson(x, 0.3)) > 0))
})

test_that("yeo_johnson_lambda reduces skew and guards degenerate input", {
  set.seed(1)
  x <- rexp(500)
  lam <- yeo_johnson_lambda(x)
  expect_false(is.na(lam))
  sk <- function(v) { v <- v[is.finite(v)]; mean((v - mean(v))^3) / sd(v)^3 }
  expect_lt(abs(sk(yeo_johnson(x, lam))), abs(sk(x)))
  expect_true(is.na(yeo_johnson_lambda(1:5)))         # < 10 finite values
  expect_true(is.na(yeo_johnson_lambda(rep(2, 50))))  # zero variance
})

test_that("compute_surrogates(transform = TRUE) appends matching _yj columns", {
  set.seed(1)
  n <- 40
  df <- data.frame(
    g000 = rnorm(n, 100, 12), g030 = rnorm(n, 160, 15), g120 = rnorm(n, 135, 15),
    i000 = rnorm(n, 12, 3),   i030 = rnorm(n, 85, 20),
    trig = rnorm(n, 150, 30), chdl = rnorm(n, 45, 8), bmi = rnorm(n, 30, 4),
    sex = rep(1:2, length.out = n)
  )
  tr <- compute_surrogates(df, verbose = FALSE, transform = TRUE)
  expect_true(all(c("homa_ir_yj", "tyg_yj", "igi_30_yj") %in% names(tr)))

  lam <- attr(tr, "surinsulot_yj_lambda")
  expect_true("homa_ir" %in% names(lam))
  expect_equal(tr$homa_ir_yj, yeo_johnson(tr$homa_ir, lam[["homa_ir"]]),
               tolerance = 1e-8)

  # default (no transform) adds no _yj columns
  base <- compute_surrogates(df, verbose = FALSE)
  expect_false(any(grepl("_yj$", names(base))))
})
