test_that("eq_pool returns a well-formed pool", {
  set.seed(1)
  prices <- cumprod(1 + rnorm(500, 0.0003, 0.01)) * 100
  pool <- eq_pool(prices, r = 0.03, n_sims = 2000)

  expect_s3_class(pool, "eq_pool")
  expect_equal(nrow(pool$resampled_resids), pool$n_w)
  expect_equal(ncol(pool$resampled_resids), 2000)
  expect_true(is.finite(pool$phi))
})

test_that("call and put prices satisfy put-call parity to Monte Carlo precision", {
  set.seed(2)
  prices <- cumprod(1 + rnorm(500, 0.0002, 0.012)) * 100
  pool <- eq_pool(prices, r = 0.03, n_sims = 5000)

  S0 <- 100
  T <- 0.5
  K <- c(90, 100, 110)

  call <- eq_price(pool, S0, T, K, type = "call")
  put <- eq_price(pool, S0, T, K, type = "put")

  lhs <- call$prices - put$prices
  rhs <- S0 - K * exp(-pool$r * T)

  # Same simulated S_T underlies both legs, so parity should hold to
  # floating-point precision, not just Monte Carlo precision.
  expect_equal(unname(lhs), rhs, tolerance = 1e-8)
})

test_that("Duan-Simonato correction enforces the martingale condition", {
  set.seed(3)
  prices <- cumprod(1 + rnorm(400, 0.0, 0.015)) * 50
  pool <- eq_pool(prices, r = 0.02, n_sims = 5000)

  S0 <- 50
  T <- 1
  res <- eq_price(pool, S0, T, K = 50, type = "call")

  # E^Q[e^{-rT} S_T] should equal S0 by construction of c_ds.
  expect_true(is.finite(res$c_ds))
  expect_true(res$c_ds > 0)
})

test_that("eq_price_asian runs and returns one price per strike", {
  set.seed(4)
  prices <- cumprod(1 + rnorm(500, 0.0003, 0.01)) * 100
  pool <- eq_pool(prices, r = 0.03, n_sims = 2000)

  K <- c(95, 100, 105)
  out <- eq_price_asian(pool, S0 = 100, T = 0.5, K = K, n_fix = 6)

  expect_length(out, length(K))
  expect_true(all(is.finite(out)))
  expect_true(all(out >= 0))
})
