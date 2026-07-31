#' Build an empirical residual pool from a historical price series
#'
#' Fits a no-intercept AR(1) to first-differenced, discounted log-prices,
#' then resamples the centered residuals with a stationary block bootstrap
#' (Politis and Romano, 1994). The resulting pool is the reusable input to
#' \code{\link{eq_price}} and \code{\link{eq_price_asian}}: build it once per
#' discount rate and price as many strikes, maturities, and payoff types as
#' needed from it.
#'
#' @param prices Numeric vector of historical prices, ordered by date, with
#'   no missing values (interpolate or drop \code{NA}s beforehand).
#' @param r Discount / funding rate (annualized, continuously compounded)
#'   used to discount \code{prices} into log-space before differencing.
#' @param dt Time step between observations, in years (default \code{1/252}
#'   for daily data).
#' @param window Optional trailing window length, in observations, to use
#'   instead of the full history (e.g. \code{252} for a trailing one-year
#'   window). \code{NULL} (default) uses the full series.
#' @param n_sims Number of bootstrap paths to draw (default \code{10000}).
#' @param seed Random seed used for the bootstrap draw, for reproducibility.
#'
#' @return An object of class \code{"eq_pool"}: a list with the resampled
#'   residual matrix (\code{n_w} rows by \code{n_sims} columns), the fitted
#'   AR(1) coefficient, and diagnostic statistics (Ljung-Box and ARCH-LM
#'   p-values on the pre-bootstrap residuals, annualized volatility, skew,
#'   and excess kurtosis).
#'
#' @references
#' Politis, D. N. and Romano, J. P. (1994). The stationary bootstrap.
#' \emph{Journal of the American Statistical Association}, 89(428), 1303-1313.
#'
#' @examples
#' set.seed(1)
#' prices <- cumprod(1 + rnorm(500, 0.0003, 0.01)) * 100
#' pool <- eq_pool(prices, r = 0.03)
#' pool
#'
#' @export
eq_pool <- function(prices, r, dt = 1 / 252, window = NULL,
                     n_sims = 10000L, seed = 123) {
  stopifnot(is.numeric(prices), all(prices > 0), !anyNA(prices))
  stopifnot(length(prices) >= 3)

  log_p <- log(prices)
  log_w <- if (!is.null(window)) utils::tail(log_p, window + 1) else log_p

  t_index <- seq(0, by = dt, length.out = length(log_w))
  log_disc <- log_w - r * t_index

  r_t <- diff(log_disc)
  n_w <- length(r_t)
  if (n_w < 3) stop("Not enough observations to fit an AR(1) after differencing.")

  X <- matrix(r_t[-length(r_t)], ncol = 1)
  yy <- r_t[-1]
  phi <- stats::.lm.fit(X, yy)$coefficients

  eps_t <- numeric(n_w)
  eps_t[-1] <- yy - X * phi
  eps_t[1] <- r_t[1] - mean(r_t)
  eps_t <- eps_t - mean(eps_t)

  set.seed(seed)
  resampled <- tseries::tsbootstrap(x = eps_t, nb = n_sims, type = "stationary")
  resampled <- resampled - rowMeans(resampled)

  lb_p <- stats::Box.test(eps_t, type = "Ljung-Box")$p.value
  arch_p <- tryCatch(
    FinTS::ArchTest(eps_t)$p.value,
    error = function(e) NA_real_
  )

  structure(
    list(
      resampled_resids = resampled,
      n_w = n_w,
      phi = phi,
      dt = dt,
      r = r,
      lb_p = lb_p,
      arch_p = arch_p,
      ann_vol = stats::sd(r_t) * sqrt(1 / dt),
      skew = mean(eps_t^3) / stats::sd(eps_t)^3,
      excess_kurt = mean(eps_t^4) / stats::sd(eps_t)^4 - 3,
      n_sims = n_sims,
      seed = seed
    ),
    class = "eq_pool"
  )
}

#' @export
print.eq_pool <- function(x, ...) {
  cat("<empirical residual pool>\n")
  cat(sprintf("  observations: %d | paths: %d | rate: %.4f\n",
              x$n_w, x$n_sims, x$r))
  cat(sprintf("  AR(1) phi: %.4f | ann. vol: %.4f\n", x$phi, x$ann_vol))
  cat(sprintf("  Ljung-Box p: %.4f | ARCH-LM p: %.2e\n", x$lb_p, x$arch_p))
  cat(sprintf("  skew: %.4f | excess kurtosis: %.4f\n", x$skew, x$excess_kurt))
  invisible(x)
}
