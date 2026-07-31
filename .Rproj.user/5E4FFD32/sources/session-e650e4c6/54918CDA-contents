#' Price European options from an empirical residual pool
#'
#' Reconstructs terminal prices at maturity \code{T} by cumulating bootstrap
#' residuals in log-space, applies a Duan-Simonato (1998) scalar martingale
#' correction, and Monte Carlo-prices the requested strikes. Calls and puts
#' are priced from the \emph{same} simulated terminal draws, so put-call
#' parity holds up to Monte Carlo noise by construction.
#'
#' @param pool An \code{"eq_pool"} object from \code{\link{eq_pool}}.
#' @param S0 Current (spot) price of the underlying.
#' @param T Maturity, in years.
#' @param K Numeric vector of strikes.
#' @param type \code{"call"} or \code{"put"}.
#' @param eps_scale Optional scalar multiplier applied to the bootstrapped
#'   residuals before reconstruction (default \code{1}); mainly a tool for
#'   sensitivity analysis, not needed in ordinary use.
#'
#' @return A list with components \code{prices} (named by strike), \code{c_ds}
#'   (the Duan-Simonato correction applied), and \code{sigma_atm} (the
#'   annualized at-the-money volatility implied by the simulated terminal
#'   distribution).
#'
#' @references
#' Duan, J.-C. and Simonato, J.-G. (1998). Empirical martingale simulation
#' for asset prices. \emph{Management Science}, 44(9), 1218-1233.
#'
#' @examples
#' set.seed(1)
#' prices <- cumprod(1 + rnorm(500, 0.0003, 0.01)) * 100
#' pool <- eq_pool(prices, r = 0.03)
#' eq_price(pool, S0 = 100, T = 0.5, K = c(95, 100, 105), type = "call")
#'
#' @export
eq_price <- function(pool, S0, T, K, type = c("call", "put"), eps_scale = 1) {
  type <- match.arg(type)
  stopifnot(inherits(pool, "eq_pool"))

  resids <- pool$resampled_resids * eps_scale
  step_idx <- max(1, min(round(T / pool$dt), pool$n_w))
  cum_resids <- colSums(resids[seq_len(step_idx), , drop = FALSE])

  logD_T <- log(S0) + cum_resids
  S_T_raw <- exp(logD_T + pool$r * T)
  disc_T <- exp(-pool$r * T)
  c_ds <- S0 / (disc_T * mean(S_T_raw))
  S_T <- S_T_raw * c_ds

  payoff <- if (type == "call") {
    function(k) pmax(S_T - k, 0)
  } else {
    function(k) pmax(k - S_T, 0)
  }
  prices <- vapply(K, function(k) disc_T * mean(payoff(k)), numeric(1))
  names(prices) <- K

  list(
    prices = prices,
    c_ds = c_ds,
    sigma_atm = stats::sd(log(S_T) - log(S0)) / sqrt(T)
  )
}
