#' Price arithmetic Asian call options from an empirical residual pool
#'
#' Extends \code{\link{eq_price}} to a path-dependent payoff by applying a
#' date-specific Duan-Simonato correction at each fixing date, then averaging
#' the corrected simulated path before applying the payoff. Requires no
#' additional inputs beyond the pool already built for European pricing.
#'
#' @param pool An \code{"eq_pool"} object from \code{\link{eq_pool}}.
#' @param S0 Current (spot) price of the underlying.
#' @param T Maturity, in years.
#' @param K Numeric vector of strikes.
#' @param n_fix Number of (equally spaced) fixing/averaging dates up to \code{T}.
#'
#' @return Named numeric vector of arithmetic Asian call prices, one per strike.
#'
#' @export
eq_price_asian <- function(pool, S0, T, K, n_fix) {
  stopifnot(inherits(pool, "eq_pool"))

  fixing_times <- seq(T / n_fix, T, length.out = n_fix)
  fixing_idx <- unique(pmin(pool$n_w, pmax(1, round(fixing_times / pool$dt))))
  max_idx <- max(fixing_idx)

  cum_path <- apply(pool$resampled_resids[seq_len(max_idx), , drop = FALSE], 2, cumsum)

  S_fix <- sapply(fixing_idx, function(idx) {
    t_k <- idx * pool$dt
    cum_resids <- cum_path[idx, ]
    S_raw <- exp(log(S0) + cum_resids + pool$r * t_k)
    disc_k <- exp(-pool$r * t_k)
    c_k <- S0 / (disc_k * mean(S_raw))
    S_raw * c_k
  })

  avgA <- rowMeans(S_fix)
  disc_T <- exp(-pool$r * T)
  prices <- vapply(K, function(k) disc_T * mean(pmax(avgA - k, 0)), numeric(1))
  names(prices) <- K
  prices
}
