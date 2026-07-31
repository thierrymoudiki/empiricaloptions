# empiricaloptions

Empirical martingale option pricing from a single historical price series —
no option-market calibration required.

Companion package to *Semi-parametric option pricing based on underlying's
historical data* (Moudiki, 2026). Given a raw price history and a discount
rate, `empiricaloptions` builds a residual pool (AR(1) filter + stationary block
bootstrap, Politis and Romano 1994) and prices European or arithmetic Asian
options against it, with a Duan-Simonato (1998) scalar correction enforcing
the martingale condition required by the Fundamental Theorem of Asset
Pricing.

## Install

```r
# install.packages("remotes")
remotes::install_github("Techtonique/empiricaloptions")
```

## Usage

```r
library(empiricaloptions)

set.seed(1)
prices <- cumprod(1 + rnorm(500, 0.0003, 0.01)) * 100

pool <- eq_pool(prices, r = 0.03)

eq_price(pool, S0 = tail(prices, 1), T = 0.5, K = c(90, 100, 110), type = "call")
eq_price(pool, S0 = tail(prices, 1), T = 0.5, K = c(90, 100, 110), type = "put")
eq_price_asian(pool, S0 = tail(prices, 1), T = 0.5, K = c(90, 100, 110), n_fix = 6)
```

Build the pool once per discount rate and reuse it across strikes,
maturities, and payoff types — the expensive step is the stationary
bootstrap draw, not the pricing itself.

## Status

Skeleton package matching the methodology and notation of the companion
paper. Not yet on CRAN.

## License

MIT
