library(hts)
library(here)
source(here("R/cov_shrink.R"))
load(here("results-empirical-traditional/tourism_base_forecasts_ets.Rdata"))

fcasts_bu <- array(, dim = c(niter, h, nc))
for (i in 1:niter) {
  bu <- fcasts_base[i, , (nc - nb + 1):nc]
  colnames(bu) <- colnames(bts_subset)
  fcasts_bu[i, ,] <- suppressMessages(allts(hts(bu, 
                                                characters = c(1, 1, 1))))
}

# OLS
fcasts_ols <- array(, dim = c(niter, h, nc))
pols <- jmat - jmat %*% t(utmat) %*% solve(utmat %*% t(utmat)) %*% utmat
for (i in 1:niter) {
  fcasts_ols[i, ,] <- as.matrix(t(s %*% pols %*% t(fcasts_base[i, ,])))
}

# WLS
fcasts_wls <- array(, dim = c(niter, h, nc))
for (i in 1:niter) {
  resid <- na.omit(original_y[i, ,] - fitted_values[i, ,])
  diagw <- Diagonal(x = colMeans(resid^2))
  pwls <- jmat - jmat %*% diagw %*% t(utmat) %*% 
    solve(utmat %*% diagw %*% t(utmat)) %*% utmat
  fcasts_wls[i, ,] <- as.matrix(t(s %*% pwls %*% t(fcasts_base[i, ,])))
}

# MinT (Shrink)
fcasts_mint <- array(, dim = c(niter, h, nc))
for (i in 1:niter) {
  resid <- na.omit(original_y[i, ,] - fitted_values[i, ,])
  tar <- Diagonal(x = colMeans(resid^2))
  w1s <- shrink_estim(resid, tar)[[1]]
  pmint <- jmat - jmat %*% w1s %*% t(utmat) %*% 
    solve(utmat %*% w1s %*% t(utmat)) %*% utmat
  fcasts_mint[i, ,] <- as.matrix(t(s %*% pmint %*% t(fcasts_base[i, ,])))
}

# EMinT-U
fcasts_emintu <- array(, dim = c(niter, h, nc))
for (i in 1:niter) {
  fitted_svd <- svd(fitted_values[i, ,])
  positive_diag <- fitted_svd$d[abs(fitted_svd$d) > sqrt(.Machine$double.eps)]
  fcasts_emintu[i, ,] <- as.matrix(fcasts_base[i, ,] %*% fitted_svd$v[, 1:length(positive_diag)] %*%
                                     Diagonal(x = 1 / positive_diag) %*% 
                                     t(fitted_svd$u[, 1:length(positive_diag)]) %*% original_y[i, ,])
}

save(list = ls(all = TRUE), file = here("results-empirical-traditional/tourism_traditional_reconciled_forecasts_ets.Rdata"))
