####################################################
# Uses rolling window forecasting approach
# Different ARIMA models are fitted through-out the process
# by minimizing the AICc
# Base forecasts were obtained from 120 onwards
# a pure hierarchy: total > states (7) > zones (27) > regions (77)
######################################################

library(hts)
library(Matrix)
library(forecast)
library(here)

# load the dataset
tourism <- read.csv(here("data/TRA-Final.csv", header = TRUE))

# extract the bottom level series 
bts <- tourism[, 3:ncol(tourism)]
bts <- ts(bts, start = c(1998, 1), frequency = 12)

# for obtaining regional level data
idx_start <- seq(1, ncol(bts), by = 4)
idx_end <- seq(4, ncol(bts), by = 4)

colnames_bts_subset <- substr(colnames(bts), 1, 3)[idx_start]
bts_subset <- matrix(, nrow(bts), ncol(bts)/4)
for (i in 1:length(colnames_bts_subset)) {
  bts_subset[, i] <- rowSums(bts[, idx_start[i]:idx_end[i]])
}

colnames(bts_subset) <- colnames_bts_subset
bts_subset <- ts(bts_subset, start = c(1998, 1), frequency = 12)
tourism_hts <- hts(bts_subset, characters = c(1, 1, 1))

ally <- allts(tourism_hts)
nc <- ncol(ally) # number of series
nr <- nrow(ally) # number of observations for each series
nb <- ncol(bts_subset) # number of bottom level series
h <- 12
nfreq <- 12L
time_attr <- tsp(ally)

# required matrices for MinT
s <- smatrix(tourism_hts)
s <- as(s, "sparseMatrix")
nagg <- nc - nb
cmat <- s[1:nagg, ]
utmat <- cbind(Diagonal(nagg), -1*cmat)
jmat <- sparseMatrix(i = 1L:nb, j = (nagg + 1L):nc, x = rep(1L, nb),
                     dims = c(nb, nc))

train_size <- 120L  # number of observations for training
niter <- nr - train_size  # total number of iterations
test_size <- nr - train_size # the available obs for testing

# placeholders for the necessary information
original_y <- array(, dim = c(niter, train_size, nc))
test_y <- array(, dim = c(niter, test_size, nc))
fcasts_base <- array(, dim = c(niter, h, nc))
fitted_values <- array(, dim = c(niter, train_size, nc))

for (i in 1L:nc) {
  # the whole sample
  single_y <- ts(ally[, i], start = c(1998, 1), frequency = nfreq)    
  
  # start rolling window forecasting
  for (j in 1:test_size) {
    starting <- time_attr[1L] + (j - 1L)/nfreq
    ending <- time_attr[1L] + (train_size + j - 2L)/nfreq
    original_y[j, 1:train_size, i] <- y_train <- window(single_y, start = starting, end = ending)
    test_y[j, 1L:(test_size - j + 1), i] <- window(single_y, start = ending + 1/nfreq, end = time_attr[2])
    
    # fitting ARIMA models
    fit <- auto.arima(y_train)
    fitted_values[j, 1:train_size, i] <- fitted(fit)
    fcasts_base[j, , i] <- forecast(fit, h = h)$mean
  }
  print(paste("done with variable", i))
}

save(list = ls(all = TRUE), file = here("results-empirical-traditional/tourism_base_forecasts_arima.Rdata"))

