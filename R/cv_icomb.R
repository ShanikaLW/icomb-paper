# a function to calculate reconciled forecasts using information combination by using rolling origin
# Arguments:
# fitted: in-sample fitted values of all the series in the structure
# actual: in-sample actual values of all the series in the structure
# forecasts: h-steps-ahead forecasts of all the series in the structure
# train_size: number of training observations
# alpha: the elasticnet mixing parameter with 0 <= alpha <= 1 (alpha = 0: ridge) and (alpha = 1: lasso)
# standardize: whether the columns of x should be standardized
# lambda: user supplied lambda sequence
# standardize_response: whether the columns of y should be standardized
# lambda_min_ratio: smallest value of lambda, as a fraction of lambda.max, "expand" uses 10^(-floor(log10(lambda_max))-2).
# nlambda: the number of lambda values
# maxit: Maximum number of passes over the data for all lambda values
# single_lambda: Once the best lambda is identified should we fit the final model
# using the lambda values leading up to the best lambda value (because glmnet uses warm start) or
# using only the best lambda value
# align_loss: whether to align the cross-validation loss function with 
# the final loss function. If TRUE, CV uses MSE on the original scale 
# when standardize = TRUE and standardize_response = TRUE. 

cv_icomb <- function (fitted, actual, forecasts, train_size, alpha = 1, standardize = FALSE, 
                      standardize_response = FALSE, intercept = TRUE, lambda = NULL,
                      lambda_min_ratio = ifelse(nobs < nvars, 0.01, 1e-04), 
                      nlambda = 100, maxit = 1e+07, single_lambda = TRUE,
                      align_loss = TRUE) {
  dimx = dim(fitted)
  if (is.null(dimx) | (dimx[2] <= 1)) 
    stop("fitted should be a matrix with 2 or more columns")
  if (any(is.na(fitted))) 
    stop("fitted has missing values")
  nobs <- dimx[1]
  nvars <- dimx[2]
  
  dimy <- dim(actual)
  if (is.null(dimy) | (dimy[2] <= 1)) 
    stop("actual should be a matrix with 2 or more columns")
  if (any(is.na(actual))) 
    stop("actual has missing values")
  
  if (dimy[1] != dimx[1]) 
    stop("number of observations in actual and fitted does not match")
  
  if (train_size > dimx[1]) 
    stop("number of training observations should be less than that for actual/fitted")
  
  dimfc <- dim(forecasts)
  if (is.null(dimfc) | (dimfc[2] <= 1)) 
    stop("forecasts should be a matrix with 2 or more columns")
  
  if (dimx[2] != dimfc[2]) 
    stop("number of variables in fitted and forecasts does not match")
  
  if (lambda_min_ratio >= 1 & is.numeric(lambda_min_ratio)) 
    stop("lambda_min_ratio should be less than 1")
  
  if (alpha > 1) {
    warning("alpha > 1; set to 1")
    alpha <- 1
  }
  
  if (alpha < 0) {
    warning("alpha < 0; set to 0")
    alpha <- 0
  }

  if (is.null(lambda)) 
    lambda <- lambda_path(fitted, actual, alpha, standardize, standardize_response, intercept, lambda_min_ratio, nlambda)

  niter <- nobs - train_size
  recon <- array(, dim = c(niter, dimy[2], length(lambda)))
  test <- actual[(train_size + 1):nobs, ]
  ysd_all <- array(, dim = c(niter, dimy[2]))
  out <-
  foreach (i=1:niter, .inorder=FALSE, .packages="glmnet") %dopar% {
    train_set <- 1:(train_size + i - 1)
    xdata <- fitted[train_set, ]
    ydata <- actual[train_set, ]
    
    xsd <- sqrt(colMeans(scale(xdata, center = TRUE, scale = FALSE)^2))
    xconst_var <- xsd < sqrt(.Machine$double.eps) # constant predictors
    
    ysd <- ysd_all[i, ] <- sqrt(colMeans(scale(ydata, center = TRUE, scale = FALSE)^2))
    yconst_var <- ysd < sqrt(.Machine$double.eps) # constant responses
    
    if (any(yconst_var)) {
      warning("There are constant responses.")
      
      if (standardize_response)
        ydata <- ydata[, !yconst_var]
    }
    
    fit <- glmnet(xdata[, !xconst_var], ydata, family = "mgaussian", standardize = standardize,
                  standardize.response = standardize_response, intercept = intercept, alpha = alpha, lambda = lambda, maxit = maxit)
    list(
      recon = predict(fit, newx = fitted[train_size + i, !xconst_var, drop = FALSE]),
      coef = coef(fit) ## currently unused, but can be also extracted
      )
  }
  ## need to extract recon from the results
  recon <- sapply(out, function(x) x$recon[1,,])
  ## sapply will flatten the result dims so have to restore them
  dim(recon) <- c(dim(out[[1]]$recon)[-1], dim(recon)[2])
  ## make the iter first dim instead of last
  recon <- aperm(recon, c(3,1,2))
  err <- sweep(recon, 1:2, test)

  if (align_loss) {
    mse <- colMeans(colMeans(err^2))
    mse_over_response <- apply(err^2, c(1, 3), mean)
  } else {
    err_scaled <- sweep(err, MARGIN = 1:2, STATS = ysd_all, FUN = "/")
    mse <- colMeans(colMeans(err_scaled^2))
    mse_over_response <- apply(err_scaled^2, c(1, 3), mean)
  }
  
  idx_min <- which.min(mse)
  lambda_min <- lambda[idx_min]
  
  se_mse <- apply(mse_over_response, 2, sd) / sqrt(niter)
  lambda_1se <- max(lambda[mse < mse[idx_min] + se_mse[idx_min]])
  idx_1se = which(lambda == lambda_1se)
  lambda_info <- c(lambda_max = max(lambda), lambda_min = lambda_min, lambda_min_idx = idx_min,
                   lambda_1se = lambda_1se, lambda_1se_idx = idx_1se)

  xsd <- sqrt(colMeans(scale(fitted, center = TRUE, scale = FALSE)^2))
  xconst_var <- xsd < sqrt(.Machine$double.eps) # constant predictors
  
  ysd <- sqrt(colMeans(scale(actual, center = TRUE, scale = FALSE)^2))
  yconst_var <- ysd < sqrt(.Machine$double.eps) # constant responses
  
  if (any(yconst_var)) {
    warning("There are constant responses.")
    
    if (standardize_response)
      actual <- actual[, !yconst_var]
  }
  
  if (single_lambda) {
    fit <- glmnet(fitted[, !xconst_var], actual, family = "mgaussian", standardize = standardize, 
                  standardize.response = standardize_response, intercept = intercept, alpha = alpha, 
                  lambda = c(lambda_min, lambda_1se), maxit = maxit)
    list(preds = predict(fit, newx = forecasts[, !xconst_var]), 
         info = list(lambda_info = lambda_info, mse_info = mse, nnzeros = fit$dfmat[1, 1]))
  } else {
    lambda_subset <- lambda[lambda >= lambda_min]
    fit <- glmnet(fitted[, !xconst_var], actual, family = "mgaussian", standardize = standardize, 
                  standardize.response = standardize_response, intercept = intercept, alpha = alpha, 
                  lambda = lambda_subset, maxit = maxit)
    list(preds = predict(fit, newx = forecasts[, !xconst_var], s = c(lambda_min, lambda_1se), exact = TRUE), 
         info = list(lambda_info = lambda_info, mse_info = mse, nnzeros = fit$dfmat[1, idx_min]))
  }
}

# Test cases
# set.seed(2024)
# nodes <- list(2, c(2, 2))
# bts <- matrix(, nrow = 100, ncol = 4)
# for (i in 1:4) {
#   bts[, i] <- arima.sim(list(order = c(1,0,0), ar = 0.7), n = 100)
# }
# eg_hts <- hts(bts, nodes)
# ally <- allts(eg_hts)
# fitted <- matrix(, nrow = 100, ncol = 7)
# fcast <- matrix(, nrow = 1, ncol = 7)
# 
# for (i in 1:7) {
#   fit <- auto.arima(ally[, i])
#   fcast[1, i] <- forecast(fit, h = 1)$mean
#   fitted[, i] <- fitted(fit)
# }
# 
# out <- cv_icomb(fitted, ally, fcast, train_size = 70)[[1]]
# utmat <- cbind(diag(3), -smatrix(eg_hts)[1:3, ])
# utmat %*% out[,,1]
