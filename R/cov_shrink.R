# Shrinking the sample covariance matrix towards a digonal matrix

# Shrunk covariance matrix - Schafer and strimmer approach
# Shrinkage intensity parameter is based on the correlation matrix
# Then DRD is used to get back the covariance matrix

shrink_estim <- function (x, tar) 
{
  if (is.matrix(x) == TRUE && is.numeric(x) == FALSE) 
    stop("The data matrix must be numeric!")
  p <- ncol(x)
  n <- nrow(x)
  covm <- crossprod(x) / n
  corm <- cov2cor(covm)
  xs <- scale(x, center = FALSE, scale = sqrt(diag(covm)))
  v <- (1/(n * (n - 1))) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
  diag(v) <- 0
  corapn <- cov2cor(tar)
  d <- (corm - corapn)^2
  lambda <- sum(v)/sum(d)
  lambda <- max(min(lambda, 1), 0)
  shrink_cov <- lambda * tar + (1 - lambda) * covm
  return(list(shrink_cov, c("The shrinkage intensity lambda is:", 
                            round(lambda, digits = 4))))
}
