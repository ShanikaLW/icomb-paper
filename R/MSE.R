# A function to calculate MSE for each level 
# Arguments
# x: an hts or gts object giving the structure of the hierarchy
# fcasts: a 3D array of forecasts (n_iter, h, series)
# actual: a 3D array of actual values (n_iter, h, series)
# row_names: row names of the resulting 2D matrix

library(hts)

CalcMSE <- function(fcasts, actual, x, row_names) {
  if (is.hts(x)) {
    nodes <- hts:::Mnodes(x$nodes)
  } else if (is.gts(x)) {
    nodes <- hts:::Mlevel(x$groups)
  } else {
    stop("Not a valid choice for x")
  }
  
  mse <- matrix(, nrow = length(nodes), ncol = h)
  res <- actual[, 1:h, ] - fcasts
  
  cs <- c(0, cumsum(nodes))
  mse[1L, ] <- colMeans(res[, , 1L]^2, na.rm = TRUE)
  for (i in 2:length(nodes)) {
    end <- cs[i + 1L]
    start <- cs[i] + 1L
    series <- seq(start, end)
    mse[i, ] <- rowMeans(colMeans(res[, , series]^2, na.rm = TRUE))
  }
  mse <- rbind(mse, rowMeans(colMeans(res^2, na.rm = TRUE)))
  mse <- cbind(mse, rowMeans(mse[, 1:(dim(fcasts)[2]/2)]), rowMeans(mse))
  colnames(mse) <- c(paste0("h = ", 1:(dim(fcasts)[2])),
                     paste0("h = 1:", dim(fcasts)[2]/2),
                     "Overall")
  rownames(mse) <- row_names
  mse
}

