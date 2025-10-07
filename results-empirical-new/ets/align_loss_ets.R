library(glmnet)
library(here)

source(here("R/cv_icomb.R"))
source(here("R/lambda_path.R"))
load(here("data/tourism_base_forecasts_ets.Rdata"))

train_size_cv <- 80 # train_size * 2 / 3
k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# unstandardized variations
out <- cv_icomb(fitted_values[k, ,], original_y[k, ,],
                forecasts = fcasts_base[k, ,], 
                standardize = FALSE, 
                standardize_response = FALSE,
                intercept = TRUE,
                alpha = 0, 
                train_size = train_size_cv,
                single_lambda = FALSE,
                align_loss = TRUE,
                nlambda = 200,
                lambda_min_ratio = "expand")
fcasts_icomb_ridge_unstd_intr <- out$preds # predictions for lambda_min and lambda_1se
info_icomb_ridge_unstd_intr <- out$info # information about chosen lambdas and mses


out <- cv_icomb(fitted_values[k, ,], original_y[k, ,],
                forecasts = fcasts_base[k, ,], 
                standardize = FALSE, 
                standardize_response = FALSE,
                intercept = FALSE,
                alpha = 0, 
                train_size = train_size_cv,
                single_lambda = FALSE,
                align_loss = TRUE,
                nlambda = 200,
                lambda_min_ratio = "expand")
fcasts_icomb_ridge_unstd_nointr <- out$preds
info_icomb_ridge_unstd_nointr <- out$info


out <- cv_icomb(fitted_values[k, ,], original_y[k, ,],
                forecasts = fcasts_base[k, ,], 
                standardize = FALSE, 
                standardize_response = FALSE,
                intercept = TRUE,
                alpha = 1, 
                train_size = train_size_cv,
                single_lambda = FALSE,
                align_loss = TRUE,
                nlambda = 200,
                lambda_min_ratio = "expand")
fcasts_icomb_lasso_unstd_intr <- out$preds
info_icomb_lasso_unstd_intr <- out$info


out <- cv_icomb(fitted_values[k, ,], original_y[k, ,],
                forecasts = fcasts_base[k, ,], 
                standardize = FALSE, 
                standardize_response = FALSE,
                intercept = FALSE,
                alpha = 1, 
                train_size = train_size_cv,
                single_lambda = FALSE,
                align_loss = TRUE,
                nlambda = 200,
                lambda_min_ratio = "expand")
fcasts_icomb_lasso_unstd_nointr <- out$preds
info_icomb_lasso_unstd_nointr <- out$info


# standardized variations
out <- cv_icomb(fitted_values[k, ,], original_y[k, ,],
                forecasts = fcasts_base[k, ,], 
                standardize = TRUE, 
                standardize_response = TRUE,
                intercept = TRUE,
                alpha = 0, 
                train_size = train_size_cv,
                single_lambda = FALSE,
                align_loss = TRUE,
                nlambda = 200,
                lambda_min_ratio = "expand")
fcasts_icomb_ridge_std_intr <- out$preds 
info_icomb_ridge_std_intr <- out$info 


out <- cv_icomb(fitted_values[k, ,], original_y[k, ,],
                 forecasts = fcasts_base[k, ,], 
                 standardize = TRUE, 
                 standardize_response = TRUE,
                 intercept = FALSE,
                 alpha = 0, 
                 train_size = train_size_cv,
                 single_lambda = FALSE,
                 align_loss = TRUE,
                 nlambda = 200,
                 lambda_min_ratio = "expand")
fcasts_icomb_ridge_std_nointr <- out$preds
info_icomb_ridge_std_nointr <- out$info


out <- cv_icomb(fitted_values[k, ,], original_y[k, ,],
                forecasts = fcasts_base[k, ,], 
                standardize = TRUE, 
                standardize_response = TRUE,
                intercept = TRUE,
                alpha = 1, 
                train_size = train_size_cv,
                single_lambda = FALSE,
                align_loss = TRUE,
                nlambda = 200,
                lambda_min_ratio = "expand")
fcasts_icomb_lasso_std_intr <- out$preds
info_icomb_lasso_std_intr <- out$info


out <- cv_icomb(fitted_values[k, ,], original_y[k, ,],
                forecasts = fcasts_base[k, ,], 
                standardize = TRUE, 
                standardize_response = TRUE,
                intercept = FALSE,
                alpha = 1, 
                train_size = train_size_cv,
                single_lambda = FALSE,
                align_loss = TRUE,
                nlambda = 200,
                lambda_min_ratio = "expand")
fcasts_icomb_lasso_std_nointr <- out$preds
info_icomb_lasso_std_nointr <- out$info

# standardize Xs but not Ys
out <- cv_icomb(fitted_values[k, ,], original_y[k, ,],
                forecasts = fcasts_base[k, ,], 
                standardize = TRUE, 
                standardize_response = FALSE,
                intercept = TRUE,
                alpha = 0, 
                train_size = train_size_cv,
                single_lambda = FALSE,
                align_loss = TRUE,
                nlambda = 200,
                lambda_min_ratio = "expand")
fcasts_icomb_ridge_stdx_intr <- out$preds 
info_icomb_ridge_stdx_intr <- out$info 


out <- cv_icomb(fitted_values[k, ,], original_y[k, ,],
                forecasts = fcasts_base[k, ,], 
                standardize = TRUE, 
                standardize_response = FALSE,
                intercept = FALSE,
                alpha = 0, 
                train_size = train_size_cv,
                single_lambda = FALSE,
                align_loss = TRUE,
                nlambda = 200,
                lambda_min_ratio = "expand")
fcasts_icomb_ridge_stdx_nointr <- out$preds
info_icomb_ridge_stdx_nointr <- out$info


out <- cv_icomb(fitted_values[k, ,], original_y[k, ,],
                forecasts = fcasts_base[k, ,], 
                standardize = TRUE, 
                standardize_response = FALSE,
                intercept = TRUE,
                alpha = 1, 
                train_size = train_size_cv,
                single_lambda = FALSE,
                align_loss = TRUE,
                nlambda = 200,
                lambda_min_ratio = "expand")
fcasts_icomb_lasso_stdx_intr <- out$preds
info_icomb_lasso_stdx_intr <- out$info


out <- cv_icomb(fitted_values[k, ,], original_y[k, ,],
                forecasts = fcasts_base[k, ,], 
                standardize = TRUE, 
                standardize_response = FALSE,
                intercept = FALSE,
                alpha = 1, 
                train_size = train_size_cv,
                single_lambda = FALSE,
                align_loss = TRUE,
                nlambda = 200,
                lambda_min_ratio = "expand")
fcasts_icomb_lasso_stdx_nointr <- out$preds
info_icomb_lasso_stdx_nointr <- out$info

save(list = ls(all = TRUE), file = paste0("results/tourism_icomb_reconciled_forecasts_ets_", k, ".Rdata"))