#' A function for scenario (event or ground motion map) selection
#'
#' This function takes the target hazard (and target hazard magnitude and/or distance distributions)
#' and the hazard (and hazard magnitude and/or distance distributions) produced by each of the
#' considered event scenarios or ground motion map scenarios to efficiently and
#' optimally select a series of subset scenarios with adjusted annual
#' occurrence rates by using LASSO regression.
#'
#' @param Y The target hazard and target hazard magnitude/distance distributions
#' @param X The hazard matrix produced by each of the considered scenarios
#' @param min_hazard The minimum hazard value in \code{Y} that is considered in regression.
#' This is used to avoid having extremely large weights. Note the weights is 1/\code{Y}.
#' The default value for \code{min_hazard} is 1e-7.
#' @param output_dir If not NULL, the results will be written into a few CSV files under the specified directory
#' @param max_rate_multiplier The maximum value for the multiplier of annual occurrence rate when running LASSO,
#' a positive numeric number. It is an additional constraint for scenario selection.
#' The rate multipliers are not straightly constrained: the rate multipliers could be much larger
#' when the number of selected scenarios is too small; and the rate multipliers will get better
#' constrained when the number of selected scenarios gets more.
#' @param num_lambda The maximum number of trail penalty parameters. A larger \code{num_lambda}
#' will search for more possible subset events.
#' @param Weight The weight vector for each element in \code{Y}. The default is 1 that indicates
#' using the default weight of scenario selection algorithm (which is 1/\code{Y}). If the other
#' weights are preferred, then a vector with the same length as \code{Y} for the desired weights
#' can be specified. Since the weight of 1/\code{Y} has been set in the function, so you'd need
#' to multiply your desired weights by \code{Y} to compensate the default weights.
#' @param min_rate_adj The minimum rate multiplier (or adjustments) for the selected events.
#' The default is 1, which means that the selected events with rate adjustments less than
#' \code{min_rate_adj} will be removed and the rate adjustments of the remaining selected
#' events will be updated. Note \code{min_rate_adj} needs to be non-negative value.
#'
#' @return A list of three elements is returned. The first element is
#' the corresponding adjusted rate multipliers for each event
#' (or each column in \code{X}), the second element is recovered hazard by the selected subset events,
#' and the last element is the logarithmic mean squared error.
#' The events with positive factors are selected events and the events with factors = 0 are
#' NOT selected.
#' \item{betas}{A matrix of the adjusted rate multipliers for each scenario for each set of selected
#' scenarios. The dimension is total number of candidate scenarios by the number of set of selected
#' scenarios. The scenarios with zero multipliers are unselected}
#' \item{errors_mat}{A matrix with four columns: the number of selected scenarios following by
#' three error metrics (mean of squared logarithmic residuals, mean of absolute logarithmic
#' relative errors, and mean of absolute arithmetic relative errors)}
#' \item{Y_hat}{The predicted hazard by each set of selected scenarios}
#'
#' @importFrom glmnet glmnet
#' @importFrom stats median
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @export
scenario_selection <- function(Y, X, min_hazard = 1e-7, output_dir = NULL, max_rate_multiplier = NULL,
                               num_lambda = 10000, Weight = 1, min_rate_adj = 1) {

  Weight <- (1 / Y) * Weight

  if (is.numeric(min_hazard) & (min_hazard >= 0)) {
    Idx_include <- which(Y >= min_hazard)  # only consider the hazards with more than min_hazard
  } else {
    Idx_include <- seq(1, length(Y))
  }

  # run LASSO
  if (is.null(max_rate_multiplier)) {
    cat('Conduct selection \n')
    lasso_res <- glmnet(x = X[Idx_include, ], y = Y[Idx_include],
                        weights = Weight[Idx_include],
                        lambda.min.ratio = 0, nlambda = num_lambda,
                        lower.limits = 0, intercept = F, trace.it = TRUE)

    lasso_res <- relax.glmnet(fit = lasso_res,
                              x = X[Idx_include, ], y = Y[Idx_include],
                              weights = Weight[Idx_include],
                              lambda.min.ratio = 0, nlambda = num_lambda,
                              lower.limits = 0, intercept = F,
                              trace.it = FALSE)
  } else if (max_rate_multiplier > 0) {
    cat('Conduct selection \n')
    lasso_res <- glmnet(x = X[Idx_include, ], y = Y[Idx_include],
                        weights = Weight[Idx_include],
                        lambda.min.ratio = 0, nlambda = num_lambda,
                        lower.limits = 0, intercept = F, trace.it = TRUE,
                        upper.limits = max_rate_multiplier)

    lasso_res <- relax.glmnet(fit = lasso_res,
                              x = X[Idx_include, ], y = Y[Idx_include],
                              weights = Weight[Idx_include],
                              lambda.min.ratio = 0, nlambda = num_lambda,
                              lower.limits = 0, intercept = F,
                              trace.it = FALSE,
                              upper.limits = max_rate_multiplier)
  } else {
    stop('max_rate_multiplier is not valid, please correct the input!')
  }

  # get the series of selected event sets
  if (min_rate_adj < 0) stop('min_rate_adj is not valid, please correct it!')
  num_increasing <- apply(lasso_res$relaxed$beta, 2, function(x) sum(x > min_rate_adj)) # only keep the variables that have multiplier > min_rate_adj
  uni_num <- sort(unique(num_increasing[num_increasing > 0]))
  idx <- as.numeric(sapply(uni_num, function(x) {
    # find the index such that the corresponding dev.ratio is max.
    tmp_idx <- which(num_increasing == x)
    tmp_idx[which.max(lasso_res$relaxed$dev.ratio[tmp_idx])]
  }))

  m_sq_logerr <- rep(0, length(idx))
  m_abs_log_rediff <- rep(0, length(idx))
  m_abs_rediff <- rep(0, length(idx))
  Y_hat <- matrix(data = 0, nrow = length(Y), ncol = length(idx))
  betas_mat <- matrix(data = 0, nrow = ncol(X), ncol = length(idx))

  cat('Re-correct and output results \n')
  pb <- txtProgressBar(min = 0, max = length(idx), initial = 0,
                       style = 3)

  for (i in 1:length(idx)) {
    betas <- lasso_res$relaxed$beta[, idx[i]]
    betas <- ifelse(betas > min_rate_adj, betas, 0)

    # correct the penalized coefficients
    re_correction_factor <- log(Y[Idx_include]) - log(X[Idx_include, ] %*% betas)
    re_correction_factor <- exp(mean(re_correction_factor[re_correction_factor != Inf]))
    if (re_correction_factor > 1)
      betas <- betas * re_correction_factor
    betas_mat[, i] <- betas

    m_sq_logerr[i] <- mean((log(Y[Idx_include]) -
                              ifelse(X[Idx_include, ] %*% betas == 0,
                                     0,
                                     log(X[Idx_include, ] %*% betas)))^2)
    m_abs_log_rediff[i] <- mean(abs((log(Y[Idx_include]) -
                                       ifelse(X[Idx_include, ] %*% betas == 0,
                                              0,
                                              log(X[Idx_include, ] %*% betas)))/log(Y[Idx_include])))
    m_abs_rediff[i] <- mean(abs((Y[Idx_include] - X[Idx_include, ] %*% betas)/Y[Idx_include]))
    Y_hat[, i] <- X %*% betas

    setTxtProgressBar(pb, i)
  }

  colnames(betas_mat) <- paste0('selected_scenario_num_', uni_num)
  rownames(betas_mat) <- paste0('rate_multiplier_scenario_', 1:ncol(X))
  errors_mat <- cbind(uni_num, m_sq_logerr, m_abs_log_rediff, m_abs_rediff)
  colnames(errors_mat) <- c('number_selected_scenario', 'mean_squared_log_error',
                            'mean_absolute_log_relative_difference',
                            'mean_absolute_relative_difference')
  colnames(Y_hat) <- paste0('selected_scenario_num_', uni_num)

  if (!is.null(output_dir)) {
    if (dir.exists(output_dir)) {
      write.csv(betas_mat, file = paste0(output_dir, '/betas_mat.csv'))
      write.csv(errors_mat, file = paste0(output_dir, '/errors_mat.csv'), row.names = F)
      write.csv(Y_hat, file = paste0(output_dir, '/Y_hat.csv'), row.names = F)
      cat('\n The results are exported! \n')
    } else {
      cat('\n The results are not written into any files, please input a valid the output_dir! \n')
    }
  } else {
    cat('\n The results are not written into any files. \n')
  }

  res <- list()
  res$betas <- betas_mat
  res$errors_mat <- errors_mat
  res$Y_hat <- Y_hat
  return(res)
}


