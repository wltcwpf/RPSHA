#' A function for event selection
#'
#' This function takes the target hazard (and target hazard magnitude and/or distance distributions)
#' and the hazard (and hazard magnitude and/or distance distributions) produced by each of the
#' considered events to efficiently and optimally select a subset events with adjusted annual
#' occurrence rates by using LASSO regression.
#'
#' @param Y The target hazard and target hazard magnitude/distance distributions
#' @param X The hazard matrix produced by each of the considered events
#' @param target_num The number of final selected event subset, a positive integer
#' @param max_rate_multiplier The maximum value for the multiplier of annual occurrence rate,
#' a positive numeric number. It is only working when \code{target_num} is not set.
#' @param min_hazard The minimum hazard value in \code{Y} that is considered in regression.
#' This is used to avoid having extremely large weights. Note the weights is 1/\code{Y}.
#' The default value for \code{min_hazard} is 1e-7.
#'
#' @return A list of three elements is returned. The first element is
#' the corresponding adjusted rate multipliers for each event
#' (or each column in \code{X}), the second element is recovered hazard by the selected subset events,
#' and the last element is the logarithmic mean squared error.
#' The events with positive factors are selected events and the events with factors = 0 are
#' NOT selected.
#'
#' @importFrom glmnet glmnet
#' @importFrom stats median
#'
#' @export
event_selection <- function(Y, X, target_num = NULL, max_rate_multiplier = 1000, min_hazard = 1e-7) {

  Weight <- 1 / Y

  if (is.numeric(min_hazard) & (min_hazard >= 0)) {
    Idx_include <- which(Y >= min_hazard)  # only consider the hazards with more than min_hazard
  } else {
    Idx_include <- seq(1, length(Y))
  }

  if (is.null(target_num) & is.null(max_rate_multiplier))
    stop('Please input a valid number for either target_num or max_rate_multiplier!')

  if (!is.null(target_num) & (target_num > 0)) {

    target_num <- round(target_num)
    # use number of final selected events to constrain
    lasso_res <- glmnet(x = X[Idx_include, ], y = Y[Idx_include],
                        weights = Weight[Idx_include], pmax = ncol(X),
                        lambda.min.ratio = 0, nlambda = 10000,
                        lower.limits = 0, intercept = F,
                        dfmax = target_num)

    betas <- lasso_res$beta[, ncol(lasso_res$beta)]

    # correct the penalized coefficients
    re_correction_factor <- log(Y[Idx_include]) - log(X[Idx_include, ] %*% betas)

    re_correction_factor <- exp(mean(re_correction_factor[re_correction_factor != Inf]))

    betas <- betas * re_correction_factor

    m_sq_logerr <- mean((log(Y[Idx_include]) - log(X[Idx_include, ] %*% betas))^2)

    m_abs_log_rediff <- mean(abs((log(Y[Idx_include]) - log(X[Idx_include, ] %*% betas))/log(Y[Idx_include])))

    m_abs_rediff <- mean(abs((Y[Idx_include] - X[Idx_include, ] %*% betas)/Y[Idx_include]))

    res <- list()

    res$betas <- betas

    res$Y_hat <- X %*% betas

    res$m_sq_logerr <- m_sq_logerr

    res$m_abs_log_rediff <- m_abs_log_rediff

    res$m_abs_rediff <- m_abs_rediff

    return(res)

  } else if (!is.null(max_rate_multiplier) & (max_rate_multiplier > 0)) {

    # use the maximum rate multipliers to constrain
    lasso_res <- glmnet(x = X[Idx_include, ], y = Y[Idx_include],
                        weights = Weight[Idx_include], pmax = ncol(X),
                        lambda.min.ratio = 0, nlambda = 10000,
                        lower.limits = 0, intercept = F,
                        upper.limits = max_rate_multiplier)

    betas <- lasso_res$beta[, ncol(lasso_res$beta)]

    res <- list()

    m_sq_logerr <- mean((log(Y[Idx_include]) - log(X[Idx_include, ] %*% betas))^2)

    m_abs_log_rediff <- mean(abs((log(Y[Idx_include]) - log(X[Idx_include, ] %*% betas))/log(Y[Idx_include])))

    m_abs_rediff <- mean(abs((Y[Idx_include] - X[Idx_include, ] %*% betas)/Y[Idx_include]))

    res$betas <- betas

    res$Y_hat <- X %*% betas

    res$m_sq_logerr <- m_sq_logerr

    res$m_abs_log_rediff <- m_abs_log_rediff

    res$m_abs_rediff <- m_abs_rediff

    return(res)

  } else {
    stop('Please input a valid number for either target_num or max_rate_multiplier!')
  }
}



