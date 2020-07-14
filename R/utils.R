#' Wald-Style Confidence Intervals
#'
#' @importFrom stats qnorm
#' @param est mean_hat
#' @param se std_hat
#' @param level p
#' @param q q
#' @keywords internal
#' @export
wald_ci <- function(est, se, level = 0.95, q = NULL) {
  if (is.null(q)) {
    q <- abs(stats::qnorm(p = (1 - level) / 2))
  }

  ci_low <- est - q * se
  ci_high <- est + q * se
  return(cbind(ci_low, ci_high))
}

#' Probability Clipping
#' 
#' @export
prob_clip = function(ps, clip=.05) {
  if (!all(ps<=1 & ps>=0)) {
    stop("All ps should within [0, 1].")
  }
  
  ps_clip <- c()
  for (p in ps) {
    if (1-p < clip) ps_clip <- c(ps_clip, 1-clip)
    else if (p < clip) ps_clip <- c(ps_clip, clip)
    else ps_clip <- c(ps_clip, p)
  }
  return(ps_clip)
}

#' L2 Loss
#' 
#' @export
l2_loss = function(x, y) {
  return(sqrt(sum((x - y)^2)) / length(x))
}

#' Expit
#' 
#' @param x x
#' @export
expit <- function(x) {
  return(1/(1+exp(-x)))
}

#' Logit
#' 
#' @param x x
#' @export
logit <- function(x) {
  return(log(x/(1-x)))
}


#' Standard Error: Mean of Logistic Regression
#' 
#' @param x covariates
#' @param y outcomes
#' @param cov covariance matrix of coefficients
#' @export
deltaMeanLogistic <- function(x, y, cov) {
  if (ncol(x)==(nrow(cov)-1)) {
    x <- cbind(rep(1, length(y)), x)
  }
  if (!all(y<=1 & y>=0)) {
    stop("All ys should within [0, 1].")
  }
  if (!(length(y)==nrow(x)) | !(ncol(x)==nrow(cov))) {
    stop("Input dimensions do not match.")
  }
  
  deriv <- colMeans(y*(1-y) * x)
  return(deriv%*%cov%*%deriv)
}

#' Standard Error: Mean of OLS Regression
#' 
#' @param x covariates
#' @param y outcomes
#' @param cov covariance matrix of coefficients
#' @export
deltaMeanOLS <- function(x, y, cov) {
  if (ncol(x)==(nrow(cov)-1)) {
    x <- cbind(rep(1, length(y)), x)
  }
  if (!(length(y)==nrow(x)) | !(ncol(x)==nrow(cov))) {
    stop("Input dimensions do not match.")
  }
  
  deriv <- colMeans(x)
  return(deriv%*%cov%*%deriv)
}

