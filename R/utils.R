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

