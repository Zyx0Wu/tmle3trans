#' Helper Functions for Survival Analysis
#'
#' Handles the W (covariates), S (location), T_tilde (time-to-event),
#' Delta (censoring indicator), t_max (the maximum time to estimate) survival data structure
#'
#' @param data a \code{data.frame}, or \code{data.table} containing data for use in estimation
#' @param node_list a list of character vectors, listing the variables that comprise each node
#' @param variable_types a list of variable types, one for each node. If missing, variable types will be guessed
#' @param tmle_task a \code{\link{tmle3_Task}} as constructed via \code{survival_task}
#' @param learner_list a list of sl3 learners, one for A and one for Y to be used for likelihood estimation
#' @param ... extra arguments.
#' @export
#' @rdname survival_trans
survival_o_npsem <- function(node_list, variable_types = NULL) {
  # make the tmle task
  
  # define censoring (lost to followup node)
  censoring <- define_node("pre_failure", node_list$pre_failure, c())
  npsem <- list(
    # TODO: causal relation, handle t_max
    define_node("W", node_list[["W"]], variable_type = variable_types[["W"]]),
    define_node("S", node_list[["S"]], c("W"), variable_type = variable_types[["S"]]),
    define_node("T.tilde", node_list[["T.tilde"]], c("W"), variable_type = variable_types[["T.tilde"]]),
    define_node("Delta", node_list[["Delta"]], variable_type = variable_types[["Delta"]]),
    censoring,
    # TODO: remove t parent, handle in get_regression
    define_node("failed", node_list[["failed"]], c("W"), variable_type = variable_types[["failed"]], censoring_node=censoring),
    define_node("censored", node_list[["censored"]], c("W"), variable_type = variable_types[["censored"]], censoring_node=censoring)   
  )
  
  return(npsem)
}

#' @export
#' @rdname survival_trans
survival_e_npsem <- function(node_list, variable_types = NULL) {
  # make the tmle task
  
  # define censoring (lost to followup node)
  censoring <- define_node("pre_failure", node_list$pre_failure, c())
  npsem <- list(
    # TODO: causal relation, handle t_max
    define_node("W", node_list[["W"]], variable_type = variable_types[["W"]]),
    define_node("S", node_list[["S"]], c("W"), variable_type = variable_types[["S"]]),
    define_node("A", node_list[["A"]], c("W"), variable_type = variable_types[["A"]]),
    define_node("T.tilde", node_list[["T.tilde"]], c("A", "W"), variable_type = variable_types[["T.tilde"]]),
    define_node("Delta", node_list[["Delta"]], variable_type = variable_types[["Delta"]]),
    censoring,
    # TODO: remove t parent, handle in get_regression
    define_node("failed", node_list[["failed"]], c("A", "W"), variable_type = variable_types[["failed"]], censoring_node=censoring),
    define_node("censored", node_list[["censored"]], c("A", "W"), variable_type = variable_types[["censored"]], censoring_node=censoring)   
  )
  
  return(npsem)
}

#' @export
#' @rdname survival_trans
survival_task <- function(data, node_list, make_npsem, variable_types = NULL, ...) {
  setDT(data)
  
  npsem <- make_npsem(node_list, variable_types)
  
  if (!is.null(node_list$id)) {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, id = node_list$id, time = node_list$time, ...)
  } else {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, ...)
  }
  
  return(tmle_task)
}

#' @export
#' @rdname survival_trans
survival_o_likelihood  <- function(tmle_task, learner_list) {
  # covariates
  # TODO: whether remove duplicates for LF_emp
  W_factor <- define_lf(LF_emp, "W")
  
  # TODO: check if necessary
  # treatment (bound likelihood away from 0 (and 1 if binary))
  S_type <- tmle_task$npsem[["S"]]$variable_type
  if (S_type$type == "continous") {
    S_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (S_type$type %in% c("binomial", "categorical")) {
    S_bound <- 0.025
  } else {
    S_bound <- NULL
  }
  
  S_factor <- define_lf(LF_fit, "S", learner = learner_list[["S"]], bound = S_bound)
  
  # TODO: modify get_regression_task and LF_fit for time variance
  # TODO: whether need bound
  outcome_bound <- 0.025
  F_factor <- define_lf(LF_fit_site, "failed", learner = learner_list[["failed"]],
                        is_time_variant = TRUE, bound = outcome_bound,
                        type = "mean")
  C_factor <- define_lf(LF_fit_site, "censored", learner = learner_list[["censored"]],
                        is_time_variant = TRUE, bound = outcome_bound,
                        type = "mean")
  
  factor_list <- list(W_factor, S_factor, F_factor, C_factor)
  
  # TODO: check
  likelihood <- Likelihood$new(factor_list)$train(tmle_task)
  return(likelihood)
}

#' @export
#' @rdname survival_trans
survival_e_likelihood  <- function(tmle_task, learner_list) {
  # covariates
  # TODO: whether remove duplicates for LF_emp
  W_factor <- define_lf(LF_emp, "W")
  
  # TODO: check if necessary
  # treatment (bound likelihood away from 0 (and 1 if binary))
  S_type <- tmle_task$npsem[["S"]]$variable_type
  if (S_type$type == "continous") {
    S_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (S_type$type %in% c("binomial", "categorical")) {
    S_bound <- 0.025
  } else {
    S_bound <- NULL
  }
  A_type <- tmle_task$npsem[["A"]]$variable_type
  if (A_type$type == "continous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }
  
  S_factor <- define_lf(LF_fit, "S", learner = learner_list[["S"]], bound = S_bound)
  A_factor <- define_lf(LF_fit_site, "A", learner = learner_list[["A"]], bound = A_bound)
  
  # TODO: modify get_regression_task and LF_fit for time variance
  # TODO: whether need bound
  outcome_bound <- 0.025
  F_factor <- define_lf(LF_fit_site, "failed", learner = learner_list[["failed"]],
                        is_time_variant = TRUE, bound = outcome_bound,
                        type = "mean")
  C_factor <- define_lf(LF_fit_site, "censored", learner = learner_list[["censored"]],
                        is_time_variant = TRUE, bound = outcome_bound,
                        type = "mean")
  
  factor_list <- list(W_factor, S_factor, A_factor, F_factor, C_factor)
  
  # TODO: check
  likelihood <- Likelihood$new(factor_list)$train(tmle_task)
  return(likelihood)
}

