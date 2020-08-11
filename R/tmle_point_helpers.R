#' Helper Functions for Point Estimators
#'
#' Handles the common W (covariates), S (location), Y (outcome) data structure
#'
#' @param data a \code{data.frame}, or \code{data.table} containing data for use in estimation
#' @param node_list a list of character vectors, listing the variables that comprise each node
#' @param variable_types a list of variable types, one for each node. If missing, variable types will be guessed
#' @param tmle_task a \code{\link{tmle3_Task}} as constructed via \code{point_task}
#' @param learner_list a list of sl3 learners, one for A and one for Y to be used for likelihood estimation
#' @param ... extra arguments.
#' 
#' @export
#' @rdname point_trans
point_o_npsem <- function(node_list, variable_types = NULL) {
  npsem <- list(
    define_node("W", node_list$W, variable_type = variable_types$W),
    define_node("S", node_list$S, c("W"), variable_type = variable_types$S),
    define_node("Y", node_list$Y, c("W"), variable_type = variable_types$Y, scale = TRUE)
  )
  
  return(npsem)
}

#' @export
#' @rdname point_trans
point_e_npsem <- function(node_list, variable_types = NULL) {
  npsem <- list(
    define_node("W", node_list$W, variable_type = variable_types$W),
    define_node("S", node_list$S, c("W"), variable_type = variable_types$S),
    define_node("A", node_list$A, c("W"), variable_type = variable_types$A),
    define_node("Y", node_list$Y, c("A", "W"), variable_type = variable_types$Y, scale = TRUE)
  )
  
  return(npsem)
}

#' @export
#' @rdname point_trans
point_task <- function(data, node_list, make_npsem, variable_types = NULL, ...) {
  setDT(data)
  
  npsem <- make_npsem(node_list, variable_types)
  
  if (!is.null(node_list$id)) {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, id = node_list$id, ...)
  } else {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, ...)
  }
  
  return(tmle_task)
}

#' @export
#' @rdname point_trans
point_o_likelihood <- function(tmle_task, learner_list) {
  # covariates
  W_factor <- define_lf(LF_emp, "W")
  
  # transportation (bound likelihood away from 0 (and 1 if binary))
  S_type <- tmle_task$npsem[["S"]]$variable_type
  if (S_type$type == "continous") {
    S_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (S_type$type %in% c("binomial", "categorical")) {
    S_bound <- 0.025
  } else {
    S_bound <- NULL
  }
  
  S_factor <- define_lf(LF_fit, "S", learner = learner_list[["S"]], bound = S_bound)
  
  # outcome
  Y_factor <- define_lf(LF_fit_site, "Y", learner = learner_list[["Y"]], type = "mean")
  
  # construct and train likelihood
  factor_list <- list(W_factor, S_factor, Y_factor)
  
  likelihood_def <- Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}

#' @export
#' @rdname point_trans
point_e_likelihood <- function(tmle_task, learner_list) {
  # covariates
  W_factor <- define_lf(LF_emp, "W")
  
  # transportation (bound likelihood away from 0 (and 1 if binary))
  S_type <- tmle_task$npsem[["S"]]$variable_type
  if (S_type$type == "continous") {
    S_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (S_type$type %in% c("binomial", "categorical")) {
    S_bound <- 0.025
  } else {
    S_bound <- NULL
  }
  
  S_factor <- define_lf(LF_fit, "S", learner = learner_list[["S"]], bound = S_bound)
  
  # treatment (bound likelihood away from 0 (and 1 if binary))
  A_type <- tmle_task$npsem[["A"]]$variable_type
  if (A_type$type == "continous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }
  
  A_factor <- define_lf(LF_fit_site, "A", learner = learner_list[["A"]], bound = A_bound)
  
  # outcome
  Y_factor <- define_lf(LF_fit_site, "Y", learner = learner_list[["Y"]], type = "mean")
  
  # construct and train likelihood
  factor_list <- list(W_factor, S_factor, A_factor, Y_factor)
  
  likelihood <- Likelihood$new(factor_list)$train(tmle_task)
  return(likelihood)
}

