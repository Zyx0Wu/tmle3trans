#' Average Outcome Transportation
#'
#' Parameter definition for the Average Outcome Transportation (AOT).
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @importFrom tmle3 Param_base define_lf make_CF_Likelihood
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_AOT, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention.
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood_onsite}}{the counterfactual likelihood for clinical trial
#'     }
#'     \item{\code{cf_likelihood_offsite}}{the counterfactual likelihood for real world
#'     }
#' }
#' @export
Param_AOT <- R6Class(
  classname = "Param_AOT",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, onsite = 1, offsite = 0, 
                          fit_s_marginal = "empirical", ..., outcome_node = "Y") {
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
      private$.onsite <- onsite
      private$.offsite <- offsite
      private$.fit_s_marginal <- fit_s_marginal
      private$.cf_likelihood_onsite <- 
        make_CF_Likelihood(observed_likelihood, 
                           define_lf(LF_static, "S", value = self$onsite))
      private$.cf_likelihood_offsite <- 
        make_CF_Likelihood(observed_likelihood, 
                           define_lf(LF_static, "S", value = self$offsite))
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      training_task <- self$observed_likelihood$training_task
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      IS1 <- self$cf_likelihood_onsite$get_likelihoods(tmle_task, "S", fold_number)

      if (self$fit_s_marginal == "empirical") {
        pS0 <- 1 - mean(training_task$get_tmle_node("S"))
      } else if (self$fit_s_marginal == "integral") {
        cf_train_offsite <- training_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$offsite, training_task$nrow)))
        pS0 <- weighted.mean(self$observed_likelihood$get_likelihood(cf_train_offsite, "S", fold_number),
                             self$observed_likelihood$get_likelihood(cf_train_offsite, "W", fold_number))
      }
      
      cf_task_onsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$onsite, tmle_task$nrow)))
      cf_task_offsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$offsite, tmle_task$nrow)))
      pS1W <- self$observed_likelihood$get_likelihood(cf_task_onsite, "S", fold_number)
      pS0W <- self$observed_likelihood$get_likelihood(cf_task_offsite, "S", fold_number)
      
      H1 <- IS1 / pS1W * pS0W / pS0
      #H1 <- IS1 / prob_clip(pS1W) * pS0W / prob_clip(pS0)
      return(list(Y = H1))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      training_task <- self$observed_likelihood$training_task
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      # clever_covariates happen here (for this param) only, but this is repeated computation
      H1 <- self$clever_covariates(tmle_task, fold_number)[[self$outcome_node]]
      Y <- tmle_task$get_tmle_node(self$outcome_node)
      EYS <- self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number)
      IS0 <- self$cf_likelihood_offsite$get_likelihoods(tmle_task, "S", fold_number)

      cf_task_onsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$onsite, tmle_task$nrow)))
      EYS1 <- self$observed_likelihood$get_likelihood(cf_task_onsite, self$outcome_node, fold_number)
      
      if (self$fit_s_marginal == "empirical") {
        pS0 <- 1 - mean(training_task$get_tmle_node("S"))
      } else if (self$fit_s_marginal == "integral") {
        cf_train_offsite <- training_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$offsite, training_task$nrow)))
        pS0 <- weighted.mean(self$observed_likelihood$get_likelihood(cf_train_offsite, "S", fold_number),
                             self$observed_likelihood$get_likelihood(cf_train_offsite, "W", fold_number))
      }
      
      psi <- mean(IS0/pS0 * EYS1)
      IC <- H1 * (Y - EYS) + IS0/pS0 * (EYS1 - psi)
      #psi <- mean(IS0/prob_clip(pS0) * EYS1)
      #IC <- H1 * (Y - EYS) + IS0/prob_clip(pS0) * (EYS1 - psi)
      return(list(psi = psi, IC = IC))
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("E[E(%s | W, trial) | reality]", self$outcome_node)
      return(param_form)
    },
    update_nodes = function() {
      return(c(self$outcome_node))
    },
    onsite = function() {
      return(private$.onsite)
    },
    offsite = function() {
      return(private$.offsite)
    },
    fit_s_marginal = function() {
      return(private$.fit_s_marginal)
    },
    cf_likelihood_onsite = function() {
      return(private$.cf_likelihood_onsite)
    },
    cf_likelihood_offsite = function() {
      return(private$.cf_likelihood_offsite)
    }
  ),
  private = list(
    .type = "TMLE_AOT",
    .onsite = NULL,
    .offsite = NULL,
    .fit_s_marginal = NULL,
    .cf_likelihood_onsite = NULL,
    .cf_likelihood_offsite = NULL
  )
)
