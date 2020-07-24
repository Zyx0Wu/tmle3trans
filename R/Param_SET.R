#' Survival Outcome Transportation
#'
#' @importFrom R6 R6Class
#' @import data.table
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_survival, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood}}{the counterfactual likelihood for this treatment
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention
#'     }
#' }
#' @export
Param_SET <- R6Class(
  classname = "Param_SET",
  portable = TRUE,
  class = TRUE,
  inherit = Param_SOT,
  public = list(
    initialize = function(observed_likelihood, intervention, 
                          onsite = 1, offsite = 0, target_times = NULL, 
                          fit_s_marginal = "empirical", ..., 
                          outcome_node = "F") {
      # TODO: check outcome_node, current I(T<=t, delta=1), need I(T=t, delta=1)
      super$initialize(observed_likelihood, target_times, onsite, offsite, fit_s_marginal, ..., outcome_node = outcome_node)
      private$.intervention <- intervention
      private$.cf_likelihood_intervention <- 
        make_CF_Likelihood(observed_likelihood, 
                           define_lf(LF_static, "A", value = self$intervention))
    },
    clever_covariates_internal = function(tmle_task = NULL, fold_number = "full", subset_times = FALSE) {
      base_covs <- super$clever_covariates_internal(tmle_task, fold_number, subset_times)
      
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, "A", fold_number)
      IA <- self$cf_likelihood_intervention$get_likelihoods(tmle_task, "A", fold_number)
      
      clever_covs <- lapply(base_covs, `*`, IA/pA)
      #clever_covs <- lapply(base_covs, `*`, IA/prob_clip(pA))
      return(clever_covs)
    }
  ),
  active = list(
    # TODO: modify
    name = function() {
      param_form <- sprintf("E[P(T > t | W, trial) | reality]")
      return(param_form)
    },
    update_nodes = function() {
      return(self$outcome_node)
    },
    intervention = function() {
      return(private$.intervention)
    },
    cf_likelihood_intervention = function() {
      return(private$.cf_likelihood_intervention)
    }
  ),
  private = list(
    .type = "TMLE_SET",
    .intervention = NULL,
    .cf_likelihood_intervention = NULL
  )
)
