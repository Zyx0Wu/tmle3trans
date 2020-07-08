#' Average Effect Transportation
#'
#' Parameter definition for the Average Effect Transportation (AET).
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
#'   \code{define_param(Param_AET, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
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
Param_AET <- R6Class(
  classname = "Param_AET",
  portable = TRUE,
  class = TRUE,
  inherit = Param_AOT,
  public = list(
    initialize = function(observed_likelihood, intervention, 
                          onsite = 1, offsite = 0,
                          ...,
                          outcome_node = "Y") {
      super$initialize(observed_likelihood, onsite, offsite, ..., 
                       outcome_node = outcome_node)
      private$.intervention <- intervention
      private$.cf_likelihood_intervention <- 
        make_CF_Likelihood(observed_likelihood, 
                           define_lf(LF_static, "A", value = self$intervention))
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      base_covs <- super$clever_covariates(tmle_task, fold_number)
      
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, "A", fold_number)
      IA <- self$cf_likelihood_intervention$get_likelihoods(tmle_task, "A", fold_number)
      
      clever_covs <- lapply(base_covs, `*`, IA/pA)
      return(clever_covs)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("E[E(%s | W, %s, trial) | reality]", 
                            self$outcome_node, self$cf_likelihood_intervention$name)
      return(param_form)
    },
    update_nodes = function() {
      return(c(self$outcome_node))
    },
    intervention = function() {
      return(private$.intervention)
    },
    cf_likelihood_intervention = function() {
      return(private$.cf_likelihood_intervention)
    }
  ),
  private = list(
    .type = "TMLE_AET",
    .intervention = NULL,
    .cf_likelihood_intervention = NULL,
  )
)
