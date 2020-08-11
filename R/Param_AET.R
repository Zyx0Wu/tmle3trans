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
                          fit_s_marginal = "empirical", ...,
                          intervention_node = "A", outcome_node = "Y") {
      super$initialize(observed_likelihood, onsite, offsite, fit_s_marginal, 
                       ..., outcome_node = outcome_node)
      private$.intervention <- intervention
      private$.intervention_node <- intervention_node
      private$.cf_likelihood_intervention <- 
        make_CF_Likelihood(observed_likelihood, 
                           define_lf(LF_static, self$intervention_node, value = self$intervention))
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      S <- tmle_task$get_tmle_node(self$site_node)
      A_trans <- tmle_task$get_tmle_node(self$intervention_node)
      A_trans[S == self$offsite] <- self$intervention
      
      tmle_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(A = A_trans))
      
      base_covs <- super$clever_covariates(tmle_task, fold_number)
      
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, self$intervention_node, fold_number)
      IA <- self$cf_likelihood_intervention$get_likelihoods(tmle_task, self$intervention_node, fold_number)
      
      clever_covs <- lapply(base_covs, `*`, IA/pA)
      #clever_covs <- lapply(base_covs, `*`, IA/prob_clip(pA))
      return(clever_covs)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      S <- tmle_task$get_tmle_node(self$site_node)
      A_trans <- tmle_task$get_tmle_node(self$intervention_node)
      A_trans[S == self$offsite] <- self$intervention
      
      tmle_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(A = A_trans))
      return(super$estimates(tmle_task, fold_number))
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
    intervention_node = function() {
      return(private$.intervention_node)
    },
    cf_likelihood_intervention = function() {
      return(private$.cf_likelihood_intervention)
    }
  ),
  private = list(
    .type = "TMLE_AET",
    .intervention = NULL,
    .intervention_node = NULL,
    .cf_likelihood_intervention = NULL
  )
)
