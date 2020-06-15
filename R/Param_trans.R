#' Average Treatment Effect
#'
#' Parameter definition for the Average Treatment Effect (ATE).
#' @importFrom R6 R6Class
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
#'   \code{define_param(Param_ATT, observed_likelihood, intervention_list, ..., outcome_node)}
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
#'     \item{\code{cf_likelihood_treatment}}{the counterfactual likelihood for the treatment
#'     }
#'     \item{\code{cf_likelihood_control}}{the counterfactual likelihood for the control
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention
#'     }
#' }
#' @export
Param_ATE <- R6Class(
  classname = "Param_ATE",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, onsite = "1", offsite = "0",
                          ...,
                          covariate_node = "W", site_node = "A", outcome_node = "Y") {
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
      private$.covariate_node <- covariate_node
      private$.site_node <- site_node
      private$.onsite <- onsite
      private$.offsite <- offsite
      private$.cf_likelihood_onsite <- make_CF_Likelihood(observed_likelihood, define_lf(LF_static, self$site_node, value = self$onsite))
      private$.cf_likelihood_offsite <- make_CF_Likelihood(observed_likelihood, define_lf(LF_static, self$site_node, value = self$offsite))
      
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      
      cf_pS_onsite <- self$cf_likelihood_onsite$get_likelihoods(tmle_task, self$site_node, fold_number)
      #cf_pS_offsite <- self$cf_likelihood_offsite$get_likelihoods(tmle_task, self$site_node, fold_number)
      
      cf_task_onsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(A = rep(self$onsite, tmle_task$nrow)))
      cf_task_offsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(A = rep(self$offsite, tmle_task$nrow)))
      p1W <- self$observed_likelihood$get_likelihood(cf_task_onsite, self$site_node, fold_number)
      p0W <- self$observed_likelihood$get_likelihood(cf_task_offsite, self$site_node, fold_number)
      
      p0 <- mean(p0W)
      HA <- cf_pS_onsite / p1W * p0W / p0
      
      return(list(Y = HA))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))
      
      # clever_covariates happen here (for this param) only, but this is repeated computation
      HA <- self$clever_covariates(tmle_task, fold_number)[[self$outcome_node]]
      
      
      # todo: make sure we support updating these params
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      
      # todo: extend for stochastic
      cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
      
      Y <- tmle_task$get_tmle_node(self$outcome_node, impute_censoring = TRUE)
      
      EY <- self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number)
      EY1 <- self$observed_likelihood$get_likelihood(cf_task_treatment, self$outcome_node, fold_number)
      EY0 <- self$observed_likelihood$get_likelihood(cf_task_control, self$outcome_node, fold_number)
      
      psi <- mean(EY1 - EY0)
      
      IC <- HA * (Y - EY) + (EY1 - EY0) - psi
      
      result <- list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("ATE[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
      return(param_form)
    },
    covariate_node = function() {
      return(private$.covariate_node)
    },
    site_node = function() {
      return(private$.site_node)
    },
  ),
  private = list(
    .type = "ATE",
    .covariate_node = NULL,
    .site_node = NULL
  )
)
