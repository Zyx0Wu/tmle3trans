#' Average Treatment Effect
#'
#' Parameter definition for the Average Treatment Effect (ATE).
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
Param_TR <- R6Class(
  classname = "Param_TR",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, onsite = 1, offsite = 0,
                          ...,
                          covariate_node = "W", site_node = "A", outcome_node = "Y") {
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
      private$.onsite <- onsite
      private$.offsite <- offsite
      private$.covariate_node <- covariate_node
      private$.site_node <- site_node
      private$.cf_likelihood_onsite <- make_CF_Likelihood(observed_likelihood, define_lf(LF_static, self$site_node, value = self$onsite))
      private$.cf_likelihood_offsite <- make_CF_Likelihood(observed_likelihood, define_lf(LF_static, self$site_node, value = self$offsite))
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      training_task <- self$observed_likelihood$training_task
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      I1 <- self$cf_likelihood_onsite$get_likelihoods(tmle_task, self$site_node, fold_number)

      cf_train_offsite <- training_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(A = rep(self$offsite, training_task$nrow)))
      p0 <- mean(self$observed_likelihood$get_likelihood(cf_train_offsite, self$site_node, fold_number))

      cf_task_onsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(A = rep(self$onsite, tmle_task$nrow)))
      cf_task_offsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(A = rep(self$offsite, tmle_task$nrow)))
      p1W <- self$observed_likelihood$get_likelihood(cf_task_onsite, self$site_node, fold_number)
      p0W <- self$observed_likelihood$get_likelihood(cf_task_offsite, self$site_node, fold_number)

      H1 <- I1 / p1W * p0W / p0

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
      EYA <- self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number)
      I0 <- self$cf_likelihood_offsite$get_likelihoods(tmle_task, self$site_node, fold_number)

      cf_task_onsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(A = rep(self$onsite, tmle_task$nrow)))
      EY1 <- self$observed_likelihood$get_likelihood(cf_task_onsite, self$outcome_node, fold_number)

      cf_train_offsite <- training_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(A = rep(self$offsite, training_task$nrow)))
      p0 <- mean(self$observed_likelihood$get_likelihood(cf_train_offsite, self$site_node, fold_number))

      psi <- mean(I0/p0 * EY1)

      IC <- H1 * (Y - EYA) + I0/p0 * (EY1 - psi)

      result <- list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("E[E(%s | %s, trial) | reality]", self$outcome_node, self$covariate_node)
      return(param_form)
    },
    observed_likelihood_onsite = function() {
      return(private$.observed_likelihood_onsite)
    },
    covariate_node = function() {
      return(private$.covariate_node)
    },
    site_node = function() {
      return(private$.site_node)
    },
    onsite = function() {
      return(private$.onsite)
    },
    offsite = function() {
      return(private$.offsite)
    },
    cf_likelihood_onsite = function() {
      return(private$.cf_likelihood_onsite)
    },
    cf_likelihood_offsite = function() {
      return(private$.cf_likelihood_offsite)
    }
  ),
  private = list(
    .type = "TMLE_TR",
    .observed_likelihood_onsite = NULL,
    .covariate_node = NULL,
    .site_node = NULL,
    .onsite = NULL,
    .offsite = NULL,
    .cf_likelihood_onsite = NULL,
    .cf_likelihood_offsite = NULL
  )
)
