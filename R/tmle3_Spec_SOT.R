#' Defines a TMLE with Survival Outcome Transportation
#'
#'
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'
#' @export
#

tmle3_Spec_SOT <- R6Class(
  classname = "tmle3_Spec_SOT",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(onsite = 1, offsite = 0, target_times = NULL, 
                          fit_s_marginal = "empirical", ...) {
      super$initialize(onsite = onsite, offsite = offsite, 
                       target_times = target_times,
                       fit_s_marginal = fit_s_marginal, ...)
    },
    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types

      tmle_task <- survival_task(data, node_list, survival_o_npsem, variable_types)

      return(tmle_task)
    },

    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      likelihood <- survival_o_likelihood(tmle_task, learner_list)
      return(likelihood)
    },

    make_params = function(tmle_task, likelihood) {
      tmle_params <- define_param(Param_SOT, likelihood,
                                  onsite = self$options$onsite,
                                  offsite = self$options$offsite,
                                  target_times = self$options$target_times,
                                  fit_s_marginal = self$options$fit_s_marginal)
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

#' Wrapper
#'
#' @importFrom sl3 make_learner Lrnr_mean
#' @param onsite value for onsite
#' @param offsite value for offsite
#' @export
tmle_SOT <- function(onsite = 1, offsite = 0, target_times = NULL, 
                     fit_s_marginal = "empirical", ...) {
  tmle3_Spec_SOT$new(onsite = onsite, offsite = offsite,
                     target_times = target_times,
                     fit_s_marginal = fit_s_marginal, ...)
}


