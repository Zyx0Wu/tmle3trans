#' Defines a TMLE with Average Effect Transportation
#'
#'
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'
#' @export
#
tmle3_Spec_AET <- R6Class(
  classname = "tmle3_Spec_AET",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(onsite = 1, offsite = 0, fit_s_marginal = "empirical", ...) {
      super$initialize(onsite = onsite, offsite = offsite, fit_s_marginal = fit_s_marginal, ...)
    },
    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types
      tmle_task <- point_task(data, node_list, point_e_npsem, variable_types)

      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      likelihood <- point_e_likelihood(tmle_task, learner_list)

      return(likelihood)
    },
    make_params = function(tmle_task, likelihood) {
      # todo: export and use sl3:::get_levels
      A_vals <- tmle_task$get_tmle_node("A")
      if (is.factor(A_vals)) {
        A_levels <- levels(A_vals)
        A_levels <- factor(A_levels, A_levels)
      } else {
        A_levels <- sort(unique(A_vals))
      }
      tmle_params <- lapply(A_levels, function(A_level) {
        tmle_param <- define_param(Param_AET, likelihood, A_level,
                                   onsite = self$options$onsite,
                                   offsite = self$options$offsite,
                                   fit_s_marginal = self$options$fit_s_marginal)
        return(tmle_param)
      })

      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

#' All Treatment Specific Means
#'
#' O=(W,S,A,Y)
#' W=Covariates
#' S=Location
#' A=Intervention
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @param onsite value for onsite
#' @param offsite value for offsite
#' @export
tmle_AET <- function(onsite = 1, offsite = 0, fit_s_marginal = "empirical", ...) {
  tmle3_Spec_AET$new(onsite = onsite, offsite = offsite, fit_s_marginal = fit_s_marginal, ...)
}


