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
    initialize = function(onsite = 1, offsite = 0, ...) {
      super$initialize(onsite = onsite, offsite = offsite, ...)
    },
    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types

      tmle_task <- point_tr_task(data, node_list, variable_types)

      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      likelihood <- point_tr_likelihood(tmle_task, learner_list)
      return(likelihood)
    },
    make_params = function(tmle_task, targeted_likelihood) {
      tmle_params <- define_param(Param_AET, targeted_likelihood,
                                  onsite = self$options$onsite,
                                  offsite = self$options$offsite)
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

#' All Treatment Specific Means
#'
#' O=(W,S,Y)
#' W=Covariates
#' S=Location
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @param onsite value for onsite
#' @param offsite value for offsite
#' @export
tmle_AET <- function(onsite = 1, offsite = 0, ...) {
  tmle3_Spec_AET$new(onsite = onsite, offsite = offsite)
}


