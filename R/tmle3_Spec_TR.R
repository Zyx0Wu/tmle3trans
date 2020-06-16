#' Defines a TMLE with site transportation
#'
#'
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'
#' @export
#
tmle3_Spec_TR <- R6Class(
  classname = "tmle3_Spec_TR",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(onsite = 1, offsite = 0, ...) {
      super$initialize(onsite = 1, offsite = 0, ...)
    },
    set_initial_likelihood = function(likelihood) {
      private$.initial_likelihood <- likelihood
    },
    make_likelihood = function(tmle_task, learner_list = NULL, initial = FALSE) {
      # produce trained likelihood when likelihood_def provided
      
      if (!is.null(self$options$likelihood_override)) {
        likelihood <- self$options$likelihood_override$train(tmle_task)
      } else {
        likelihood <- point_tx_likelihood(tmle_task, learner_list)
      }
      
      if (initial) {
        private$.initial_likelihood <- likelihood
      }
      return(likelihood)
    },
    make_params = function(tmle_task, targeted_likelihood_onsite) {
      tmle_params <- define_param(Param_TR, targeted_likelihood_onsite,
                                  initial_likelihood_all = self$initial_likelihood,
                                  onsite = self$options$onsite,
                                  offsite = self$options$offsite)
      return(tmle_params)
    }
  ),
  active = list(
    initial_likelihood = function() {
      return(private$.initial_likelihood)
    }
  ),
  private = list(
    .initial_likelihood = NULL
  )
)

#' All Treatment Specific Means
#'
#' O=(W,A,Y)
#' W=Covariates
#' A=Treatment (continuous)
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @param modification s(W)
#' @param mapping f(A, modification)
#' @param n_samples number of samples to draw for each observation
#' @export
tmle_TR <- function(initial_likelihood, onsite = 1, offsite = 0, ...) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_TR$new(initial_likelihood, onsite = 1, offsite = 0)
}