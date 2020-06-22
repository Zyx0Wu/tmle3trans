#' Likelihood Factor Estimated from Data using sl3.
#'
#' Uses an \code{sl3} learner to estimate a likelihood factor from data.
#' Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @importFrom tmle3 LF_fit
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_lf(LF_fit, name, learner, ..., type = "density")}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{learner}}{An sl3 learner to be used to estimate the factor
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{site}}{Value that indicates an observation being a trial data
#'     }
#'     \item{\code{site_node}}{Node of the trial indicator
#'     }
#'     \item{\code{type}}{character, either "density", for conditional density or, "mean" for conditional mean
#'     }
#'     }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{learner}}{The learner or learner fit object}
#'     }
#'
#' @export
LF_fit_site <- R6Class(
  classname = "LF_fit_site",
  portable = TRUE,
  class = TRUE,
  inherit = LF_fit,
  public = list(
    initialize = function(name, learner, ..., site = 1, site_node = "A", type = "density") {
      super$initialize(name, learner, ..., type = type)
      private$.site <- site
      private$.site_node <- site_node
    },
    delayed_train = function(tmle_task) {
      # just return prefit learner if that's what we have
      # otherwise, make a delayed fit and return that
      if (self$learner$is_trained) {
        return(self$learner)
      }
      
      outcome_node <- self$name
      
      # fit scaled task for bounded continuous
      site_data <- tmle_task$get_tmle_node(self$site_node)
      tmle_task_site <- tmle_task$subset_task(site_data == self$site)
      learner_task_site <- tmle_task_site$get_regression_task(outcome_node, scale = TRUE)
      learner_fit <- delayed_learner_train(self$learner, learner_task_site)
      return(learner_fit)
    }
  ),
  active = list(
    site = function() {
      return(private$.site)
    },
    site_node = function() {
      return(private$.site_node)
    }
  ),
  private = list(
    .site = NULL,
    .site_node = NULL
  )
)
