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
Param_SOT <- R6Class(
  classname = "Param_SOT",
  portable = TRUE,
  class = TRUE,
  inherit = Param_AOT,
  public = list(
    initialize = function(observed_likelihood, target_times = NULL, 
                          onsite = 1, offsite = 0, 
                          fit_s_marginal = "empirical", ..., 
                          outcome_node = "F") {
      # TODO: check outcome_node, current I(T<=t, delta=1), need I(T=t, delta=1)
      super$initialize(observed_likelihood, onsite, offsite, fit_s_marginal, ..., outcome_node = outcome_node)
      private$.target_times <- target_times
    },
    clever_covariates_internal = function(tmle_task = NULL, fold_number = "full", subset_times = FALSE) {
      training_task <- self$observed_likelihood$training_task
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      IS1 <- self$cf_likelihood_onsite$get_likelihoods(tmle_task, "S", fold_number)
      
      if (self$fit_s_marginal == "empirical") {
        pS0 <- 1 - mean(training_task$data$S)
      } else if (self$fit_s_marginal == "integral") {
        cf_train_offsite <- training_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$offsite, training_task$nrow)))
        pS0 <- weighted.mean(self$observed_likelihood$get_likelihood(cf_train_offsite, "S", fold_number),
                             self$observed_likelihood$get_likelihood(cf_train_offsite, "W", fold_number))
      }
      
      cf_task_onsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$onsite, tmle_task$nrow)))
      cf_task_offsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$offsite, tmle_task$nrow)))
      pS1W <- self$observed_likelihood$get_likelihood(cf_task_onsite, "S", fold_number)
      pS0W <- self$observed_likelihood$get_likelihood(cf_task_offsite, "S", fold_number)
      
      hF <- self$observed_likelihood$get_likelihoods(tmle_task, "F", fold_number)
      # TODO: make bound configurable
      hF <- bound(hF, 0.005)
      
      hC <- self$observed_likelihood$get_likelihoods(tmle_task, "C", fold_number)
      
      time <- tmle_task$time
      id <- tmle_task$id
      
      t_mat <- long_to_mat(time,id,time)
      hF_mat <- long_to_mat(hF,id,time)
      hC_mat <- long_to_mat(hC,id,time)
      
      sF_mat <- hm_to_sm(hF_mat)
      sC_mat <- hm_to_sm(hC_mat)
      
      ks <- sort(unique(time))
      hk_all <- lapply(ks,function(k){
        Ikt <- k <= t_mat
        sF_mat_k <- matrix(sF_mat[,k+1],nrow=nrow(t_mat),ncol=ncol(t_mat))
        sC_mat_k <- matrix(sC_mat[,k],nrow=nrow(t_mat),ncol=ncol(t_mat))
        hk <- -1 * (Ikt/sC_mat_k)*(sF_mat[,-1]/sF_mat_k)
      })
      
      # TODO: this might need to be reordered
      H1 <- (IS1/pS1W)*(pS0W/pS0) * do.call(rbind, hk_all)
      
      if(subset_times & !is.null(self$target_times)){
        H1[,!(ks%in%self$target_times)] = 0
      }
      return(list(F = H1))
    },
    clever_covariates = function(tmle_task, fold_number = "full"){
      self$clever_covariates_internal(tmle_task, fold_number, subset_times = TRUE)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      training_task <- self$observed_likelihood$training_task
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      
      # TODO: return format
      # TODO: share work between this and the IC code
      H1 <- self$clever_covariates_internal(tmle_task, fold_number, subset_times = FALSE)[[self$outcome_node]]
      IS0 <- self$cf_likelihood_offsite$get_likelihoods(tmle_task, "S", fold_number)
      
      cf_task_onsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$onsite, tmle_task$nrow)))
      hFS1 <- self$observed_likelihood$get_likelihood(cf_task_onsite, self$outcome_node, fold_number)
      
      if (self$fit_s_marginal == "empirical") {
        pS0 <- 1 - mean(training_task$data$S)
      } else if (self$fit_s_marginal == "integral") {
        cf_train_offsite <- training_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$offsite, training_task$nrow)))
        pS0 <- weighted.mean(self$observed_likelihood$get_likelihood(cf_train_offsite, "S", fold_number),
                             self$observed_likelihood$get_likelihood(cf_train_offsite, "W", fold_number))
      }
      
      time <- tmle_task$time
      id <- tmle_task$id
      
      # TODO: make bound configurable
      hFS1 <- bound(hFS1, 0.005)
      hFS1_mat <- long_to_mat(hFS1,id,time)
      sFS1_mat <- hm_to_sm(hFS1_mat)
      psi <- colMeans(sFS1_mat[,-1])
      T_tilde <- tmle_task$get_tmle_node("T")
      Delta <- tmle_task$get_tmle_node("D")
      k <- time
      Fail <- (T_tilde == k) & (Delta==1)
      ITk <- (T_tilde >= k)
      
      D1_tk <- H1 * as.vector(Fail - (ITk * hFS1))
      # zero out entries that don't contribute to sum
      
      ts <- sort(unique(k))
      t_mat <- matrix(ts,nrow=nrow(D1_tk),ncol=ncol(D1_tk),byrow = TRUE)
      Itk <- (k<=t_mat)
      D1_tk <- D1_tk * Itk
      D1_tk_dt <- data.table(id=id, k=time, D1_tk)
      
      # sum to IC for 1:N ids
      D1 <- D1_tk_dt[,lapply(.SD,sum),by=list(id),.SDcols=as.character(ts)]
      D1 <- as.matrix(D1[,-1,with=FALSE])
      
      psi_mat <- matrix(psi,nrow=nrow(D1),ncol=ncol(D1),byrow=TRUE)
      D2 <- (IS0/pS0)*(sFS1_mat - psi_mat)

      IC_mat <- D1 + D2
      
      # copy IC to make it match the observation structure
      # TODO: consider if this is the best approach
      IC_id <- sort(unique(id))
      IC <- IC_mat[match(id,IC_id),]
      result <- list(psi = psi, IC = IC)
      return(result)
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
    target_times = function() {
      return(private$.target_times)
    }
  ),
  private = list(
    .type = "TMLE_SOT",
    .supports_outcome_censoring = TRUE,
    .target_times = NULL
  )
)
