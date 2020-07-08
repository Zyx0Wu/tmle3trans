#' Survival Outcome Transportation
#'
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
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, target_times = NULL, 
                          onsite = 1, offsite = 0, 
                          ..., 
                          outcome_node = "N") {
      # TODO: check outcome_node, current I(T<=t, delta=1), need I(T=t, delta=1)
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
      private$.onsite <- onsite
      private$.offsite <- offsite
      private$.cf_likelihood_onsite <- make_CF_Likelihood(observed_likelihood, define_lf(LF_static, "S", value = self$onsite))
      private$.cf_likelihood_offsite <- make_CF_Likelihood(observed_likelihood, define_lf(LF_static, "S", value = self$offsite))
      private$.target_times <- target_times
    },
    long_to_mat = function(x,id, time){
      dt <- data.table(id=id,time=time,x=as.vector(x))
      wide <- dcast(dt, id~time, value.var="x")
      mat <- as.matrix(wide[,-1,with=FALSE])
      return(mat)
    },
    hm_to_sm = function(hm){
      # TODO: check
      sm <- t(apply(1-hm,1,cumprod))
      sm <- cbind(1,sm[,-ncol(sm)])
      return(sm)
    },
    clever_covariates_internal = function(tmle_task = NULL, fold_number = "full", subset_times = FALSE) {
      training_task <- self$observed_likelihood$training_task
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      I1 <- self$cf_likelihood_onsite$get_likelihoods(tmle_task, "S", fold_number)
      
      cf_train_offsite <- training_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$offsite, training_task$nrow)))
      p0 <- mean(self$observed_likelihood$get_likelihood(cf_train_offsite, "S", fold_number))
      
      cf_task_onsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$onsite, tmle_task$nrow)))
      cf_task_offsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$offsite, tmle_task$nrow)))
      p1W <- self$observed_likelihood$get_likelihood(cf_task_onsite, "S", fold_number)
      p0W <- self$observed_likelihood$get_likelihood(cf_task_offsite, "S", fold_number)
      
      pN <- self$observed_likelihood$get_likelihoods(tmle_task, "N", fold_number)
      # TODO: make bound configurable
      pN <- bound(pN, 0.005)
      
      pA_c <- self$observed_likelihood$get_likelihoods(tmle_task, "A_c", fold_number)
      
      time <- tmle_task$time
      id <- tmle_task$id
      long_order <- order(id,time)
      
      I1_mat <- self$long_to_mat(I1,id,time)
      p0_mat <- self$long_to_mat(p0,id,time)
      p1W_mat <- self$long_to_mat(p1W,id,time)
      p0W_mat <- self$long_to_mat(p0W,id,time)
      t_mat <- self$long_to_mat(time,id,time)
      
      pN_mat <- self$long_to_mat(pN,id,time)
      pA_c_mat <- self$long_to_mat(pA_c,id,time)
      SN_mat <- self$hm_to_sm(pN_mat)
      SA_c_mat <- self$hm_to_sm(pA_c_mat)
      
      ks <- sort(unique(time))
      
      hk_all <- lapply(ks,function(k){
        Ikt <- k <= t_mat
        SN_mat_k <- matrix(SN_mat[,k],nrow=nrow(t_mat),ncol=ncol(t_mat))
        SA_c_mat_k <- matrix(SA_c_mat[,k],nrow=nrow(t_mat),ncol=ncol(t_mat))
        hk <- -1 * ((I1*p0W*Ikt)/(p1W*p0*SA_c_mat_k))*(SN_mat/SN_mat_k)
      })
      
      # TODO: this might need to be reordered
      H1 <- do.call(rbind, hk_all)
      
      if(subset_times & !is.null(self$target_times)){
        H1[,!(ks%in%self$target_times)] = 0
      }
      return(list(N = H1))
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
      I0 <- self$cf_likelihood_offsite$get_likelihoods(tmle_task, "S", fold_number)
      
      cf_task_onsite <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$onsite, tmle_task$nrow)))
      pN1 <- self$observed_likelihood$get_likelihood(cf_task_onsite, self$outcome_node, fold_number)
      
      cf_train_offsite <- training_task$generate_counterfactual_task(UUIDgenerate(), new_data = data.table(S = rep(self$offsite, training_task$nrow)))
      p0 <- mean(self$observed_likelihood$get_likelihood(cf_train_offsite, "S", fold_number))
      
      time <- tmle_task$time
      id <- tmle_task$id
      
      # TODO: make bound configurable
      pN1 <- bound(pN1, 0.005)
      pN1_mat <- self$long_to_mat(pN1,id,time)
      SN1_mat <- self$hm_to_sm(pN1_mat)
      psi <- colMeans(SN1_mat)
      T_tilde <- tmle_task$get_tmle_node("T_tilde")
      Delta <- tmle_task$get_tmle_node("Delta")
      k <- time
      Ittkd <- (T_tilde == k) & (Delta==1)
      Ittk <- (T_tilde >= k)
      
      resid <- as.vector(Ittkd - (Ittk * pN1))
      D1_tk <- H1*resid
      
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
      I0_mat <- self$long_to_mat(I0,id,time)
      p0_mat <- self$long_to_mat(p0,id,time)
      D2 <- (I0/p0)*(SN1_mat - psi_mat)
      
      IC <- D1+D2
      
      # copy IC to make it match the observation structure
      # TODO: consider if this is the best approach
      IC_id <- sort(unique(id))
      IC_long <- IC[match(id,IC_id),]
      result <- list(psi = psi, IC = IC_long)
      return(result)
    }
  ),
  active = list(
    # TODO: modify
    name = function() {
      param_form <- sprintf("E[E(%s | W, trial) | reality]", self$outcome_node)
      return(param_form)
    },
    cf_likelihood = function() {
      return(private$.cf_likelihood)
    },
    intervention_list = function() {
      return(self$cf_likelihood$intervention_list)
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
    .cf_likelihood = NULL,
    .supports_outcome_censoring = TRUE,
    .target_times = NULL
  )
)
