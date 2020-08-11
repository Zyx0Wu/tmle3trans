library(simcausal)
library(testthat)
library(data.table)
library(ggplot2)

library(sl3)
library(tmle3)
library(tmle3trans)

simulate_data <- function(n_sim = 2e2) {
  D <- DAG.empty()
  D <- D +
    node("S", distr = "rbinom", size = 1, prob = .6) +
    node("W1", distr = "rbinom", size = 1, prob = .5) +
    node("W", distr = "runif", min = .1 - .1 * S, max = 1.6 - .1 * S) +
    node("Trexp", distr = "rexp", rate = .2 + .7 * W^2) +
    node("Cweib", distr = "rweibull", shape = 1 + .5 * W, scale = 75) +
    node("T", distr = "rconst", const = round(Trexp * 2) + 1) +
    node("C", distr = "rconst", const = round(Cweib * 2) + 1) +
    # Observed random variable (follow-up time):
    node("T.tilde", distr = "rconst", const = ifelse(T <= C, T, C)) +
    # Observed random variable (censoring indicator, 1 - failure event, 0 - censored):
    node("Delta", distr = "rconst", const = ifelse(T <= C, 1, 0))
  setD <- set.DAG(D)
  dat <- sim(setD, n = n_sim)
  # only grab W's, S, T.tilde, Delta
  Wname <- grep("W", colnames(dat), value = TRUE)
  dat <- dat[, c(Wname, "S", "T.tilde", "Delta")]
  return(list(dat = dat))
}

set.seed(4321)
n_sim <- 2e4
simulated <- simulate_data(n_sim = n_sim)
data <- simulated$dat

# TODO: check
while (all(data$Delta == 1)) {
  simulated <- simulate_data(n_sim = n_sim)
  data <- simulated$dat
}

tmax <- max(data$T.tilde)
times <- 1:tmax

node_list <- list(W = c("W", "W1"), S = "S", T.tilde = "T.tilde", Delta = "Delta")

transformed <- transform_data(data, node_list)
long_data <- transformed$long_data
long_node_list <- transformed$long_node_list

# observed survival curve
max_time <- 20
sF <- 1
for (i in times) {
  sF <- c(sF, sF[length(sF)] - mean(long_data[((i-1)*n_sim+1):(i*n_sim),][data$S==0,]$Failed))
}
sF <- sF[-1][1:max_time]

seed_fit = function(seed, n_sim) {
  set.seed(seed)
  simulated <- simulate_data(n_sim = n_sim)
  data <- simulated$dat
  
  # TODO: check
  while (all(data$Delta == 1) | max(data$T.tilde) < max_time) {
    simulated <- simulate_data(n_sim = n_sim)
    data <- simulated$dat
  }
  
  tmax <- max(data$T.tilde)
  times <- 1:tmax
  
  node_list <- list(W = c("W", "W1"), S = "S", T.tilde = "T.tilde", Delta = "Delta")
  
  transformed <- transform_data(data, node_list)
  long_data <- transformed$long_data
  long_node_list <- transformed$long_node_list
  
  ### sl3 ###
  lrnr_mean <- make_learner(Lrnr_mean)
  lrnr_glm <- make_learner(Lrnr_glm)
  lrnr_sl <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm))
  learner_list <- list(S = lrnr_sl, failed = lrnr_sl, censored = lrnr_sl)
  
  ### tmle3 ###
  # TODO: check
  # only target time point 1
  tmle_spec <- tmle_SOT(1, 0, target_times = times[1:min(5, tmax)])
  tmle_task <- tmle_spec$make_tmle_task(long_data, long_node_list)
  
  initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)
  '
  updater <- tmle3_Update$new(maxit=5, verbose = TRUE)
  
  updater <- tmle3_Update$new(
    constrain_step = TRUE, one_dimensional = TRUE, 
    delta_epsilon = 3e-2, verbose = TRUE, 
    convergence_type = "scaled_var", maxit = 10
  )
  '
  updater <- tmle3_Update_survival$new(
    maxit = 1e2, 
    cvtmle = TRUE,
    convergence_type = "scaled_var",
    delta_epsilon = 1e-2,
    fit_method = "l2",
    use_best = TRUE,
    verbose = TRUE
  )
  
  targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = updater)
  tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
  
  tmle_fit_manual <- fit_tmle3(
    tmle_task, targeted_likelihood, tmle_params,
    targeted_likelihood$updater
  )
  
  # extract results
  tmle_summary <- tmle_fit_manual$summary
  init_psi <- tmle_summary$init_est[1:max_time]
  tmle_psi <- tmle_summary$tmle_est[1:max_time]
  
  # loss
  l2_diff_sl <- (init_psi - sF)^2
  l2_diff_tl <- (tmle_psi - sF)^2
  
  return(list(init_psi, tmle_psi, l2_diff_sl, l2_diff_tl))
}

reps <- 50
obs <- 2e2
sl_psis <- c()
tl_psis <- c()
diff_sl <- c()
diff_tl <- c()
for (i in 1:reps) {
  fits <- seed_fit(1234+i, obs)
  sl_psis <- cbind(sl_psis, fits[[1]])
  tl_psis <- cbind(tl_psis, fits[[2]])
  diff_sl <- cbind(diff_sl, fits[[3]])
  diff_tl <- cbind(diff_tl, fits[[4]])
}

