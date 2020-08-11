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
    node("A", distr = "rbinom", size = 1, prob = .15 + .5 * as.numeric(W > .75)) +
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
  # only grab S, W's, A, T.tilde, Delta
  Wname <- grep("W", colnames(dat), value = TRUE)
  dat <- dat[, c("S", Wname, "A", "T.tilde", "Delta")]
  return(list(dat = dat))
}

set.seed(4321)
n_sim <- 2e2
simulated <- simulate_data(n_sim = n_sim)
data <- simulated$dat

# TODO: check
while (all(data$Delta == 1)) {
  simulated <- simulate_data(n_sim = n_sim)
  data <- simulated$dat
}

tmax <- max(data$T.tilde)
times <- 1:tmax

node_list <- list(S = "S", W = c("W", "W1"), A = "A", 
                  T.tilde = "T.tilde", Delta = "Delta")

transformed <- transform_data(data, node_list)
long_data <- transformed$long_data
long_node_list <- transformed$long_node_list

### sl3 ###
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_sl <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm))
learner_list <- list(S = lrnr_sl, A = lrnr_sl, failed = lrnr_sl, censored = lrnr_sl)

### tmle3 ###
# TODO: check
# survival curve for treatment 1 only targeting time point 1
tmle_spec <- tmle_SET(1, target_times = intersect(1, times))
tmle_task <- tmle_spec$make_tmle_task(long_data, long_node_list)

initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

updater <- tmle3_Update$new(maxit=5, verbose = TRUE)
'
updater <- tmle3_Update$new(
  constrain_step = TRUE, one_dimensional = TRUE, 
  delta_epsilon = 3e-2, verbose = TRUE, 
  convergence_type = "scaled_var", maxit = 10
)

updater <- tmle3_Update_survival$new(
  maxit = 3e1, 
  cvtmle = TRUE,
  convergence_type = "scaled_var",
  delta_epsilon = 1e-2,
  fit_method = "l2",
  use_best = TRUE,
  verbose = TRUE
)
'
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = updater)
tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
# initial mean IC
mean(tmle_params$estimates(tmle_task, "validation")$IC[,1])

tmle_fit_manual <- fit_tmle3(
  tmle_task, targeted_likelihood, tmle_params,
  targeted_likelihood$updater
)
# tmle mean IC
mean(tmle_params$estimates(tmle_task, "validation")$IC[,1])

