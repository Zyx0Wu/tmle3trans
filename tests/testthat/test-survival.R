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
  # only grab ID, W's, S, T.tilde, Delta
  Wname <- grep("W", colnames(dat), value = TRUE)
  dat <- dat[, c("ID", Wname, "S", "T.tilde", "Delta")]
  # input: scalar q, W vector. computes for all W, the S(q|S,W)
  true_surv_w <- function(q, W) sapply(W, function(w) {
    1 - pexp(q, rate = .2 + .7 * w^2)
  })
  # input: vector q. mean(S(q|S=s,W)), average out W. loop over q
  true_surv_s <- function(q_grid, surv_fn, s) {
    W_grid <- seq(.1 - .1 * s, 1.6 - .1 * s, .01)
    survout <- numeric()
    for (q in q_grid) survout <- c(survout, mean(surv_fn(q = q / 2, W = W_grid)))
    return(survout)
  }
  truth_surv1 <- function(q) true_surv_s(q_grid = q, surv_fn = true_surv_w, s = 1)
  truth_surv0 <- function(q) true_surv_s(q_grid = q, surv_fn = true_surv_w, s = 0)
  return(list(dat = dat, true_surv1 = truth_surv1, true_surv0 = truth_surv0))
}

set.seed(1234)
n_sim <- 1e2
simulated <- simulate_data(n_sim = n_sim)
df <- simulated$dat

# TODO: check
while (all(df$Delta == 1)) {
  simulated <- simulate_data(n_sim = n_sim)
  df <- simulated$dat
}

###########
# tlverse #
###########
tmax <- max(df$T.tilde)
k_grid <- 1:tmax
all_times <- lapply(k_grid, function(t) df_time(df, t))

df_long <- rbindlist(all_times)

node_list <- list(id ="ID", W = c("W", "W1"), S = "S", T = "T.tilde", D = "Delta", 
                  time = "t", F = "Failed", C = "Censored", pre_failure = "pre_failure")

### sl3 ###
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_gam <- make_learner(Lrnr_gam)
lrnr_sl <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm, lrnr_gam))
learner_list <- list(S = lrnr_sl, F = lrnr_sl, C = lrnr_sl)

### tmle3 ###
# TODO: check
var_types <- list(T = Variable_Type$new("continuous"), t = Variable_Type$new("continuous"), 
                  D = Variable_Type$new("binomial"))
tmle_spec <- tmle_SOT(1, 0, target_times = intersect(1, k_grid), variable_types = var_types)
tmle_task <- tmle_spec$make_tmle_task(df_long, node_list)

initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

updater <- tmle3_Update_survival$new(
  maxit = 5, 
  cvtmle = TRUE,
  convergence_type = "scaled_var",
  delta_epsilon = 1e-2,
  fit_method = "l2",
  use_best = TRUE,
  verbose = TRUE
)
'
updater <- tmle3_Update$new(
  constrain_step = TRUE, one_dimensional = TRUE, 
  delta_epsilon = 3e-2, verbose = TRUE, 
  convergence_type = "scaled_var", maxit = 10
)
'
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = updater)

tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)

tmle_fit_manual <- fit_tmle3(
  tmle_task, targeted_likelihood, tmle_params,
  targeted_likelihood$updater
)

psi0_0 <- tmle_fit_manual$summary$init_est
psi0_n <- tmle_fit_manual$summary$tmle_est

# TODO: check
t_surv1 <- simulated$true_surv1(k_grid)
t_surv0 <- simulated$true_surv0(k_grid)
l2_loss(psi0_0, t_surv0)
l2_loss(psi0_n, t_surv0)

