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
    node("W1", distr = "rbinom", size = 1, prob = .5) +
    node("W", distr = "runif", min = 0, max = 1.5) +
    node("S", distr = "rbinom", size = 1, prob = .15 + .5 * as.numeric(W > .75)) +
    node("Trexp", distr = "rexp", rate = 1 + .7 * W^2 - .8 * S) +
    node("Cweib", distr = "rweibull", shape = 1 + .5 * W, scale = 75) +
    node("T", distr = "rconst", const = round(Trexp * 2)) +
    node("C", distr = "rconst", const = round(Cweib * 2)) +
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
  true_surv <- function(q, W, S) sapply(W, function(w) {
    1 - pexp(q, rate = 1 + .7 * w^2 - .8 * S)
  })
  # input: vector q. mean(S(q|S=s,W)), average out W. loop over q
  true_surv <- function(q_grid, surv_fn, s) {
    W_grid <- seq(0, 1.5, .01)
    survout <- numeric()
    for (q in q_grid) survout <- c(survout, mean(surv_fn(q = q / 2, W = W_grid, S = s)))
    return(survout)
  }
  truth_surv0 <- function(q) true_surv(q_grid = q, surv_fn = true_surv_one, s = 0)
  return(list(dat = dat, true_surv0 = truth_surv0))
}

set.seed(1234)
n_sim <- 1e2
simulated <- simulate_data(n_sim = n_sim)
df <- simulated$dat
true_surv <- simulated$true_surv0

# TODO: check
while (all(df$Delta == 1)) {
  simulated <- simulate_data(n_sim = n_sim)
  df <- simulated$dat
  true_surv <- simulated$true_surv0
}

###########
# tlverse #
###########
tmax <- max(df$T.tilde)
all_times <- lapply(seq_len(tmax), function(t) df_time(df, t))

df_long <- rbindlist(all_times)

node_list <- list(id ="ID", W = c("W", "W1"), S = "S", T = "T.tilde", D = "Delta", 
                  time = "t", F = "Failed", C = "Censored", pre_failure = "pre_failure")

### sl3 ###
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_gam <- make_learner(Lrnr_gam)
#lrnr_earth <- make_learner(Lrnr_earth)
#lrnr_sl <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm, lrnr_gam, lrnr_earth))
lrnr_sl <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm, lrnr_gam))
learner_list <- list(S = lrnr_sl, A = lrnr_sl, F = lrnr_sl, C = lrnr_sl)

### tmle3 ###
# TODO: check
#var_types <- list(T = Variable_Type$new("continuous"), t = Variable_Type$new("continuous"), 
#                  D = Variable_Type$new("binomial"))
#tmle_spec <- tmle_SOT(1, 0, variable_types = var_types)
tmle_spec <- tmle_SOT(1, 0)
tmle_task <- tmle_spec$make_tmle_task(df_long, node_list)

initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

updater <- tmle3_Update_survival$new(
  maxit = 3e1, 
  cvtmle = TRUE,
  convergence_type = "scaled_var",
  delta_epsilon = 1e-2,
  fit_method = "l2",
  use_best = TRUE,
  verbose=FALSE
)
#updater <- tmle3_Update$new(
#  constrain_step = TRUE, one_dimensional = TRUE, 
#  delta_epsilon = 3e-2, verbose = TRUE, 
#  convergence_type = "scaled_var", maxit = 10)

targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = updater)

tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)

tmle_fit_manual <- fit_tmle3(
  tmle_task, targeted_likelihood, tmle_params,
  targeted_likelihood$updater
)

tmle_fit_manual

