library(simcausal)
library(testthat)
library(data.table)
library(ggplot2)

library(sl3)
library(tmle3)
library(tmle3tr)

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
n_sim <- 1e3
simulated <- simulate_data(n_sim = n_sim)
df <- simulated$dat
true_surv <- simulated$true_surv0

# TODO: check
while (all(df$Delta == 1)) {
  simulated <- simulate_data(n_sim = n_sim)
  df <- simulated$dat
  true_surv <- simulated$true_surv0
}

################################################################################
# tlverse
tmax <- max(df$T.tilde)
all_times <- lapply(seq_len(tmax), function(t) df_time(df, t))

df_long <- rbindlist(all_times)

node_list <- list(id ="ID", W = c("W", "W1"), S = "S", T = "T.tilde", D = "Delta", 
                  time = "t", F = "Failed", C = "Censored", pre_failure = "pre_failure")

lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_gam <- make_learner(Lrnr_gam)
# lrnr_earth <- make_learner(Lrnr_earth)
# lrnr_sl <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm, lrnr_gam, lrnr_earth))
lrnr_sl <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm, lrnr_gam))
learner_list <- list(A = lrnr_sl, F = lrnr_sl, C = lrnr_sl)

# TODO: check
var_types <- list(T = Variable_Type$new("continuous"), t = Variable_Type$new("continuous"), 
                  D = Variable_Type$new("binomial"))
tmle_spec <- tmle_SOT(treatment_level = 1, control_level = 0, variable_types = var_types)
tmle_task <- tmle_spec$make_tmle_task(df_long, node_list)
F_task <- tmle_task$get_regression_task("F",is_time_variant = TRUE, drop_censored = TRUE)
F_task$nrow
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# TODO: check
# <<<<<<< HEAD
# # up <- tmle3_Update_survival$new(maxit = 1e2, clipping = 1e-1 / sqrt(n_sim))
# # up <- tmle3_Update_survival$new(
# #     one_dimensional = TRUE, constrain_step = TRUE,
# #     maxit = 1e2, cvtmle = TRUE,
# #     convergence_type = "sample_size",
# #     delta_epsilon = 1e-2,
# #     fit_method = "classic"
# #   )
up <- tmle3_Update_survival$new(
  # TODO: check
  # one_dimensional = TRUE, constrain_step = TRUE,
  maxit = 3e1, cvtmle = TRUE,
  convergence_type = "sample_size",
  # delta_epsilon = 1e-2,
  fit_method = "l2",
  clipping = 1e-2
)
# =======
#up <- tmle3_Update_survival$new(maxit = 2e1, clipping = 1e-2)
up <- tmle3_Update$new(constrain_step = TRUE, one_dimensional = TRUE, 
                       delta_epsilon = 3e-2, verbose = TRUE,
                       convergence_type = "scaled_var", maxit = 10)
# up <- tmle3_Update$new(verbose = TRUE)
# debugonce(up$fit_submodel)
# debugonce(up$generate_submodel_data)
# debugonce(up$apply_submodel)
# >>>>>>> 74c33415cee164d7c7513d05c46747ef6792f661
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = up)
# targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)

ps <- tmle_params[[1]]
# debugonce(ps$estimates)
max(abs(colMeans(ps$estimates(tmle_task)$IC)))
# HA <- ps$clever_covariates(tmle_task)$N

# TODO: initial
ps <- tmle_params[[1]]
cf_task <- ps$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]
pN1 <- ps$observed_likelihood$get_likelihoods(cf_task, "N")
time <- tmle_task$time
id <- tmle_task$id
pN1_mat <- ps$long_to_mat(pN1,id,time)
SN1_mat <- ps$hm_to_sm(pN1_mat)
psi_tl_initial <- colMeans(SN1_mat)
psi_tl_initial <- c(1, psi_tl_initial[seq(1, length(psi_tl_initial) - 1)])


# tlverse update process
# debugonce(tmle_params[[1]]$estimates)
# tmle_params[[1]]$estimates(tmle_task)
tmle_fit_manual <- fit_tmle3(
  tmle_task, targeted_likelihood, tmle_params,
  targeted_likelihood$updater
)

# cf_task <- ps$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]
# ps$estimates(cf_task)
# pN1 <- ps$observed_likelihood$get_likelihoods(cf_task, "N")
# pS_N1 <- ps$hazards_to_survival(pN1, tmax)
# r_pN1 <- reshape_long_data(pN1, tmax)
# r_pS_N1 <- reshape_long_data(pS_N1, tmax)
# sum(sl_fit$density_failure_1$survival - r_pS_N1)

# TODO: check
t_surv1 <- simulated$true_surv1(k_grid - 1)
# survival_truth_1 <- survival_curve$new(t = k_grid, survival = t_surv1)

psi0_moss <- colMeans(sl_fit$density_failure_1$survival)
psi0_tl <- ps$get_psi(pS_N1, tmax)
l2_loss(psi0_moss, t_surv1)
l2_loss(psi0_tl, t_surv1)

dt <- data.table(psi0_moss, psi0_tl, t_surv1)
dt[,t:=.I]
long <- melt(dt,id="t")
ggplot(long,aes(x=t, y=value, color=variable))+geom_line()+theme_bw()

################################################################################
# moss hazard submodel
moss_hazard_l2 <- MOSS_hazard$new(
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1,
  k_grid = k_grid
)
# TODO: check
rs_moss <- moss_hazard_l2$iterate_onestep(
  method = "l2", epsilon = 1e-1 / sqrt(n_sim), max_num_interation = 1e2, verbose = FALSE
)
psi1_moss <- rs_moss$psi_n
eic_moss <- rs_moss$eic_list



rs <- tmle_fit_manual$estimates[[1]]
psi1_tl <- rs$psi
psi1_tl <- c(1, psi1_tl[seq(1, length(psi1_tl) - 1)])

# sum(psi1_moss - psi1_tl)

# TODO: check
eic_tl <- up$update(targeted_likelihood, tmle_task)

l2_loss(psi1_moss, t_surv1)
l2_loss(psi1_tl, t_surv1)

dt <- data.table(psi1_moss, psi1_tl, t_surv1)
dt[,t:=.I]
long <- melt(dt,id="t")
ggplot(long,aes(x=t, y=value, color=variable))+geom_line()+theme_bw()
