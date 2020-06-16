context("Test - cpp data")

library(sl3)
library(tmle3)
library(tmle3tr)
library(uuid)
library(assertthat)
library(data.table)
library(future)

# note: S <=> A, site variable

# setup data for test
data(cpp)
data <- as.data.table(cpp)
data$parity01 <- as.numeric(data$parity > 0)
data$parity01_fac <- factor(data$parity01)
data$haz01 <- as.numeric(data$haz > 0)
data[is.na(data)] <- 0
node_list <- list(
  W = c("sexn"),
  A = "parity01",
  Y = "haz01"
)

qlib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

glib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

logit_metalearner <- make_learner(
  Lrnr_solnp, metalearner_logistic_binomial,
  loss_loglik_binomial
)
Q_learner <- make_learner(Lrnr_sl, qlib, logit_metalearner)
g_learner <- make_learner(Lrnr_sl, glib, logit_metalearner)
learner_list <- list(Y = Q_learner, A = g_learner)

tmle_spec <- tmle_TR(1, 0)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)
tmle_task_onsite <- tmle_spec$make_tmle_task(data_onsite, node_list)

# define likelihood
tmle_spec$make_initial_likelihood(tmle_task, learner_list, initial=TRUE)
initial_likelihood_onsite <- tmle_spec$make_initial_likelihood(tmle_task_onsite, learner_list)

# define update method (submodel + loss function)
# disable cvtmle for this test to compare with tmle package
updater <- tmle3_Update$new(cvtmle = FALSE, convergence_type = "sample_size")
targeted_likelihood_onsite <- Targeted_Likelihood$new(initial_likelihood_onsite, updater)

# define parameter
tmle_params <- tmle_spec$make_params(targeted_likelihood_onsite)
updater$tmle_params <- tmle_params
TR <- tmle_params[[1]]

# fit
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood_onsite, list(TR), updater)

# extract results
tmle_summary <- tmle_fit$summary
tmle_psi <- tmle_summary$tmle_est
tmle_se <- tmle_summary$se
tmle_epsilon <- updater$epsilons[[1]]$Y
