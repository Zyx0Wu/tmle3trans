context("Test - Simulation")

library(sl3)
library(tmle3)
library(tmle3tr)
library(uuid)
library(assertthat)
library(data.table)
library(future)

# note: S <=> A, site variable

# generate data
set.seed(1234)
n <- 1000
W <- sample(1:4, n, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05))
A <- rbinom(n, 1, expit(1.4 - 0.6 * W))
#Y <- rnorm(n, 1 + .5 * W * sin(W + 8) + .2 * sqrt(abs(-W^3 + exp(W))), .4)
Y <- rnorm(n, 2 + .5 * W, 1)

data <- data.table(W, A, Y)
data_onsite <- data[A == 1,]
node_list <- list(W = "W", A = "A", Y = "Y")

tmle_spec <- tmle_TR(1, 0)
tmle_task <- tmle_spec$make_tmle_task(data, node_list)
tmle_task_onsite <- tmle_spec$make_tmle_task(data_onsite, node_list)

lrnr_glm <- make_learner(Lrnr_glm_fast)
learner_list <- list(Y = lrnr_glm, A = lrnr_glm)

factor_list <- list(
  define_lf(LF_emp, "W"),
  define_lf(LF_fit, "A", learner = learner_list[["A"]]),
  define_lf(LF_fit, "Y", learner = learner_list[["Y"]], type = "mean")
)

initial_likelihood_all <- Likelihood$new(factor_list)$train(tmle_task)
tmle_spec$set_initial_likelihood(initial_likelihood_all)
initial_likelihood_onsite <- Likelihood$new(factor_list)$train(tmle_task)

### 1. non-parametric estimator ###
pAW <- initial_likelihood_all$get_likelihood(tmle_task, "A")
Y1W <- initial_likelihood_onsite <- Likelihood$new(factor_list)$train(tmle_task)


### 2. TMLE ###
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
