context("Test - Numeric Check")

library(sl3)
library(tmle3)
library(tmle3tr)
library(uuid)
library(assertthat)
library(testthat)
library(data.table)
library(future)

# generate data
set.seed(4321)
n <- 10000
W <- sample(1:4, n, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05))
S <- rbinom(n, 1, expit(1.4 - 0.6 * W))
Y <- rnorm(n, 2 + .5 * W, 1)
'
# truth
w = 1:4
pw = c(0.1, 0.2, 0.65, 0.05)
pw0 = (1 - expit(1.4 - 0.6 * w)) * pw
pw0 = pw0 / sum(pw0)
Ey0 = sum((2 + .5 * w) * pw0)
'

## TMLE should not update in this case and
## should numerically match the non-parametric estimator

### 1. non-parametric estimator ###
psi_0 = function(b) mean(W == b & S == 0)
psi_1 = function(b) mean(W == b & S == 1)
psi_2 = function(b) mean((W == b & S == 1) * Y)
psi_s = mean(S == 0)

f_n = sum(sapply(levels(factor(W)),
               function(b) psi_0(b)/psi_s * psi_2(b)/ifelse(psi_1(b), psi_1(b), 1)))
D_n = rowSums(sapply(levels(factor(W)),
                     function(b) (psi_0(b) * (W == b & S == 1))/(psi_s * ifelse(psi_1(b), psi_1(b), 1)) * (Y - psi_2(b)/ifelse(psi_1(b), psi_1(b), 1)) +
                       (W == b & S == 0)/psi_s * (psi_2(b)/ifelse(psi_1(b), psi_1(b), 1) - f_n)))
s_n = sd(D_n) / sqrt(n)

### 2. TMLE ###
data <- data.table(W=factor(W), S=factor(S), Y)
node_list <- list(W = "W", S = "S", Y = "Y")

tmle_spec <- tmle_AOT(1, 0)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define learners
lrnr_glm <- make_learner(Lrnr_glm)
learner_list <- list(Y = lrnr_glm, S = lrnr_glm)

# estimate likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# define update method (submodel + loss function)
updater <- tmle3_Update$new()
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)

# define parameter
tmle_param <- tmle_spec$make_params(tmle_task, targeted_likelihood)

# fit
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_param, updater)

# extract results
tmle_summary <- tmle_fit$summary
tmle_psi <- tmle_summary$tmle_est
tmle_se <- tmle_summary$se
tmle_epsilon <- updater$epsilons[[1]]$Y

### tests ###
tol <- 1 / sqrt(n)
test_that("psi results match", {
  expect_equal(tmle_psi, f_n, tol = tol)
})
test_that("se results match", {
  expect_equal(tmle_se, s_n, tol = tol)
})
'
Ey0 >= tmle_summary$lower && Ey0 <= tmle_summary$upper
'

