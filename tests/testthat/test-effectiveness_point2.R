context("Test - Effectiveness")

library(sl3)
library(tmle3)
library(tmle3trans)
library(uuid)
library(assertthat)
library(data.table)
library(future)

# setup data for test
set.seed(1234)
n <- 1000
W1 <- sample(2:4, n, replace=TRUE, prob=c(0.3, 0.65, 0.05))
W2 <- sample(c(0.2, 0.9), n, replace=TRUE, prob=c(0.2, 0.8))
W3 <- sample(5:6, n, replace=TRUE, prob=c(0.55, 0.45))
#S <- rbinom(n, 1, expit(1.4 - 0.6 * W1 - 2 * W2 + 0.7 * W3))
S <- rbinom(n, 1, expit(1.4 / acosh(2 * W1) - 0.6 * W2 * (pmax(0, cos(W2 + 2)) + 1) * (W3 - 2.5)^2))
Y <- rnorm(n, -1 + .5 * W1 * sin(W3 + 8) + .2 * sqrt(abs(-W2^3 + exp(W2/(W3-3.5)))), .4)

data <- data.table(W1,W2,W3,S,Y)
node_list <- list(W = c("W1", "W2", "W3"), S = "S", Y = "Y")

data0 <- data[data[[node_list$S]] == 0, ]
data1 <- data[data[[node_list$S]] == 1, ]

### observed mean ###
YS0 <- data0[ ,colnames(data0) %in% node_list$Y, with=FALSE]
mean <- mean(as.matrix(YS0))
std <- sd(as.matrix(YS0)) / sqrt(length(as.matrix(YS0)))

### 1. Naive ###
WS1 <- data1[ ,colnames(data1) %in% node_list$W, with=FALSE]
YS1 <- data1[ ,colnames(data1) %in% node_list$Y, with=FALSE]
WS0 <- data0[ ,colnames(data0) %in% node_list$W, with=FALSE]

fit_y <- glm(paste(node_list$Y, "~", paste(node_list$W, collapse = " + ")),
             family = gaussian(), data = cbind(WS1, YS1))
beta_y_cov <- as.matrix(vcov(fit_y))

psis <- predict(fit_y, newdata = WS0, type = 'response')

psi <- mean(psis)
# delta method:
#se <- sqrt(deltaMeanOLS(WS0, psis, beta_cov))
# non-parametric:
se <- sd(psis)/sqrt(length(psis))
CI95 <- wald_ci(psi, se)

### 2. Non-parametric ###
psi_0 = function(b) mean(W1 == b[1] & W2 == b[2] & W3 == b[3] & S == 0)
psi_1 = function(b) mean(W1 == b[1] & W2 == b[2] & W3 == b[3] & S == 1)
psi_2 = function(b) mean((W1 == b[1] & W2 == b[2] & W3 == b[3] & S == 1) * Y)
psi_s = mean(S == 0)

f_n = sum(apply(expand.grid(levels(factor(W1)), levels(factor(W2)), levels(factor(W3))), 1,
                function(b) psi_0(b)/psi_s * psi_2(b)/ifelse(psi_1(b), psi_1(b), 1)))
D_n = rowSums(apply(expand.grid(levels(factor(W1)), levels(factor(W2)), levels(factor(W3))), 1,
                    function(b) (psi_0(b) * (W1 == b[1] & W2 == b[2] & W3 == b[3] & S == 1))/(psi_s * ifelse(psi_1(b), psi_1(b), 1)) * 
                      (Y - psi_2(b)/ifelse(psi_1(b), psi_1(b), 1)) +
                      (W1 == b[1] & W2 == b[2] & W3 == b[3] & S == 0)/psi_s * 
                      (psi_2(b)/ifelse(psi_1(b), psi_1(b), 1) - f_n)))
s_n = sd(D_n) / sqrt(n)

### 3. TML ###
tmle_spec <- tmle_AOT(1, 0)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define learners
qlib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast",
  "Lrnr_xgboost"
)

glib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast",
  "Lrnr_xgboost"
)

ls_metalearner <- make_learner(Lrnr_nnls)
bn_metalearner <- make_learner(
  Lrnr_solnp, metalearner_logistic_binomial,
  loss_loglik_binomial
)
Q_learner <- make_learner(Lrnr_sl, qlib, ls_metalearner)
g_learner <- make_learner(Lrnr_sl, glib, bn_metalearner)
#g_learner <- make_learner(Lrnr_glm)
learner_list <- list(Y = Q_learner, S = g_learner)

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
sl_psi <- tmle_summary$init_est
tmle_psi <- tmle_summary$tmle_est
tmle_se <- tmle_summary$se
tmle_epsilon <- updater$epsilons[[1]]$Y

tmle_CI95 <- wald_ci(tmle_psi, tmle_se)

