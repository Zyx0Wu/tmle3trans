context("Test - Robustness")

library(sl3)
library(tmle3)
library(tmle3trans)
library(uuid)
library(assertthat)
library(data.table)
library(future)

# setup data for test
'
set.seed(1234)
n <- 1000
W1 <- sample(1:4, n, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05))
W2 <- rnorm(n, 0.7, 1)
W3 <- rpois(n, 3)
S <- rbinom(n, 1, expit(-1 + .5 * W1 * sin(W3 + 8) + .2 * sqrt(abs(-W2^3))))
Y <- rnorm(n, 1.4 - 0.6 * W1 - 2 * W2 + 0.7 * W3, .4)
'
set.seed(4321)
n <- 1000
W1 <- sample(1:4, n, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05))
W2 <- rnorm(n, 0.7, 1)
W3 <- rpois(n, 3)
S <- rbinom(n, 1, expit(1.4 - 0.6 * W1 - 2 * W2 + 0.7 * W3))
Y <- rnorm(n, -1 + .5 * W1 * sin(W3 + 8) + .2 * sqrt(abs(-W2^3)), .4)

data <- data.table(W1,W2,W3,S,Y)
node_list <- list(W = c("W1", "W2", "W3"), S = "S", Y = "Y")
W <- data[ ,colnames(data) %in% node_list$W, with=FALSE]

data0 <- data[data[[node_list$S]] == 0, ]
data1 <- data[data[[node_list$S]] == 1, ]

### observed mean ###
YS0 <- data0[ ,colnames(data0) %in% node_list$Y, with=FALSE]
mean <- mean(as.matrix(YS0))
sd <- sd(as.matrix(YS0)) / sqrt(length(as.matrix(YS0)))

### glm fit ###
WS1 <- data1[ ,colnames(data1) %in% node_list$W, with=FALSE]
YS1 <- data1[ ,colnames(data1) %in% node_list$Y, with=FALSE]
WS0 <- data0[ ,colnames(data0) %in% node_list$W, with=FALSE]

fit_s <- glm(paste(node_list$S, "~", paste(node_list$W, collapse = " + ")),
             family = binomial(), data = cbind(W, S))
fit_y <- glm(paste(node_list$Y, "~", paste(node_list$W, collapse = " + ")),
             family = gaussian(), data = cbind(WS1, YS1))
beta_s_cov <- as.matrix(vcov(fit_s))
beta_y_cov <- as.matrix(vcov(fit_y))

### 1. Naive ###
# note: not Plug-In
IS0EYW <- predict(fit_y, newdata = WS0, type = 'response')

naive_psi <- mean(IS0EYW)
# delta method:
#naive_se <- sqrt(deltaMeanOLS(WS0, IS0EYW, beta_y_cov))
# by EIC:
naive_se <- sd(IS0EYW)/sqrt(length(IS0EYW))
naive_CI95 <- wald_ci(naive_psi, naive_se)

### 2. IPTW ###
pS1W <- predict(fit_s, newdata = W, type = 'response')
IS1 <- S == 1

pS0W <- 1 - pS1W
pS0 <- mean(S == 0)

iptw_psis <- IS1/prob_clip(pS1W) * pS0W/prob_clip(pS0) * Y

iptw_psi <- mean(iptw_psis)
# delta method:
#TODO: IPTW delta method
#iptw_se <- ???
# by EIC:
iptw_se <- sd(iptw_psis)/sqrt(n)
iptw_CI95 <- wald_ci(iptw_psi, iptw_se)

### 3. TML ###
tmle_spec <- tmle_AOT(1, 0)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define learners
lrnr_glm <- make_learner(Lrnr_glm_fast)
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
init_psi <- tmle_summary$init_est
tmle_psi <- tmle_summary$tmle_est
tmle_se <- tmle_summary$se
tmle_epsilon <- updater$epsilons[[1]]$Y

tmle_CI95 <- wald_ci(tmle_psi, tmle_se)

