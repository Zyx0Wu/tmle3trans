context("Test - Robustness")

library(sl3)
library(tmle3)
library(tmle3tr)
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
set.seed(1234)
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
Y0 <- data0[ ,colnames(data0) %in% node_list$Y, with=FALSE]
mean <- mean(as.matrix(Y0))
std <- sd(as.matrix(Y0)) / sqrt(length(as.matrix(Y0)))

### 1. naive ###
W1 <- data1[ ,colnames(data1) %in% node_list$W, with=FALSE]
Y1 <- data1[ ,colnames(data1) %in% node_list$Y, with=FALSE]
W0 <- data0[ ,colnames(data0) %in% node_list$W, with=FALSE]

fit <- glm(paste(node_list$Y, "~", paste(node_list$W, collapse = " + ")),
           data = cbind(W1, Y1))
est <- predict(fit, newdata = W0, type = 'response', se.fit = TRUE)
psis <- predict(fit, newdata = W0, type = 'response', se.fit = TRUE)$fit
ses <- predict(fit, newdata = W0, type = 'response', se.fit = TRUE)$se.fit

psi <- mean(psis)
se <- sqrt(mean(ses^2)) / sqrt(length(ses))
CI95 <- sprintf("(%f, %f)", psi - 1.96*se, psi + 1.96*se)

### 2. TMLE ###
tmle_spec <- tmle_AET(1, 0)

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
sl_psi <- tmle_summary$init_est
tmle_psi <- tmle_summary$tmle_est
tmle_se <- tmle_summary$se
tmle_epsilon <- updater$epsilons[[1]]$Y

tmle_CI95 <- sprintf("(%f, %f)", tmle_psi - 1.96*tmle_se, tmle_psi + 1.96*tmle_se)

