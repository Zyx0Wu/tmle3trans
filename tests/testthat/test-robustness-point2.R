context("Test - Robustness")

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
W1 <- sample(1:4, n, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05))
W2 <- rnorm(n, 0.7, 1)
W3 <- rpois(n, 3)
S <- rbinom(n, 1, prob_clip(expit(1.4 - 0.6 * W1 - 2 * W2 + 0.7 * W3)))
Y <- rnorm(n, -1 - .5 * W1 + .8 * W2 + .2 * W3, .4)

data <- data.table(W1,W2,W3,S,Y)
node_list <- list(W = c("W1", "W2", "W3"), S = "S", Y = "Y")
W <- data[ ,colnames(data) %in% node_list$W, with=FALSE]

data0 <- data[data[[node_list$S]] == 0, ]
data1 <- data[data[[node_list$S]] == 1, ]

### observed mean ###
YS0 <- data0[ ,colnames(data0) %in% node_list$Y, with=FALSE]
mean <- mean(as.matrix(YS0))
sd <- sd(as.matrix(YS0)) / sqrt(length(as.matrix(YS0)))

### fit ###
WS1 <- data1[ ,colnames(data1) %in% node_list$W, with=FALSE]
YS1 <- data1[ ,colnames(data1) %in% node_list$Y, with=FALSE]
WS0 <- data0[ ,colnames(data0) %in% node_list$W, with=FALSE]

s_dens <- function(w, bias=FALSE) prob_clip(expit(1.4 - 0.6 * w[1] - 2 * w[2] + 0.7 * w[3] + bias * n^(-0.3)))
y_dens <- function(w, bias=FALSE) -1 - .5 * w[1] + .8 * w[2] + .2 * w[3] + bias * n^(-0.3)

fit_s <- function(W) apply(W, 1, s_dens, bias=FALSE)
fit_y <- function(W) apply(W, 1, y_dens, bias=TRUE)
#fit_s <- function(W) apply(W, 1, s_dens, bias=TRUE)
#fit_y <- function(W) apply(W, 1, y_dens, bias=FALSE)

### 1. Naive ###
# note: not Plug-In
IS0EYW <- fit_y(WS0)

naive_psi <- mean(IS0EYW)
naive_se <- sd(IS0EYW)
naive_CI95 <- sprintf("(%f, %f)", naive_psi - 1.96*naive_se, naive_psi + 1.96*naive_se)

### 2. IPTW ###
pS1W <- fit_s(W)
IS1 <- S == 1

pS0W <- 1 - pS1W
pS0 <- mean(pS0W)

iptw_psi <- mean(IS1/prob_clip(pS1W) * pS0W/prob_clip(pS0) * Y)
iptw_se <- sd(IS1/prob_clip(pS1W) * pS0W/prob_clip(pS0) * Y)
iptw_CI95 <- sprintf("(%f, %f)", iptw_psi - 1.96*iptw_se, iptw_psi + 1.96*iptw_se)

### 3. TML ###
tmle_spec <- tmle_AOT(1, 0)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
g_dens <- function(task) apply(task$get_node("covariates"), 1, s_dens, bias=FALSE)
Q_mean <- function(task) apply(task$get_node("covariates"), 1, y_dens, bias=TRUE)
#g_dens <- function(task) apply(task$get_tmle_node("W"), 1, s_dens, bias=TRUE)
#Q_mean <- function(task) apply(task$get_tmle_node("W"), 1, y_dens, bias=FALSE)

factor_list <- list(
  define_lf(LF_emp, "W"),
  define_lf(LF_known, "S", density_fun = g_dens),
  define_lf(LF_known, "Y", mean_fun = Q_mean, type = "mean")
)

# estimate likelihood
initial_likelihood <- Likelihood$new(factor_list)$train(tmle_task)

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

