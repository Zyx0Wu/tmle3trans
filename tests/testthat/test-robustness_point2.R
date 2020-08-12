context("Test - Robustness")

library(sl3)
library(tmle3)
library(tmle3trans)
library(uuid)
library(assertthat)
library(data.table)
library(future)

gendata_trans = function(n) {
  W1 <- sample(1:4, n, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05))
  W2 <- rnorm(n, 0.7, 1)
  W3 <- rpois(n, 3)
  S <- rbinom(n, 1, expit(1.4 - 0.6 * W1 - 2 * W2 + 0.7 * W3))
  Y <- rnorm(n, -1 - .5 * W1 + .8 * W2 + .2 * W3, .4)
  data <- data.table(W1,W2,W3,S,Y)
  data
}

n <- 1000
biasS <- TRUE
biasY <- TRUE
wt <- TRUE
seed <- NULL

if (!is.null(seed)) {set.seed(seed)}

data <- gendata_trans(n)
node_list <- list(W = c("W1", "W2", "W3"), S = "S", Y = "Y")

W <- select(data, node_list$W)

data0 <- filter(data, S==0)
data1 <- filter(data, S==1)

### fit ###
WS1 <- select(data1, node_list$W)
YS1 <- data1$Y
WS0 <- select(data0, node_list$W)

s_dens <- function(w, bias=FALSE, n)
  prob_clip(expit(1.4 - 0.6 * w[1] - 2 * w[2] + 0.7 * w[3] +
                    bias * n^(-0.3)))
y_dens <- function(w, bias=FALSE, n) -1 - .5 * w[1] +
  .8 * w[2] + .2 * w[3] + bias * n^(-0.3)

### 1. Plug-In ###
Q <- apply(W, 1, y_dens, bias = biasY, n = n)
pS0 <- mean(data$S==0)

S <- data$S
Y <- data$Y

H <- (1-S)/pS0

plugin_psi <- mean(H*Q)
plugin_se <- sd(H*(Q-plugin_psi))/sqrt(n)

### 2. IPTW & S-IPTW ###
pS1W <- apply(W, 1, s_dens,n=n, bias=T)

pS0W <- 1 - pS1W

H1 <- S*pS0W/(pS0*pS1W)

iptw_psi <- mean(H1*Y)
iptw_se <- sd(H1*(Y-iptw_psi))/sqrt(n)

siptw_psi <- iptw_psi/(mean(H1))
siptw_se <- sd(H1/mean(H1)*(Y-siptw_psi))/sqrt(n)
'
### 3. TML ###

tmle_psi_init = mean(Q[S==0])

if (wt) {
  tmlefit = glm(Y[S==1]~1+offset(Q[S==1]),
                family = gaussian, weights=H1[S==1])
  Qstar = Q+ tmlefit$coefficients
} else {
  tmlefit = glm(Y[S==1]~-1+H1[S==1] + offset(Q[S==1]),
                family = gaussian)
  Qstar = tmlefit$coefficients*pS0W/(mean(1-S)*pS1W)+Q
}

D = S*pS0W/(mean(1-S)*pS1W)*(Y-Qstar) +
  (1-S)/mean(1-S)*(Qstar-mean(Qstar[S==0]))

mean(D)

tmle_psi = mean(Qstar[S==0])
tmle_psi_init
tmle_psi
tmle_se <- sd(D)/sqrt(n)

tmle_CI95 <- c(tmle_psi, tmle_psi - 1.96*tmle_se, tmle_psi + 1.96*tmle_se)
tmle_CI95
'
### 4. TML3 ###

tmle3_spec <- tmle_AOT(1, 0)

# define data
tmle3_task <- tmle3_spec$make_tmle_task(data, node_list)

# define likelihood
g_dens <- function(task) apply(task$get_node("covariates"), 1, s_dens, bias=TRUE, n)
Q_mean <- function(task) apply(task$get_node("covariates"), 1, y_dens, bias=TRUE, n)

factor_list <- list(
  define_lf(LF_emp, "W"),
  define_lf(LF_known, "S", density_fun = g_dens),
  define_lf(LF_known, "Y", mean_fun = Q_mean, type = "mean")
)

# estimate likelihood
initial_likelihood <- Likelihood$new(factor_list)$train(tmle3_task)

# define update method (submodel + loss function)
updater <- tmle3_Update$new(maxit = 10)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)

# define parameter
tmle3_param <- tmle3_spec$make_params(tmle3_task, targeted_likelihood)

# fit
tmle3_fit <- fit_tmle3(tmle3_task, targeted_likelihood, tmle3_param, updater)

# extract results
tmle3_summary <- tmle3_fit$summary
sl_psi <- tmle3_summary$init_est
tmle3_psi <- tmle3_summary$tmle_est
tmle3_se <- tmle3_summary$se

