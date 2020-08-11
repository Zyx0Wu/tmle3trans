library(sl3)
library(tmle3)
library(tmle3trans)

library(uuid)
library(future)
library(assertthat)
library(data.table)
library(ggplot2)
library(scales)
library(plyr)
library(gridExtra)
library(grid)
library(gtable)

set.seed(1234)
n <- 1e5
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

# setup data for test
seed_fit = function(seed) {
  set.seed(seed)
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
  
  ### 1. Plug-In ###
  WS1 <- data1[ ,colnames(data1) %in% node_list$W, with=FALSE]
  YS1 <- data1[ ,colnames(data1) %in% node_list$Y, with=FALSE]
  WS0 <- data0[ ,colnames(data0) %in% node_list$W, with=FALSE]
  
  fit_y <- glm(paste(node_list$Y, "~", paste(node_list$W, collapse = " + ")),
               family = gaussian(), data = cbind(WS1, YS1))
  beta_y_cov <- as.matrix(vcov(fit_y))
  
  plugin_psis <- predict(fit_y, newdata = WS0, type = 'response')
  
  plugin_psi <- mean(plugin_psis)
  # delta method:
  plugin_se <- sqrt(deltaMeanOLS(WS0, plugin_psis, beta_y_cov))
  # empirical:
  #plugin_se <- sd(plugin_psis)/sqrt(length(plugin_psis))
  #plugin_CI95 <- wald_ci(plugin_psi, plugin_se)
  
  ### 2. Non-parametric ###
  psi_0 = function(b) mean(W1 == b[1] & W2 == b[2] & W3 == b[3] & S == 0)
  psi_1 = function(b) mean(W1 == b[1] & W2 == b[2] & W3 == b[3] & S == 1)
  psi_2 = function(b) mean((W1 == b[1] & W2 == b[2] & W3 == b[3] & S == 1) * Y)
  psi_s = mean(S == 0)
  
  nonpar_psi = sum(apply(expand.grid(levels(factor(W1)), levels(factor(W2)), levels(factor(W3))), 1,
                         function(b) psi_0(b)/psi_s * psi_2(b)/ifelse(psi_1(b), psi_1(b), 1)))
  # influence curve:
  nonpar_ses = rowSums(apply(expand.grid(levels(factor(W1)), levels(factor(W2)), levels(factor(W3))), 1,
                             function(b) (psi_0(b) * (W1 == b[1] & W2 == b[2] & W3 == b[3] & S == 1))/(psi_s * ifelse(psi_1(b), psi_1(b), 1)) * 
                               (Y - psi_2(b)/ifelse(psi_1(b), psi_1(b), 1)) +
                               (W1 == b[1] & W2 == b[2] & W3 == b[3] & S == 0)/psi_s * 
                               (psi_2(b)/ifelse(psi_1(b), psi_1(b), 1) - nonpar_psi)))
  nonpar_se = sd(nonpar_ses) / sqrt(n)
  
  nonpar_CI95 <- wald_ci(nonpar_psi, nonpar_se)
  
  ### 3. SL + TML ###
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
  
  sl_psis <- initial_likelihood$get_likelihood(tmle_task, "Y")[S==0]
  sl_psi <- mean(sl_psis)
  sl_se <- sd(sl_psis)/sqrt(length(sl_psis))
  
  sl_CI95 <- wald_ci(sl_psi, sl_se)
  
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
  
  tmle_CI95 <- wald_ci(tmle_psi, tmle_se)
  
  # loss
  l2_diff_plugin <- (plugin_psi - mean)^2
  l2_diff_nonpar <- (nonpar_psi - mean)^2
  l2_diff_sl <- (sl_psi - mean)^2
  l2_diff_tl <- (tmle_psi - mean)^2
  
  return(c(plugin_psi, nonpar_psi, sl_psi, tmle_psi, 
           l2_diff_plugin, l2_diff_nonpar, l2_diff_sl, l2_diff_tl))
}

reps <- 100
plugin_psis <- c()
nonpar_psis <- c()
sl_psis <- c()
tl_psis <- c()
diff_plugin <- c()
diff_nonpar <- c()
diff_sl <- c()
diff_tl <- c()
for (i in 1:reps) {
  fits <- seed_fit(i)
  plugin_psis <- c(plugin_psis, fits[1])
  nonpar_psis <- c(nonpar_psis, fits[2])
  sl_psis <- c(sl_psis, fits[3])
  tl_psis <- c(tl_psis, fits[4])
  diff_plugin <- c(diff_plugin, fits[5])
  diff_nonpar <- c(diff_nonpar, fits[6])
  diff_sl <- c(diff_sl, fits[7])
  diff_tl <- c(diff_tl, fits[8])
}

# estimator
dat_psis <- data.frame(Method=rep(c("NP", "TL"), each=reps), 
                       Estimator=c(nonpar_psis, tl_psis))

mu <- ddply(dat_psis, "Method", summarise, grp.mean=mean(Estimator))
plt_hist <- ggplot(dat_psis, aes(x=Estimator, color=Method)) +
  geom_histogram(fill="white", position="dodge") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Method),
             linetype="dashed") +
  geom_vline(aes(xintercept=mean, color='Truth'))
theme(legend.position="top")

# loss
dat_diff <- data.frame(Simulation_No.=rep(1:reps, 2), 
                       Method=rep(c("NP", "TL"), each=reps), 
                       Loss_l2=c(diff_nonpar, diff_tl))

mse <- ddply(dat_diff, "Method", summarise, grp.mean=mean(Loss_l2))
plt_plot <- ggplot(data=dat_diff, aes(x=Simulation_No., y=Loss_l2, group=Method)) +
  geom_point(aes(color=Method)) +
  geom_hline(data=mse, aes(yintercept=grp.mean, color=Method),
             linetype="dashed")

# summary
truth <- mean
sample_mean <- c(mean(nonpar_psis), mean(tl_psis))
sample_sd <- c(sd(nonpar_psis), sd(tl_psis))
mse <- c(mean(diff_nonpar), mean(diff_tl))
dat_summary <- data.frame(Method=c("NP", "TL"), Truth=truth, 
                          Sample_Mean=sample_mean, Sample_Sd=sample_sd,
                          MSE=mse)

g <- tableGrob(dat_summary, rows = NULL)
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))
grid.draw(g)

