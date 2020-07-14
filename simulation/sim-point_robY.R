context("Test - Robustness Y")

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

# setup data for test
seed_fit = function(seed) {
  set.seed(seed)
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
  
  ### 1. Benchmark ###
  WS1 <- data1[ ,colnames(data1) %in% node_list$W, with=FALSE]
  YS1 <- data1[ ,colnames(data1) %in% node_list$Y, with=FALSE]
  WS0 <- data0[ ,colnames(data0) %in% node_list$W, with=FALSE]
  
  fit_y <- glm(paste(node_list$Y, "~", paste(node_list$W, collapse = " + ")),
               family = gaussian(), data = cbind(WS1, YS1))
  beta_y_cov <- as.matrix(vcov(fit_y))
  
  IS0EYW <- predict(fit_y, newdata = WS0, type = 'response')
  
  psi <- mean(IS0EYW)
  # delta method:
  #se <- sqrt(deltaMeanOLS(WS0, IS0EYW, beta_y_cov))
  # by EIC:
  se <- sd(IS0EYW)
  CI95 <- sprintf("(%f, %f)", psi - 1.96*se, psi + 1.96*se)
  
  ### 2. TMLE ###
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
  sl_psi <- tmle_summary$init_est
  tmle_psi <- tmle_summary$tmle_est
  tmle_se <- tmle_summary$se
  tmle_epsilon <- updater$epsilons[[1]]$Y
  
  tmle_CI95 <- sprintf("(%f, %f)", tmle_psi - 1.96*tmle_se, tmle_psi + 1.96*tmle_se)
  
  # MSE, coverage
  l2_diff_bench <- (psi - mean)^2
  l2_diff_tmle <- (tmle_psi - mean)^2
  coverage_bench <- pnorm(psi + 1.96*se, mean = mean, sd = sd) - 
    pnorm(psi - 1.96*se, mean = mean, sd = sd)
  coverage_tmle <- pnorm(tmle_psi + 1.96*tmle_se, mean = mean, sd = sd) - 
    pnorm(tmle_psi - 1.96*tmle_se, mean = mean, sd = sd)
  
  return(c(l2_diff_bench, l2_diff_tmle, coverage_bench, coverage_tmle))
}

reps <- 100
diff_bench <- c()
diff_tmle <- c()
coverage_bench <- c()
coverage_tmle <- c()
for (i in 1:reps) {
  fits <- seed_fit(i)
  diff_bench <- c(diff_bench, fits[1])
  diff_tmle <- c(diff_tmle, fits[2])
  coverage_bench <- c(coverage_bench, fits[3])
  coverage_tmle <- c(coverage_tmle, fits[4])
}

# loss
dat_diff <- data.frame(Simulation_No.=rep(1:reps, 2), 
                       Method=rep(c("Plug-In", "TML"), each=reps), 
                       Loss_l2=c(diff_bench, diff_tmle))

plt_plot <- ggplot(data=dat_diff, aes(x=Simulation_No., y=Loss_l2, group=Method)) +
  geom_point(aes(color=Method))

mu <- ddply(dat_diff, "Method", summarise, grp.mean=mean(Loss_l2))
plt_hist <- ggplot(dat_diff, aes(x=Loss_l2, color=Method)) +
  geom_histogram(fill="white", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Method),
             linetype="dashed")+
  theme(legend.position="top")

# coverage
tol <- 1e-2
coverage <- c(sum(abs(1-coverage_bench)<=tol)/reps, sum(abs(1-coverage_tmle)<=tol)/reps)
dat_coverage <- data.frame(Method=c("Plug-In", "TML"), 
                           Coverage_95=label_percent()(coverage))

g <- tableGrob(dat_coverage, rows = NULL)
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))
grid.draw(g)

