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
library(dplyr)
library(gridExtra)
library(grid)
library(gtable)

gendata_trans = function(n) {
  W1 <- sample(1:4, n, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05))
  W2 <- rnorm(n, 0.7, 1)
  W3 <- rpois(n, 3)
  S <- rbinom(n, 1, expit(1.4 - 0.6 * W1 - 2 * W2 + 0.7 * W3))
  Y <- rnorm(n, -1 - .5 * W1 + .8 * W2 + .2 * W3, .4)
  data <- data.table(W1,W2,W3,S,Y)
  data
}

### fit ###
s_dens <- function(w, bias=FALSE) expit(1.4 - 0.6 * w[1] - 2 * w[2] + 0.7 * w[3] + bias * n^(-0.4))
y_dens <- function(w, bias=FALSE) -1 - .5 * w[1] + .8 * w[2] + .2 * w[3] + bias * n^(-0.3)

fit_s <- function(W, bias=TRUE) apply(W, 1, s_dens, bias=bias)
fit_y <- function(W, bias=TRUE) apply(W, 1, y_dens, bias=bias)

### data ###
set.seed(1234)
n <- 1e5
data <- gendata_trans(n)
node_list <- list(W = c("W1", "W2", "W3"), S = "S", Y = "Y")

W <- data[ ,colnames(data) %in% node_list$W, with=FALSE]

### observed true parameter ###
IS0 <- data$S == 0
pS0 <- mean(data$S == 0)

means <- IS0/pS0 * data$Y

mean <- mean(means)
sd <- sd(means) / sqrt(n)

### observed true IC ###
pS1W <- fit_s(W, bias=FALSE)
EYW <- fit_y(W, bias=FALSE)
IS1 <- data$S == 1
IS0 <- data$S == 0
pS0W <- 1 - pS1W
pS0 <- mean(data$S == 0)

D = IS1/pS1W*pS0W/pS0 * (data$Y-EYW) + IS0/pS0 * (EYW-mean)

# setup data for test
seed_fit = function(seed, n) {
  set.seed(seed)
  
  data = gendata_trans(n)
  node_list <- list(W = c("W1", "W2", "W3"), S = "S", Y = "Y")
  
  W <- data[ ,colnames(data) %in% node_list$W, with=FALSE]
  
  ### 1. Naive ###
  EYW <- fit_y(W)
  IS0 <- data$S == 0
  pS0 <- mean(data$S == 0)
  
  naive_psis <- IS0/pS0 * EYW
  
  naive_psi <- mean(naive_psis)
  naive_se <- sd(naive_psis)/sqrt(n)
  naive_CI95 <- wald_ci(naive_psi, naive_se)
  
  ### 2. IPTW ###
  pS1W <- fit_s(W)
  IS1 <- data$S == 1
  pS0W <- 1 - pS1W
  pS0 <- mean(data$S == 0)
  
  Hs <- IS1/prob_clip(pS1W) * pS0W/pS0
  iptw_psis <- Hs/mean(Hs) * data$Y
  
  iptw_psi <- mean(iptw_psis)
  iptw_se <- sd(iptw_psis)/sqrt(n)
  iptw_CI95 <- wald_ci(iptw_psi, iptw_se)
  
  ### 3. TML ###
  tmle_spec <- tmle_AOT(1, 0)
  
  # define data
  tmle_task <- tmle_spec$make_tmle_task(data, node_list)
  
  # define likelihood
  g_dens <- function(task) apply(task$get_node("covariates"), 1, s_dens, bias=TRUE)
  Q_mean <- function(task) apply(task$get_node("covariates"), 1, y_dens, bias=TRUE)
  
  factor_list <- list(
    define_lf(LF_emp, "W"),
    define_lf(LF_known, "S", density_fun = g_dens),
    define_lf(LF_known, "Y", mean_fun = Q_mean, type = "mean")
  )
  
  # estimate likelihood
  initial_likelihood <- Likelihood$new(factor_list)$train(tmle_task)
  
  # define update method (submodel + loss function)
  updater <- tmle3_Update$new(maxit = 10)
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
  
  return(c(naive_psi, iptw_psi, tmle_psi, naive_se, iptw_se, tmle_se))
}

reps <- 500
obs <- 10000
naive_psis <- c()
iptw_psis <- c()
tmle_psis <- c()
naive_ses <- c()
iptw_ses <- c()
tmle_ses <- c()
for (i in 1:reps) {
  fits <- seed_fit(1234+i, obs)
  naive_psis <- c(naive_psis, fits[1])
  iptw_psis <- c(iptw_psis, fits[2])
  tmle_psis <- c(tmle_psis, fits[3])
  naive_ses <- c(naive_ses, fits[4])
  iptw_ses <- c(iptw_ses, fits[5])
  tmle_ses <- c(tmle_ses, fits[6])
}
diff_naive <- (naive_psis - mean)^2
diff_iptw <- (iptw_psis - mean)^2
diff_tmle <- (tmle_psis - mean)^2
naive_CI95 <- wald_ci(naive_psis, naive_ses)
iptw_CI95 <- wald_ci(iptw_psis, iptw_ses)
tmle_CI95 <- wald_ci(tmle_psis, tmle_ses)
coverage_naive <- pnorm(naive_CI95[,2], mean = mean, sd = sd) - 
  pnorm(naive_CI95[,1], mean = mean, sd = sd)
coverage_iptw <- pnorm(iptw_CI95[,2], mean = mean, sd = sd) - 
  pnorm(iptw_CI95[,1], mean = mean, sd = sd)
coverage_tmle <- pnorm(tmle_CI95[,2], mean = mean, sd = sd) - 
  pnorm(tmle_CI95[,1], mean = mean, sd = sd)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# estimator
dat_psis <- data.frame(Method=rep(c("Naive", "IPTW", "TML"), each=reps), 
                       Estimator=c(naive_psis, iptw_psis, tmle_psis))

mu <- ddply(dat_psis, "Method", summarise, grp.mean=mean(Estimator))
plt_hist <- ggplot(dat_psis, aes(x=Estimator, color=Method)) +
  geom_histogram(fill="white", position="dodge") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Method),
             linetype="dashed") +
  geom_vline(aes(xintercept=mean, color='Truth'))
theme(legend.position="top")
ggsave("point2_robust_estimator_histogram.png", 
       plot = plt_hist, path = "../plot")

# loss
dat_diff <- data.frame(Simulation_No.=rep(1:reps, 3), 
                       Method=rep(c("Naive", "IPTW", "TML"), each=reps), 
                       Loss_l2=c(diff_naive, diff_iptw, diff_tmle))

mse <- ddply(dat_diff, "Method", summarise, grp.mean=mean(Loss_l2))
plt_plot <- ggplot(data=dat_diff, aes(x=Simulation_No., y=Loss_l2, group=Method)) +
  geom_point(aes(color=Method)) +
  geom_hline(data=mse, aes(yintercept=grp.mean, color=Method),
             linetype="dashed")
ggsave("point2_robust_loss_point.png", 
       plot = plt_plot, path = "../plot")

# summary
truth <- mean
sample_mean <- c(mean(naive_psis), mean(iptw_psis), mean(tmle_psis))
sample_sd <- c(sd(naive_psis), sd(iptw_psis), sd(tmle_psis))
coverage <- c(mean(coverage_naive), mean(coverage_iptw), mean(coverage_tmle))
dat_summary <- data.frame(Method=c("Naive", "IPTW", "TML"), 
                          Truth=truth, Sample_Mean=sample_mean, Sample_Sd=sample_sd, 
                          MSE=mse[match(mse$Method, c("Naive", "IPTW", "TML")),]$grp.mean,
                          Coverage_95=label_percent()(coverage))

theme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.5)),
  colhead = list(fg_params=list(cex = 0.5)),
  rowhead = list(fg_params=list(cex = 0.5)))

g <- tableGrob(dat_summary, rows = NULL, theme = theme)
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))
grid.draw(g)

# efficiency
mean(tmle_ses) * sqrt(obs) / sd(D)

