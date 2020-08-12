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

set.seed(1234)

gendata_trans = function(n) {
  W1 <- sample(1:4, n, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05))
  W2 <- rnorm(n, 0.7, 1)
  W3 <- rpois(n, 3)
  S <- rbinom(n, 1, prob_clip(expit(1.4 - 0.6 * W1 - 2 * W2 + 0.7 * W3)))
  Y <- rnorm(n, -1 - .5 * W1 + .8 * W2 + .2 * W3, .4)
  data <- data.table(W1,W2,W3,S,Y)
  data
}

# compute the truth for this DGP:
pop = gendata_trans(5e6)
pop
popS0 =  filter(pop, S==0)
popS1 =  filter(pop, S==1)
Psi0 = mean(popS0$Y)
mean(popS1$Y)
Ypop = pop$Y
Spop = pop$S

### fit ###
node_list <- list(W = c("W1", "W2", "W3"), S = "S", Y = "Y")
popW = select(pop, node_list$W)
W <- select(pop, node_list$W)

### 1. Plug-In ###
Qpop = mutate(W, Q = -1 - .5 * W1 +
                .8 * W2 + .2 * W3)$Q
PSpop = mutate(W, g=prob_clip(expit(1.4 - 0.6 * W1 - 2 * W2 +
                                      0.7 * W3)))$g
PSpop
Psi0

Dpop = Spop*(1-PSpop)/PSpop/mean(1-Spop)*(Ypop - Qpop)+
  (1-Spop)/mean(1-Spop)*(Qpop-mean(Qpop[Spop==0]))

se0 = sd(Dpop)
se0
###
###
# This is your simulation function except there is not a method
# for giving CI for plugin_psi--use delta method as you have coded
###
###

simJL = function(n, biasS=TRUE, biasY=TRUE, wt=TRUE, seed=NULL) {
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
  
  return(c(plugin_psi, iptw_psi, siptw_psi, tmle3_psi, 
           plugin_se, iptw_se, siptw_se, tmle3_se))
}

#simJL(1000,T,T,F)

#cl = makeCluster(8)
reps <- 50
obs <- 1000

fits <- lapply(1:reps, function(i) simJL(obs,T,T,F,i))
fits <- do.call(rbind, fits)
plugin_psis <- fits[,1]
iptw_psis <- fits[,2]
siptw_psis <- fits[,3]
tmle_psis <- fits[,4]
plugin_ses <- fits[,5]
iptw_ses <- fits[,6]
siptw_ses <- fits[,7]
tmle_ses <- fits[,8]

diff_plugin <- (plugin_psis - Psi0)^2
diff_iptw <- (iptw_psis - Psi0)^2
diff_siptw <- (siptw_psis - Psi0)^2
diff_tmle <- (tmle_psis - Psi0)^2
plugin_CI95 <- wald_ci(plugin_psis, plugin_ses)
iptw_CI95 <- wald_ci(iptw_psis, iptw_ses)
siptw_CI95 <- wald_ci(siptw_psis, siptw_ses)
tmle_CI95 <- wald_ci(tmle_psis, tmle_ses)
coverage_plugin <- Psi0 >= plugin_CI95[,1] & Psi0 <= plugin_CI95[,2]
coverage_iptw <- Psi0 >= iptw_CI95[,1] & Psi0 <= iptw_CI95[,2]
coverage_siptw <- Psi0 >= siptw_CI95[,1] & Psi0 <= siptw_CI95[,2]
coverage_tmle <- Psi0 >= tmle_CI95[,1] & Psi0 <= tmle_CI95[,2]

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# estimator
dat_psis <- data.frame(Method=rep(c("Plug-In", "IPTW", "S-IPTW", "TML"), each=reps), 
                       Estimator=c(plugin_psis, iptw_psis, siptw_psis, tmle_psis))

mu <- ddply(dat_psis, "Method", summarise, grp.mean=mean(Estimator))
plt_hist <- ggplot(dat_psis, aes(x=Estimator, color=Method)) +
  geom_histogram(fill="white", position="dodge") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Method),
             linetype="dashed") +
  geom_vline(aes(xintercept=Psi0, color='Truth'))
theme(legend.position="top")
ggsave("point2_robust_estimator_histogram.png", 
       plot = plt_hist, path = "../plot")

# loss
dat_diff <- data.frame(Simulation_No.=rep(1:reps, 4), 
                       Method=rep(c("Plug-In", "IPTW", "S-IPTW", "TML"), each=reps), 
                       Loss_l2=c(diff_plugin, diff_iptw, diff_siptw, diff_tmle))

mse <- ddply(dat_diff, "Method", summarise, grp.mean=mean(Loss_l2))
plt_plot <- ggplot(data=dat_diff, aes(x=Simulation_No., y=Loss_l2, group=Method)) +
  geom_point(aes(color=Method)) +
  geom_hline(data=mse, aes(yintercept=grp.mean, color=Method),
             linetype="dashed")
ggsave("point2_robust_loss_point.png", 
       plot = plt_plot, path = "../plot")

# summary
sample_mean <- c(mean(plugin_psis), mean(iptw_psis), mean(siptw_psis), mean(tmle_psis))
sample_sd <- c(sd(plugin_psis), sd(iptw_psis), sd(siptw_psis), sd(tmle_psis))
coverage <- c(mean(coverage_plugin), mean(coverage_iptw), mean(coverage_siptw), mean(coverage_tmle))
dat_summary <- data.frame(Method=c("Plug-In", "IPTW", "S-IPTW", "TML"), 
                          Truth=Psi0, Sample_Mean=sample_mean, Sample_Sd=sample_sd, 
                          MSE=mse[match(mse$Method, c("Plug-In", "IPTW", "S-IPTW", "TML")),]$grp.mean,
                          Coverage_95=label_percent()(coverage))

g <- tableGrob(dat_summary, rows = NULL)
'
theme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.5)),
  colhead = list(fg_params=list(cex = 0.5)),
  rowhead = list(fg_params=list(cex = 0.5)))
g <- tableGrob(dat_summary, rows = NULL, theme = theme)
'
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))
grid.draw(g)

# efficiency
mean(tmle_ses * sqrt(obs) / se0)

