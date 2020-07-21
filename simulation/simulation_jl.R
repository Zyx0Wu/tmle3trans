context("Test - Robustness")

library(sl3)
library(tmle3)
library(tmle3trans)
library(uuid)
library(assertthat)
library(data.table)
library(future)
library(dplyr)
library(foreach)
library(parallel)
library(ggplot2)

# setup data for test
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
pop = gendata_trans(1e5)
pop
popS0 =  filter(pop, S==0)
popS1 =  filter(pop, S==1)
Psi0 = mean(popS0$Y)
mean(popS1$Y)
Ypop = pop$Y
Spop = pop$S

### fit ###
popW = select(pop, node_list$W)
node_list <- list(W = c("W1", "W2", "W3"), S = "S", Y = "Y")
W <- select(pop, node_list$W)

### 1. Naive ###
# note: not Plug-In
bias = T
Qpop = mutate(W, Q = -1 - .5 * W1 +
                .8 * W2 + .2 * W3 + bias * n^(-0.3))$Q
Qpop
PSpop = mutate(W, g=prob_clip(expit(1.4 - 0.6 * W1 - 2 * W2 +
                                    0.7 * W3 + bias * n^(-0.3))))$g
PSpop
Psi0

Dpop = Spop*(1-PSpop)/PSpop/mean(1-Spop)*(Ypop - Qpop)+
  (1-Spop)/mean(1-Spop)*(Qpop-mean(Qpop[Spop==0]))

se0 = sd(Dpop)

###
###
# This is your simulation function except there is not a method
# for giving CI for naive_psi--use delta method as you have coded
###
###

simJL = function(n, biasS, biasY, wt=TRUE) {
  # n=100
  # biasS=biasY=T
  # wt=T
  data = gendata_trans(n)
  node_list <- list(W = c("W1", "W2", "W3"), S = "S", Y = "Y")

  W <- select(data, node_list$W)

  s_dens <- function(w, bias=FALSE, n) prob_clip(expit(1.4 - 0.6 * w[1] - 2 * w[2] + 0.7 * w[3] + bias * n^(-0.3)))
  y_dens <- function(w, bias=FALSE, n) -1 - .5 * w[1] +
    .8 * w[2] + .2 * w[3] + bias * n^(-0.3)

  ### 1. Naive ###
  # note: not Plug-In
  Q = apply(W, 1, y_dens, bias = biasY, n = n)

  S = data$S
  Y = data$Y

  ### 2. IPTW ###
  pS1W <- apply(W, 1, s_dens,n=n, bias=biasS)

  pS0W <- 1 - pS1W
  #pS0 <- mean(data$S==0)
  pS0 <- 1 - mean(pS1W)

  H1 = S*pS0W/(pS0*pS1W)

  iptw_psi_init <- mean(H1*Y)
  iptw_psi = iptw_psi_init/(mean(H1))
  iptw_se <- sd(H1*(Y-iptw_psi_init))/sqrt(n)
  iptw_CI95 <- c(iptw_psi, iptw_psi - 1.96*iptw_se, iptw_psi + 1.96*iptw_se)
  iptw_CI95

  ### 3. TML ###
  
  #tmle_psi_init = mean(Q[S==0])
  tmle_psi_init = mean((S == 0)/pS0 * Q)
  
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

  #tmlejl_psi = mean(Qstar[S==0])
  tmlejl_psi = mean((S == 0)/pS0 * Qstar)
  tmle_psi_init
  tmlejl_psi
  tmlejl_se <- sd(D)/sqrt(n)

  tmlejl_CI95 <- c(tmlejl_psi, tmlejl_psi - 1.96*tmlejl_se, tmlejl_psi + 1.96*tmlejl_se)
  tmlejl_CI95
  
  ### 3. TML3 ###
  
  tmle_spec <- tmle_AOT(1, 0, fit_s_marginal = "integral")
  
  # define data
  tmle_task <- tmle_spec$make_tmle_task(data, node_list)
  
  # define likelihood
  g_dens <- function(task) apply(task$get_node("covariates"), 1, s_dens, n=n, bias=TRUE)
  Q_mean <- function(task) apply(task$get_node("covariates"), 1, y_dens, n=n, bias=TRUE)
  
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
  browser()
  

  return(list(tmle =  tmle_CI95,
              iptw = iptw_CI95,
              naive_psi_init = tmle_psi_init))
}

simJL(1000,T,T,F)
cl = makeCluster(8)

n=1000
B=5000
res = foreach(i=1:B, .errorhandling = "remove") %do%
  simJL(n,T,T,T)

getres = function(res) {

  res = do.call(rbind,lapply(res, function(L) unlist(L)))
  res = data.frame(res)
  head(res)
  tmlecov = mean(res$tmle2 <= Psi0 & res$tmle3 >= Psi0)
  tmleeff = mean((res$tmle3- res$tmle2)/(2*1.96)/se0_n)
  iptwcov = mean(res$iptw2 <= Psi0 & res$iptw3 >= Psi0)
  iptweff = mean((res$iptw3- res$iptw2)/(2*1.96)/se0_n)

  coverage = c(tmle = tmlecov,iptw = iptwcov)
  efficiency = c(tmle = tmleeff, iptw = iptweff)
  cov_eff = rbind(coverage, efficiency)
  cov_eff

  mse = apply(res[,c(1,4,7)], 2, function(x) mean((x-Psi0)^2))
  mse

  ests = unlist(res[,c(1,4,7)])
  names(ests) = NULL
  ests
  L = length(ests)/3
  type = c(rep("tmle",L),
           rep("iptw",L), rep("naive_psi_init",L))
  type
  plotdf = data.frame(ests = ests, type = type)
  return(list(mse = mse, cov_eff = cov_eff, plotdf = plotdf))
}

se0_n = se0/sqrt(n)
Psi0
results = getres(res)

results[1:2]

sampling_dists = ggplot(results$plotdf, aes(x=ests, fill = type))+
  geom_density(alpha=.4)+geom_vline(xintercept=Psi0)+
  xlim(-1.3,-.4)

sampling_dists

