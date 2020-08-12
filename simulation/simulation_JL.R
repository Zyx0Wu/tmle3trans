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

### 1. Naive ###
# note: not Plug-In
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
  
  data0 <- filter(data, S==0)
  data1 <- filter(data, S==1)
  
  ### fit ###
  WS1 = select(data1, node_list$W)
  YS1 <- data1$Y
  WS0 <- select(data0, node_list$W)
  
  s_dens <- function(w, bias=FALSE, n)
    prob_clip(expit(1.4 - 0.6 * w[1] - 2 * w[2] + 0.7 * w[3] +
                      bias * n^(-0.3)))
  y_dens <- function(w, bias=FALSE, n) -1 - .5 * w[1] +
    .8 * w[2] + .2 * w[3] + bias * n^(-0.3)
  
  ### 1. Naive ###
  # note: not Plug-In
  Q = apply(W, 1, y_dens, bias = biasY, n = n)
  
  S = data$S
  Y = data$Y
  
  naive_psi <- mean(Q[S==0])
  naive_psi
  
  ###
  ###
  # use the delta method for inference here for naive_psi
  ###
  ###
  
  ### 2. IPTW ###
  pS1W <- apply(W, 1, s_dens,n=n, bias=T)
  
  pS0W <- 1 - pS1W
  pS0 <- mean(data$S==0)
  
  H1 = S*pS0W/(pS0*pS1W)
  
  iptw_psi_init <- mean(H1*Y)
  iptw_psi = iptw_psi_init/(mean(H1))
  iptw_se <- sd(H1*(Y-iptw_psi_init))/sqrt(n)
  iptw_CI95 <- c(iptw_psi, iptw_psi - 1.96*iptw_se, iptw_psi + 1.96*iptw_se)
  iptw_CI95
  
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
  
  return(list(tmle =  tmle_CI95,
              iptw = iptw_CI95,
              tmle_psi_init = tmle_psi_init,
              naive_psi = naive_psi))
}

simJL(1000,T,T,T)
cl = makeCluster(8)

N=1000
B=5000
res = foreach(i=1:B, .errorhandling = "remove") %do%
  simJL(N,T,T,T)

getres = function(res, se0,boots, N) {
  
  se0_n = se0/sqrt(N)
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
  
  mse = apply(res[,c(1,4,7,8)], 2, function(x) mean((x-Psi0)^2))
  mse
  
  ests = unlist(res[,c(1,4,7,8)])
  names(ests) = NULL
  ests
  type = c(rep("tmle",boots),
           rep("iptw",boots), rep("tmle_init",boots),
           rep("naive_psi", boots))
  type
  plotdf = data.frame(ests = ests, type = type)
  return(list(mse = mse, cov_eff = cov_eff, plotdf = plotdf))
}


Psi0
results = getres(res, se0,B, N)

# efficiency and coverage
results[1:2]
plotdf = results$plotdf

sampling_dists = ggplot(plotdf, aes(x=ests, fill = type))+
  geom_density(alpha=.3)+geom_vline(xintercept=Psi0)+
  xlim(-1.2,-.1)

sampling_dists

