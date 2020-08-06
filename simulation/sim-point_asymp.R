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

# generate data
set.seed(1234)
n <- 1e5
W <- sample(1:4, n, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05))
S <- rbinom(n, 1, expit(1.4 - 0.6 * W))
Y <- rnorm(n, 2 + .5 * W, 1)

# truth
w = 1:4
pw = c(0.1, 0.2, 0.65, 0.05)
pw0 = (1 - expit(1.4 - 0.6 * w)) * pw
pw0 = pw0 / sum(pw0)
Ey0 = sum((2 + .5 * w) * pw0)

seed_fit = function(seed) {
  # generate data
  set.seed(seed)
  n <- 1000
  W <- sample(1:4, n, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05))
  S <- rbinom(n, 1, expit(1.4 - 0.6 * W))
  Y <- rnorm(n, 2 + .5 * W, 1)
  
  ## TMLE should not update in this case and
  ## should numerically match the non-parametric estimator
  
  ### non-parametric estimator ###
  psi_0 = function(b) mean(W == b & S == 0)
  psi_1 = function(b) mean(W == b & S == 1)
  psi_2 = function(b) mean((W == b & S == 1) * Y)
  psi_s = mean(S == 0)
  
  f_n = sum(sapply(levels(factor(W)),
                   function(b) psi_0(b)/psi_s * psi_2(b)/ifelse(psi_1(b), psi_1(b), 1)))
  D_n = rowSums(sapply(levels(factor(W)),
                       function(b) (psi_0(b) * (W == b & S == 1))/(psi_s * ifelse(psi_1(b), psi_1(b), 1)) * (Y - psi_2(b)/ifelse(psi_1(b), psi_1(b), 1)) +
                         (W == b & S == 0)/psi_s * (psi_2(b)/ifelse(psi_1(b), psi_1(b), 1) - f_n)))
  s_n = sd(D_n) / sqrt(n)
  CI95 <- wald_ci(f_n, s_n)
  
  return(Ey0 >= CI95[1] && Ey0 <= CI95[2])
}

reps <- 1e5
coverages <- c()
for (i in 1:reps) {
  coverages <- c(coverages, seed_fit(1234+i))
}

mean(coverages)
plt_plot <- ggplot(data=data.table(Simulation_No.=1:100, coverages=coverages[1:100]),
                   aes(x=Simulation_No., y=coverages)) +
  geom_point(size=0.1, aes(color=coverages))

