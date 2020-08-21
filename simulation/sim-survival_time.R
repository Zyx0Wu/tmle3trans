library(simcausal)
library(testthat)
library(data.table)
library(ggplot2)
library(pammtools)

library(sl3)
library(tmle3)
library(tmle3trans)

simulate_data <- function(n_sim = 2e2) {
  D <- DAG.empty()
  D <- D +
    node("S", distr = "rbinom", size = 1, prob = .6) +
    node("W1", distr = "rbinom", size = 1, prob = .5) +
    node("W", distr = "runif", min = .1 - .1 * S, max = 1.6 - .1 * S) +
    node("Trexp", distr = "rexp", rate = .2 + .7 * W^2) +
    node("Cweib", distr = "rweibull", shape = 1 + .5 * W, scale = 75) +
    node("T", distr = "rconst", const = round(Trexp * 2) + 1) +
    node("C", distr = "rconst", const = round(Cweib * 2) + 1) +
    # Observed random variable (follow-up time):
    node("T.tilde", distr = "rconst", const = ifelse(T <= C, T, C)) +
    # Observed random variable (censoring indicator, 1 - failure event, 0 - censored):
    node("Delta", distr = "rconst", const = ifelse(T <= C, 1, 0))
  setD <- set.DAG(D)
  dat <- sim(setD, n = n_sim)
  # only grab W's, S, T.tilde, Delta
  Wname <- grep("W", colnames(dat), value = TRUE)
  dat <- dat[, c(Wname, "S", "T.tilde", "Delta")]
  return(list(dat = dat))
}

set.seed(4321)
n_sim <- 2e4
simulated <- simulate_data(n_sim = n_sim)
data <- simulated$dat

# TODO: check
while (all(data$Delta == 1)) {
  simulated <- simulate_data(n_sim = n_sim)
  data <- simulated$dat
}

tmax <- max(data$T.tilde)
times <- 1:tmax

node_list <- list(W = c("W", "W1"), S = "S", T.tilde = "T.tilde", Delta = "Delta")

transformed <- transform_data(data, node_list)
long_data <- transformed$long_data
long_node_list <- transformed$long_node_list

# observed survival curve
max_time <- 20
sF <- 1
for (i in times) {
  sF <- c(sF, sF[length(sF)] - mean(long_data[((i-1)*n_sim+1):(i*n_sim),][data$S==0,]$Failed))
}
sF <- sF[-1][1:max_time]

seed_fit = function(seed, n_sim) {
  set.seed(seed)
  simulated <- simulate_data(n_sim = n_sim)
  data <- simulated$dat
  
  # TODO: check
  while (all(data$Delta == 1) | max(data$T.tilde) < max_time) {
    simulated <- simulate_data(n_sim = n_sim)
    data <- simulated$dat
  }
  
  tmax <- max(data$T.tilde)
  times <- 1:tmax
  
  node_list <- list(W = c("W", "W1"), S = "S", T.tilde = "T.tilde", Delta = "Delta")
  
  transformed <- transform_data(data, node_list)
  long_data <- transformed$long_data
  long_node_list <- transformed$long_node_list
  
  ### sl3 ###
  lrnr_mean <- make_learner(Lrnr_mean)
  lrnr_glm <- make_learner(Lrnr_glm)
  lrnr_sl <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm))
  learner_list <- list(S = lrnr_sl, failed = lrnr_sl, censored = lrnr_sl)
  
  ### tmle3 ###
  # TODO: check
  # only target time point 1
  tmle_spec <- tmle_SOT(1, 0, target_times = times[1:min(5, tmax)])
  tmle_task <- tmle_spec$make_tmle_task(long_data, long_node_list)
  
  initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)
  '
  updater <- tmle3_Update$new(maxit=5, verbose = TRUE)
  
  updater <- tmle3_Update$new(
    constrain_step = TRUE, one_dimensional = TRUE, 
    delta_epsilon = 3e-2, verbose = TRUE, 
    convergence_type = "scaled_var", maxit = 10
  )
  '
  updater <- tmle3_Update_survival$new(
    maxit = 1e2, 
    cvtmle = TRUE,
    convergence_type = "scaled_var",
    delta_epsilon = 1e-2,
    fit_method = "l2",
    use_best = TRUE,
    verbose = TRUE
  )
  
  targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = updater)
  tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
  
  tmle_fit_manual <- fit_tmle3(
    tmle_task, targeted_likelihood, tmle_params,
    targeted_likelihood$updater
  )
  
  # extract results
  tmle_summary <- tmle_fit_manual$summary
  init_psi <- tmle_summary$init_est[1:max_time]
  tmle_psi <- tmle_summary$tmle_est[1:max_time]
  
  # loss
  l2_diff_sl <- (init_psi - sF)^2
  l2_diff_tl <- (tmle_psi - sF)^2
  
  return(list(init_psi, tmle_psi, l2_diff_sl, l2_diff_tl))
}

reps <- 50
obs <- 2e2
sl_psis <- c()
tl_psis <- c()
diff_sl <- c()
diff_tl <- c()
for (i in 1:reps) {
  fits <- seed_fit(1234+i, obs)
  sl_psis <- cbind(sl_psis, fits[[1]])
  tl_psis <- cbind(tl_psis, fits[[2]])
  diff_sl <- cbind(diff_sl, fits[[3]])
  diff_tl <- cbind(diff_tl, fits[[4]])
}
'
psis_df <- data.frame(time=0:max_time,
                      sl=c(1,apply(sl_psis, 1, mean)),
                      sl_upper=c(1,apply(sl_psis, 1, mean)+apply(sl_psis, 1, sd)),
                      sl_lower=c(1,apply(sl_psis, 1, mean)-apply(sl_psis, 1, sd)),
                      tl_10=c(1,apply(tl_psis_10, 1, mean)),
                      tl_10_upper=c(1,apply(tl_psis_10, 1, mean)+apply(tl_psis_10, 1, sd)),
                      tl_10_lower=c(1,apply(tl_psis_10, 1, mean)-apply(tl_psis_10, 1, sd)),
                      tl_30=c(1,apply(tl_psis_30, 1, mean)),
                      tl_30_upper=c(1,apply(tl_psis_30, 1, mean)+apply(tl_psis_30, 1, sd)),
                      tl_30_lower=c(1,apply(tl_psis_30, 1, mean)-apply(tl_psis_30, 1, sd)),
                      tl_100=c(1,apply(tl_psis_100, 1, mean)),
                      tl_100_upper=c(1,apply(tl_psis_100, 1, mean)+apply(tl_psis_100, 1, sd)),
                      tl_100_lower=c(1,apply(tl_psis_100, 1, mean)-apply(tl_psis_100, 1, sd)),
                      truth=c(1,sF))

plt_step <- ggplot(data=psis_df, aes(time)) + 
  #geom_stepribbon(aes(ymin=sl_lower, ymax=sl_upper), fill = "grey") + 
  #geom_stepribbon(aes(ymin=tl_100_lower, ymax=tl_100_upper), fill = "grey") + 
  geom_step(aes(y=sl), color="red", alpha=0.8) + 
  geom_step(aes(y=tl_100), color="orange", alpha=0.8) + 
  geom_step(aes(y=truth), color="blue", alpha=0.8) + 
  ylab("survival")

plt_line <- ggplot(data=data.frame(step=1:99, ic=ic_step)) +
  geom_line(aes(x=step, y=ic)) + 
  ylab("validation ic")

plt_line <- ggplot(data=data.frame(time=1:5, 
                                   mse_10=apply((tl_psis_10-sF)^2, 1, mean)[1:5], 
                                   mse_30=apply((tl_psis_30-sF)^2, 1, mean)[1:5], 
                                   mse_100=apply((tl_psis_100-sF)^2, 1, mean)[1:5])) +
  geom_line(aes(x=time, y=mse_10), color="green") + 
  geom_line(aes(x=time, y=mse_30), color="brown") + 
  geom_line(aes(x=time, y=mse_100), color="purple") + 
  ylab("mse")
'
  