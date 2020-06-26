context("Test - Effectiveness")

library(sl3)
library(tmle3)
library(tmle3tr)
library(uuid)
library(assertthat)
library(data.table)
library(future)

# note: S <=> A, site variable

# setup data for test
washb_data <- fread("https://raw.githubusercontent.com/tlverse/tlverse-data/master/wash-benefits/washb_data.csv", stringsAsFactors = TRUE)
node_list <- list(
  W = c(
    "momheight", "aged", "momage", "Nlt18", "Ncomp", "watmin",
    "elec", "floor", "walls", "roof", "asset_wardrobe",
    "asset_table", "asset_chair", "asset_khat",
    "asset_chouki", "asset_refrig", "asset_bike",
    "asset_moto", "asset_sewmach", "asset_mobile"
  ),
  A = "asset_tv", # pretested insignificant variable to the outcome
  Y = "whz"
)
processed <- process_missing(washb_data, node_list)
data <- processed$data
node_list <- processed$node_list

### 1. naive ###
W <- data[ ,colnames(data) %in% node_list$W, with=FALSE]
Y <- data[ ,colnames(data) %in% node_list$Y, with=FALSE]

W1 <- W[data[[node_list$A]] == 1, ]
Y1 <- Y[data[[node_list$A]] == 1, ]
W0 <- W[data[[node_list$A]] == 0, ]
Y0 <- Y[data[[node_list$A]] == 0, ]

fit <- glm(paste(node_list$Y, "~", paste(node_list$W, collapse = " + ")),
           data = cbind(W1, Y1))
est <- predict(fit, newdata = W0, type = 'response', se.fit = TRUE)
psis <- predict(fit, newdata = W0, type = 'response', se.fit = TRUE)$fit
stds <- predict(fit, newdata = W0, type = 'response', se.fit = TRUE)$se.fit

psi <- mean(psis)
std <- sqrt(mean(stds^2)) / sqrt(length(stds))
CI95 <- sprintf("(%f, %f)", psi - 1.96*std, psi + 1.96*std)

### 2. TMLE ###
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
mn_metalearner <- make_learner(
  Lrnr_solnp, metalearner_linear_multinomial,
  loss_loglik_multinomial
)
Q_learner <- make_learner(Lrnr_sl, qlib, ls_metalearner)
g_learner <- make_learner(Lrnr_sl, glib, mn_metalearner)
learner_list <- list(Y = Q_learner, A = g_learner)

tmle_spec <- tmle_TR(1, 0)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
factor_list <- list(
  define_lf(LF_emp, "W"),
  define_lf(LF_fit, "A", learner = learner_list[["A"]]),
  define_lf(LF_fit_site, "Y", learner = learner_list[["Y"]], type = "mean")
)

initial_likelihood <- Likelihood$new(factor_list)$train(tmle_task)

# define update method (submodel + loss function)
# disable cvtmle for this test to compare with tmle package
updater <- tmle3_Update$new(cvtmle = FALSE, convergence_type = "sample_size")
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

tmle_CI95 <- sprintf("(%f, %f)", tmle_psi - 1.96*tmle_se, tmle_psi + 1.96*tmle_se)
