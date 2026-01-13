### Script to simulate type I error results for one-sample log-rank test
### with reference curves estimated via MLE under Weibull distribution
###
### Script can be run from command line with 3 arguments:
### 1st argument: n_a (sample size of control group)
### 2nd argument: n_b (sample size of experimental group)
### 3rd argument: shape parameter of Weibull distribution
### (Values default to n_a = 100, n_b = 100, shape = 1 if not provided)
###
### Alternatively, values can also be set directly in the script

passed_args <- commandArgs(trailingOnly = TRUE)
n_a <- 100
n_b <- 100
my_shape <- 1
hr <- 1

if (length(passed_args) >= 4) {
  hr <- as.numeric(passed_args[[4]])
}
if (length(passed_args) >= 3) {
  my_shape <- as.numeric(passed_args[[3]])
}
if (length(passed_args) >= 2) {
  n_b <- as.integer(passed_args[[2]])
}
if (length(passed_args) >= 1) {
  n_a <- as.integer(passed_args[[1]])
}

require(flexsurv)
require(nph)

alloc_ratio <- n_b / n_a

# Define remaining parameters of survival time distribution
one_year_surv <- 0.5
my_scale <- (-log(one_year_surv))^(-1 / my_shape)
# Scale parameter in treatment group
my_scale_tr <- my_scale / (hr^(1 / my_shape))

cum_hazard_fct <- function(time, shape, scale) (time / scale)^shape
gradient_cum_haz_fct <- function(time, shape, scale) {
  shape_derivative <- log(time / scale) * cum_hazard_fct(time, shape, scale)
  scale_derivative <- (-1) *
    shape *
    (time / scale)^(shape - 1) *
    (time / scale^2)
  return(c(shape_derivative, scale_derivative))
}
gradient_cum_haz_fct <- Vectorize(gradient_cum_haz_fct, vectorize.args = "time")

accrual <- 2
follow_up <- 3

simulation_runs <- 100000
weibull_mles <- matrix(0, ncol = 2, nrow = simulation_runs)

hazard_mle <- rep(0, simulation_runs)
est_hazard_random_error <- rep(0, simulation_runs)
oslr_test_statistic <- rep(0, simulation_runs)

# Save p-value of log-rank test
lr_p <- rep(0, simulation_runs)

# estimate variances
est_var_oslr_qv <- rep(0, simulation_runs)
est_var_oslr_pqv <- rep(0, simulation_runs)
est_var_oslr_pqv_mle <- rep(0, simulation_runs)
est_var_mle_error <- rep(0, simulation_runs)

for (i in 1:simulation_runs) {
  raw_time_control <- rweibull(n_a, shape = my_shape, scale = my_scale)
  censoring_time_control <- runif(
    n_a,
    min = follow_up,
    max = accrual + follow_up
  )

  time_control <- pmin(raw_time_control, censoring_time_control)
  event_control <- raw_time_control <= censoring_time_control

  control_data <- data.frame(time = time_control, event = event_control)

  ### Estimation via flexsurv
  mle_fit <- flexsurvreg(
    Surv(time_control, event_control) ~ 1,
    data = control_data,
    dist = "weibull"
  )
  # Note: Estimator from flexsurv is on log-scale
  mle_estimate <- exp(mle_fit$coefficients)
  # Apply delta method to obtain variance estimate on the original scale
  mle_var <- n_a * (diag(mle_estimate) %*% vcov(mle_fit) %*% diag(mle_estimate))

  weibull_mles[i, ] <- mle_estimate

  raw_time_exp <- rweibull(n_b, shape = my_shape, scale = my_scale_tr)
  censoring_time_exp <- runif(n_b, min = follow_up, max = accrual + follow_up)

  time_exp <- pmin(raw_time_exp, censoring_time_exp)
  event_exp <- raw_time_exp <= censoring_time_exp

  exp_data <- data.frame(time = time_exp, event = event_exp)

  all_data <- rbind(control_data, exp_data)
  all_data$group <- c(rep(0, n_a), rep(1, n_b))
  lr_test <- logrank.test(
    time = all_data$time,
    event = all_data$event,
    group = all_data$group,
    alternative = "less"
  )
  lr_p[i] <- lr_test$test$p

  est_hazard_random_error[i] <- 1 /
    sqrt(n_b) *
    sum(
      cum_hazard_fct(time_exp, my_shape, my_scale) -
        cum_hazard_fct(time_exp, mle_estimate[1], mle_estimate[2])
    )

  oslr_test_statistic[i] <- 1 /
    sqrt(n_b) *
    sum(event_exp - cum_hazard_fct(time_exp, my_shape, my_scale))

  est_var_oslr_qv[i] <- sum(event_exp) / n_b
  est_var_oslr_pqv[i] <- sum(cum_hazard_fct(time_exp, my_shape, my_scale)) / n_b
  est_var_oslr_pqv_mle[i] <- sum(cum_hazard_fct(
    time_exp,
    mle_estimate[1],
    mle_estimate[2]
  )) /
    n_b

  mean_gradient_cum_haz <- apply(
    gradient_cum_haz_fct(time_exp, mle_estimate[1], mle_estimate[2]),
    MARGIN = 1,
    FUN = mean
  )
  est_var_mle_error[i] <- alloc_ratio *
    (mean_gradient_cum_haz %*% mle_var %*% mean_gradient_cum_haz)
}

results <- data.frame(
  oslr = oslr_test_statistic,
  hazard_diff = est_hazard_random_error,
  est_var_oslr_qv = est_var_oslr_qv,
  est_var_oslr_pqv = est_var_oslr_pqv,
  est_var_oslr_pqv_mle = est_var_oslr_pqv_mle,
  est_var_mle_error = est_var_mle_error,
  lr_p = lr_p
)

results$n_b <- n_b
results$alloc_ratio <- alloc_ratio
results$shape <- my_shape
results$hr <- hr

save(
  results,
  file = paste(
    "results/single_scenarios/t2e_raw_KAPPA",
    sub(x = my_shape, pattern = "\\.", replacement = "dec"),
    "_NA",
    n_a,
    "_NB",
    n_b,
    "_HR",
    sub(x = round(hr, digits = 4), pattern = "\\.", replacement = "dec"),
    ".Rda",
    sep = ""
  )
)
