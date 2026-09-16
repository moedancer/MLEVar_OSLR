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

hpc_version <- FALSE # Is this the version running on PALMA II?
if (hpc_version) {
  .libPaths("/home/d/danzerm/R/library")
}

require(survival)
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

cum_hazard_fct_exp <- function(time, rate) time * rate
gradient_cum_haz_fct_exp <- function(time, rate) time
gradient_cum_haz_fct_exp <- Vectorize(
  gradient_cum_haz_fct_exp,
  vectorize.args = "time"
)

accrual <- 2
follow_up <- 3

simulation_runs <- 100000
weibull_mles <- matrix(0, ncol = 2, nrow = simulation_runs)
if (my_shape == 1) {
  exp_mles <- rep(0, simulation_runs)
}

est_hazard_random_error <- rep(0, simulation_runs)
if (my_shape == 1) {
  est_hazard_random_error_exp <- rep(0, simulation_runs)
}
est_hazard_random_error_np <- rep(0, simulation_runs)
oslr_test_statistic <- rep(0, simulation_runs)

# Save p-value of log-rank test
lr_p <- rep(0, simulation_runs)

# estimate variances
est_var_oslr_qv <- rep(0, simulation_runs)
est_var_oslr_pqv <- rep(0, simulation_runs)
est_var_oslr_pqv_mle <- rep(0, simulation_runs)
if (my_shape == 1) {
  est_var_oslr_pqv_mle_exp <- rep(0, simulation_runs)
}
est_var_mle_error <- rep(0, simulation_runs)
if (my_shape == 1) {
  est_var_mle_error_exp <- rep(0, simulation_runs)
}

# Addition: Estimate variances for non-parametric reference curve
est_var_oslr_pqv_na <- rep(0, simulation_runs)
est_var_na_error <- rep(0, simulation_runs)

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

  # Repeat the procedure with estimation by exponential distribution if shape parameter is 1
  if (my_shape == 1) {
    mle_fit_exp <- flexsurvreg(
      Surv(time_control, event_control) ~ 1,
      data = control_data,
      dist = "exp"
    )
    # Note: Estimator from flexsurv is on log-scale
    mle_estimate_exp <- exp(mle_fit_exp$coefficients)
    # Apply delta method to obtain variance estimate on the original scale
    mle_var_exp <- n_a * mle_estimate_exp^2 * vcov(mle_fit_exp)

    exp_mles[i] <- mle_estimate_exp
  }

  # Addition: Also use non-parametric procedure to fit reference curve
  np_fit <- survfit(Surv(time_control, event_control) ~ 1, data = control_data)

  # Generate data for experimental group
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

  # Compute difference between real hazards and estimated hazards from reference curve...
  # ...in parametric setting
  est_hazard_random_error[i] <- 1 /
    sqrt(n_b) *
    sum(
      cum_hazard_fct(time_exp, my_shape, my_scale) -
        cum_hazard_fct(time_exp, mle_estimate[1], mle_estimate[2])
    )
  if (my_shape == 1) {
    est_hazard_random_error_exp[i] <- 1 /
      sqrt(n_b) *
      sum(
        cum_hazard_fct(time_exp, my_shape, my_scale) -
          cum_hazard_fct_exp(time_exp, mle_estimate_exp) # Note the different parametrizations of exponential and Weibull distribution
      )
  }
  # ...in non-parametric setting
  reference_evaluation_na <- summary(
    np_fit,
    times = exp_data$time,
    extend = TRUE
  )
  reference_na_values <- reference_evaluation_na$cumhaz
  est_hazard_random_error_np[i] <- 1 /
    sqrt(n_b) *
    sum(
      cum_hazard_fct(time_exp, my_shape, my_scale) -
        reference_na_values
    )

  # Compute one-sample log-rank statistics based on true curve
  oslr_test_statistic[i] <- 1 /
    sqrt(n_b) *
    sum(event_exp - cum_hazard_fct(time_exp, my_shape, my_scale))

  # Compute naive variance estimates
  # ...based on true cumulative hazard function
  est_var_oslr_qv[i] <- sum(event_exp) / n_b
  est_var_oslr_pqv[i] <- sum(cum_hazard_fct(time_exp, my_shape, my_scale)) / n_b
  # ...based on estimated parametric hazard function
  est_var_oslr_pqv_mle[i] <- sum(cum_hazard_fct(
    time_exp,
    mle_estimate[1],
    mle_estimate[2]
  )) /
    n_b
  if (my_shape == 1) {
    est_var_oslr_pqv_mle_exp[i] <- sum(cum_hazard_fct_exp(
      time_exp,
      mle_estimate_exp
    )) /
      n_b
  }
  # ...based on estimated non-parametric hazard function
  est_var_oslr_pqv_na <- sum(reference_na_values)

  # Estimate additional variance...
  # ...for parametric reference curve
  mean_gradient_cum_haz <- apply(
    gradient_cum_haz_fct(time_exp, mle_estimate[1], mle_estimate[2]),
    MARGIN = 1,
    FUN = mean
  )
  est_var_mle_error[i] <- alloc_ratio *
    (mean_gradient_cum_haz %*% mle_var %*% mean_gradient_cum_haz)

  if (my_shape == 1) {
    mean_gradient_cum_haz_exp <- mean(gradient_cum_haz_fct_exp(
      time_exp,
      mle_estimate_exp
    ))
    est_var_mle_error_exp[i] <- alloc_ratio *
      (mean_gradient_cum_haz_exp %*% mle_var_exp %*% mean_gradient_cum_haz_exp)
  }
  # ...for non-parametric reference curve
  pairwise_min_times <- rep(
    sort(exp_data$time),
    times = seq(from = 2 * n_b - 1, to = 1, by = -2)
  )
  reference_evaluation_gw <- summary(
    np_fit,
    times = pairwise_min_times,
    extend = TRUE
  )
  reference_gw_values <- reference_evaluation_gw$std.chaz^2
  est_var_na_error[i] <- sum(reference_gw_values) / n_b
}

if (my_shape != 1) {
  results <- data.frame(
    oslr = oslr_test_statistic,
    hazard_diff = est_hazard_random_error,
    hazard_diff_np = est_hazard_random_error_np,
    est_var_oslr_qv = est_var_oslr_qv,
    est_var_oslr_pqv = est_var_oslr_pqv,
    est_var_oslr_pqv_mle = est_var_oslr_pqv_mle,
    est_var_oslr_pqv_na = est_var_oslr_pqv_na,
    est_var_mle_error = est_var_mle_error,
    est_var_na_error = est_var_na_error,
    lr_p = lr_p
  )
} else {
  results <- data.frame(
    oslr = oslr_test_statistic,
    hazard_diff = est_hazard_random_error,
    hazard_diff_exp = est_hazard_random_error_exp,
    hazard_diff_np = est_hazard_random_error_np,
    est_var_oslr_qv = est_var_oslr_qv,
    est_var_oslr_pqv = est_var_oslr_pqv,
    est_var_oslr_pqv_mle = est_var_oslr_pqv_mle,
    est_var_oslr_pqv_mle_exp = est_var_oslr_pqv_mle_exp,
    est_var_oslr_pqv_na = est_var_oslr_pqv_na,
    est_var_mle_error = est_var_mle_error,
    est_var_mle_error_exp = est_var_mle_error_exp,
    est_var_na_error = est_var_na_error,
    lr_p = lr_p
  )
}

results$n_b <- n_b
results$alloc_ratio <- alloc_ratio
results$shape <- my_shape
results$hr <- hr

save(
  results,
  file = paste(
    "results/single_scenarios/t2e_raw_w_np_KAPPA",
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
