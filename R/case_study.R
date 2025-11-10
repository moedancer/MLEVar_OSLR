### Script to reproduce case study results
### Note: Due to copyright reasons, the digitised data points and KM-plots
###       from the original publication cannot be shared publicly.
###       Users need to digitise the KM-plots themselves to reproduce the results.
###       Hence, we provide purely simulated data.
###       Parameters from simulation are inspired by reconstructed data

library(IPDfromKM)
library(pracma)
library(flexsurv)

### Fix known details from control and experimental cohort
n_b <- 91
n_a <- 136
alloc_ratio <- n_b / n_a
events_b <- 68

### For both groups (entrol and experimental) we retrieve data in the following way:
# 1. Check if .csv files with reconstructed data already exist
# 2. If not, check if .Rda files with digitised points exist
# 3. If not, digitise points from .jpg files and save them as .Rda files

if (file_test("-f", "data/reconstructed_data_control.csv")) {
  control_ipd <- read.csv("data/reconstructed_data_control.csv")
} else if (file_test("-f", "control_points.Rda")) {
  load("control_points.Rda")
  control_pre <- preprocess(control_points, totalpts = n_a)
  control_ipd <- getIPD(control_pre, armID = 0)$IPD
  write.csv(control_ipd, "reconstructed_data_control.csv")
} else if (file_test("-f", "data/control_group.jpg")) {
  control_points <- getpoints(
    "control_group.jpg",
    x1 = 0,
    x2 = 8,
    y1 = 0,
    y2 = 100
  )
  save(control_points, file = "data/control_points.Rda")
  control_pre <- preprocess(control_points, totalpts = n_a)
  control_ipd <- getIPD(control_pre, armID = 0)$IPD
  write.csv(control_ipd, "reconstructed_data_control.csv")
} else {
  stop("No data found for control group.")
}

if (file_test("-f", "data/reconstructed_data_exp.csv")) {
  exp_ipd <- read.csv("data/reconstructed_data_exp.csv")
} else if (file_test("-f", "exp_points.Rda")) {
  load("exp_points.Rda")
  exp_pre <- preprocess(exp_points, totalpts = n_b)
  exp_ipd <- getIPD(exp_pre, tot.events = events_b, armID = 1)$IPD
  write.csv(exp_ipd, "reconstructed_data_exp.csv")
} else if (file_test("-f", "data/exp_group.jpg")) {
  exp_points <- getpoints(
    "exp_group.jpg",
    x1 = 0,
    x2 = 5,
    y1 = 0,
    y2 = 100
  )
  save(exp_points, file = "exp_points.Rda")
  exp_pre <- preprocess(exp_points, totalpts = n_b)
  exp_ipd <- getIPD(exp_pre, tot.events = events_b, armID = 1)$IPD
  write.csv(exp_ipd, "reconstructed_data_exp.csv")
} else {
  stop("No data found for experimental group.")
}

# Apply additional cut-off at 4.5 years
control_ipd$time <- pmin(control_ipd$time, 4.5)

all_data <- rbind(control_ipd, exp_ipd)
raw_comp <- survfit(Surv(time, status) ~ treat, data = all_data)
lr_test <- survdiff(Surv(time, status) ~ treat, data = all_data)

### Make MLE fit for control data
mle_fit <- flexsurvreg(
  Surv(time, status) ~ 1,
  data = control_ipd,
  dist = "llogis"
)
# Note: Estimator from flexsurv is on log-scale
mle_estimate <- exp(mle_fit$coefficients)
# Apply delta method to obtain variance estimate on the original scale
mle_var <- n_a * (diag(mle_estimate) %*% vcov(mle_fit) %*% diag(mle_estimate))

### Create plot for manuscript to show KM-curves and parametric fit
mle_surv_fct <- function(t) {
  pllogis(
    t,
    shape = mle_estimate[1],
    scale = mle_estimate[2],
    lower.tail = FALSE
  )
}
pdf("results/case_study/surv_curves_casestudy.pdf", width = 8, height = 6)
plot(
  raw_comp,
  col = c("blue", "orange"),
  lwd = 2,
  xlab = "Time",
  ylab = "Survival probability"
)
curve(mle_surv_fct, lty = 2, col = "blue", lwd = 2, add = TRUE)
legend(
  x = "topright",
  legend = c(
    "experimental group",
    "historic control",
    "historic control (parametric fit)"
  ),
  col = c("orange", "blue", "blue"),
  lty = c(1, 1, 2),
  lwd = 2
)
dev.off()

### Define hazard function and gradient w.r.t. parameter components
cum_hazard_fct <- function(t) {
  -log(pllogis(
    t,
    shape = mle_estimate[1],
    scale = mle_estimate[2],
    lower.tail = FALSE
  ))
}
gradient_cum_haz_fct <- function(t) {
  shape <- mle_estimate[1]
  scale <- mle_estimate[2]
  t_resc <- t / scale

  shape_derivative <- t_resc^shape * log(t_resc) / ((1 + t_resc)^shape)
  scale_derivative <- shape *
    t *
    t_resc^(shape - 1) /
    (scale^2 * (1 + t_resc)^shape)
  return(c(shape_derivative, scale_derivative))
}
gradient_cum_haz_fct <- Vectorize(gradient_cum_haz_fct, vectorize.args = "t")

### Compute unstandardised test statistic
test_statistic <- 1 /
  sqrt(n_b) *
  sum(exp_ipd$status - cum_hazard_fct(exp_ipd$time))

### Compute OSLR variance estimates
est_var_oslr_qv <- sum(exp_ipd$status) / n_b
est_var_oslr_pqv <- sum(cum_hazard_fct(exp_ipd$time)) / n_b
est_var_oslr_Wu <- (est_var_oslr_qv + est_var_oslr_pqv) / 2

### Compute additional variance from estimation uncertainty
mean_gradient_cum_haz <- apply(
  gradient_cum_haz_fct(exp_ipd$time),
  MARGIN = 1,
  FUN = mean
)
est_var_mle_error <- alloc_ratio *
  (mean_gradient_cum_haz %*% mle_var %*% mean_gradient_cum_haz)

### Compute test statistics and p-values
oslr_std_test_pqv <- test_statistic / sqrt(est_var_oslr_pqv)
pnorm(oslr_std_test_pqv)
oslr_std_test_Wu <- test_statistic / sqrt(est_var_oslr_pqv)
pnorm(oslr_std_test_Wu)
std_test_pqv <- test_statistic / sqrt(est_var_oslr_pqv + est_var_mle_error)
pnorm(std_test_pqv)
std_test_Wu <- test_statistic / sqrt(est_var_oslr_Wu + est_var_mle_error)
pnorm(std_test_Wu)
