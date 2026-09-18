library(tidyverse)
library(patchwork)

# Load and combine results from all single scenario files under null hypothesis
results_all <- NULL
for (tmp_file in list.files(
  path = "results/single_scenarios/w_np/"
)) {
  load(paste("results/single_scenarios/w_np/", tmp_file, sep = ""))
  results_all <- bind_rows(results_all, results)
}

# Compute missing sample size
results_all$n_a <- results_all$n_b / results_all$alloc_ratio
# Set hazard ratio to 1 for simulations under H0
results_all$hr[is.na(results_all$hr)] <- 1

### Compute different variance estimates
## Compute Wu's variance with true parameters
results_all$est_var_oslr_Wu <- 0.5 *
  results_all$est_var_oslr_qv +
  0.5 * results_all$est_var_oslr_pqv
## Compute Wu's variance with estimated parameters
# Version 1: Use Weibull estimate
results_all$est_var_oslr_Wu_mle <- 0.5 *
  results_all$est_var_oslr_qv +
  0.5 * results_all$est_var_oslr_pqv_mle
# Version 2: Use exponential estimate
results_all$est_var_oslr_Wu_mle_exp <- 0.5 *
  results_all$est_var_oslr_qv +
  0.5 * results_all$est_var_oslr_pqv_mle_exp
# Version 3: Use non-parametric estimate
results_all$est_var_oslr_Wu_na <- 0.5 *
  results_all$est_var_oslr_qv +
  0.5 * results_all$est_var_oslr_pqv_na

## Compute total variance (including reference curve uncertainty)
# Version 1: Use Weibull estimates
results_all$est_var_total_qv <- results_all$est_var_oslr_qv +
  results_all$est_var_mle_error
results_all$est_var_total_pqv <- results_all$est_var_oslr_pqv_mle +
  results_all$est_var_mle_error
results_all$est_var_total_Wu <- results_all$est_var_oslr_Wu_mle +
  results_all$est_var_mle_error
# Version 2: Use exponential estimates
results_all$est_var_total_qv_exp <- results_all$est_var_oslr_qv +
  results_all$est_var_mle_error_exp
results_all$est_var_total_pqv_exp <- results_all$est_var_oslr_pqv_mle_exp +
  results_all$est_var_mle_error_exp
results_all$est_var_total_Wu_exp <- results_all$est_var_oslr_Wu_mle_exp +
  results_all$est_var_mle_error_exp
# Version 3: Use non-parametric estimates
results_all$est_var_total_qv_na <- results_all$est_var_oslr_qv +
  results_all$est_var_na_error
results_all$est_var_total_pqv_na <- results_all$est_var_oslr_pqv_na +
  results_all$est_var_na_error
results_all$est_var_total_Wu_na <- results_all$est_var_oslr_Wu_na +
  results_all$est_var_na_error

# Compute test statistics and p-values for true (but in practice unknown) reference hazard
results_all$oslr_std_qv <- results_all$oslr / sqrt(results_all$est_var_oslr_qv)
results_all$pos_oslr_qv <- pnorm(results_all$oslr_std_qv)
results_all$pts_oslr_qv <- 2 * pnorm(-abs(results_all$oslr_std_qv))
results_all$oslr_std_pqv <- results_all$oslr /
  sqrt(results_all$est_var_oslr_pqv)
results_all$pos_oslr_pqv <- pnorm(results_all$oslr_std_pqv)
results_all$pts_oslr_pqv <- 2 * pnorm(-abs(results_all$oslr_std_pqv))
results_all$oslr_std_Wu <- results_all$oslr / sqrt(results_all$est_var_oslr_Wu)
results_all$pos_oslr_Wu <- pnorm(results_all$oslr_std_Wu)
results_all$pts_oslr_Wu <- 2 * pnorm(-abs(results_all$oslr_std_Wu))

## Compute raw test statistics with estimates instead of true value
# Version 1: Use Weibull estimate
results_all$raw_test_mle <- results_all$oslr + results_all$hazard_diff
# Version 2: Use exponential estimate
results_all$raw_test_mle_exp <- results_all$oslr + results_all$hazard_diff_exp
# Version 2: Use non-parametric estimate
results_all$raw_test_na <- results_all$oslr + results_all$hazard_diff_np

## Compute test statistics without correction for MLE
# Version 1: Use Weibull estimates
results_all$mle_std_test_uncor_qv <- results_all$raw_test_mle /
  sqrt(results_all$est_var_oslr_qv)
results_all$pos_mle_uncor_qv <- pnorm(results_all$mle_std_test_uncor_qv)
results_all$pts_mle_uncor_qv <- 2 *
  pnorm(-abs(results_all$mle_std_test_uncor_qv))
results_all$mle_std_test_uncor_pqv <- results_all$raw_test_mle /
  sqrt(results_all$est_var_oslr_pqv_mle)
results_all$pos_mle_uncor_pqv <- pnorm(results_all$mle_std_test_uncor_pqv)
results_all$pts_mle_uncor_pqv <- 2 *
  pnorm(-abs(results_all$mle_std_test_uncor_pqv))
results_all$mle_std_test_uncor_Wu <- results_all$raw_test_mle /
  sqrt(results_all$est_var_oslr_Wu_mle)
results_all$pos_mle_uncor_Wu <- pnorm(results_all$mle_std_test_uncor_Wu)
results_all$pts_mle_uncor_Wu <- 2 *
  pnorm(-abs(results_all$mle_std_test_uncor_Wu))
# Version 2: Use exponential estimates
results_all$mle_exp_std_test_uncor_qv <- results_all$raw_test_mle_exp /
  sqrt(results_all$est_var_oslr_qv)
results_all$pos_mle_exp_uncor_qv <- pnorm(results_all$mle_exp_std_test_uncor_qv)
results_all$pts_mle_exp_uncor_qv <- 2 *
  pnorm(-abs(results_all$mle_exp_std_test_uncor_qv))
results_all$mle_exp_std_test_uncor_pqv <- results_all$raw_test_mle_exp /
  sqrt(results_all$est_var_oslr_pqv_mle_exp)
results_all$pos_mle_exp_uncor_pqv <- pnorm(
  results_all$mle_exp_std_test_uncor_pqv
)
results_all$pts_mle_exp_uncor_pqv <- 2 *
  pnorm(-abs(results_all$mle_exp_std_test_uncor_pqv))
results_all$mle_exp_std_test_uncor_Wu <- results_all$raw_test_mle_exp /
  sqrt(results_all$est_var_oslr_Wu_mle_exp)
results_all$pos_mle_exp_uncor_Wu <- pnorm(results_all$mle_exp_std_test_uncor_Wu)
results_all$pts_mle_exp_uncor_Wu <- 2 *
  pnorm(-abs(results_all$mle_exp_std_test_uncor_Wu))
# Version 3: Use non-parametric estimates
results_all$na_std_test_uncor_qv <- results_all$raw_test_na /
  sqrt(results_all$est_var_oslr_qv)
results_all$pos_na_uncor_qv <- pnorm(results_all$na_std_test_uncor_qv)
results_all$pts_na_uncor_qv <- 2 *
  pnorm(-abs(results_all$na_std_test_uncor_qv))
results_all$na_std_test_uncor_pqv <- results_all$raw_test_na /
  sqrt(results_all$est_var_oslr_pqv_na)
results_all$pos_na_uncor_pqv <- pnorm(
  results_all$na_std_test_uncor_pqv
)
results_all$pts_na_uncor_pqv <- 2 *
  pnorm(-abs(results_all$na_std_test_uncor_pqv))
results_all$na_std_test_uncor_Wu <- results_all$raw_test_na /
  sqrt(results_all$est_var_oslr_Wu_na)
results_all$pos_na_uncor_Wu <- pnorm(results_all$na_std_test_uncor_Wu)
results_all$pts_na_uncor_Wu <- 2 *
  pnorm(-abs(results_all$na_std_test_uncor_Wu))

## Compute test statistics with correction for MLE
# Version 1: Use Weibull estimates
results_all$mle_std_test_cor_qv <- results_all$raw_test_mle /
  sqrt(results_all$est_var_total_qv)
results_all$pos_mle_cor_qv <- pnorm(results_all$mle_std_test_cor_qv)
results_all$pts_mle_cor_qv <- 2 * pnorm(-abs(results_all$mle_std_test_cor_qv))
results_all$mle_std_test_cor_pqv <- results_all$raw_test_mle /
  sqrt(results_all$est_var_total_pqv)
results_all$pos_mle_cor_pqv <- pnorm(results_all$mle_std_test_cor_pqv)
results_all$pts_mle_cor_pqv <- 2 * pnorm(-abs(results_all$mle_std_test_cor_pqv))
results_all$mle_std_test_cor_Wu <- results_all$raw_test_mle /
  sqrt(results_all$est_var_total_Wu)
results_all$pos_mle_cor_Wu <- pnorm(results_all$mle_std_test_cor_Wu)
results_all$pts_mle_cor_Wu <- 2 * pnorm(-abs(results_all$mle_std_test_cor_Wu))
# Version 2: Use exponential estimates
results_all$mle_exp_std_test_cor_qv <- results_all$raw_test_mle_exp /
  sqrt(results_all$est_var_total_qv_exp)
results_all$pos_mle_exp_cor_qv <- pnorm(results_all$mle_exp_std_test_cor_qv)
results_all$pts_mle_exp_cor_qv <- 2 *
  pnorm(-abs(results_all$mle_exp_std_test_cor_qv))
results_all$mle_exp_std_test_cor_pqv <- results_all$raw_test_mle_exp /
  sqrt(results_all$est_var_total_pqv_exp)
results_all$pos_mle_exp_cor_pqv <- pnorm(results_all$mle_exp_std_test_cor_pqv)
results_all$pts_mle_exp_cor_pqv <- 2 *
  pnorm(-abs(results_all$mle_exp_std_test_cor_pqv))
results_all$mle_exp_std_test_cor_Wu <- results_all$raw_test_mle_exp /
  sqrt(results_all$est_var_total_Wu_exp)
results_all$pos_mle_exp_cor_Wu <- pnorm(results_all$mle_exp_std_test_cor_Wu)
results_all$pts_mle_exp_cor_Wu <- 2 *
  pnorm(-abs(results_all$mle_exp_std_test_cor_Wu))
# Version 3: Use non-parametric estimates
results_all$na_std_test_cor_qv <- results_all$raw_test_na /
  sqrt(results_all$est_var_total_qv_na)
results_all$pos_na_cor_qv <- pnorm(results_all$na_std_test_cor_qv)
results_all$pts_na_cor_qv <- 2 *
  pnorm(-abs(results_all$na_std_test_cor_qv))
results_all$na_std_test_cor_pqv <- results_all$raw_test_na /
  sqrt(results_all$est_var_total_pqv_na)
results_all$pos_na_cor_pqv <- pnorm(results_all$na_std_test_cor_pqv)
results_all$pts_na_cor_pqv <- 2 *
  pnorm(-abs(results_all$na_std_test_cor_pqv))
results_all$na_std_test_cor_Wu <- results_all$raw_test_na /
  sqrt(results_all$est_var_total_Wu_na)
results_all$pos_na_cor_Wu <- pnorm(results_all$na_std_test_cor_Wu)
results_all$pts_na_cor_Wu <- 2 *
  pnorm(-abs(results_all$na_std_test_cor_Wu))

# Compute two-sided p-values for two-sample log-rank test
# NOTE: Check if direction of tests correspond
results_all$pos_lr <- 1 - results_all$lr_p
results_all$pts_lr <- ifelse(
  results_all$pos_lr < 0.5,
  2 * results_all$pos_lr,
  2 * (1 - results_all$pos_lr)
)

# Set significance levels and compute confidence intervals
os_alpha <- 0.025
ts_alpha <- 2 * os_alpha

runs <- 100000

os_alpha_lb_ci <- os_alpha -
  qnorm(0.975) * sqrt(os_alpha * (1 - os_alpha) / runs)
os_alpha_ub_ci <- os_alpha +
  qnorm(0.975) * sqrt(os_alpha * (1 - os_alpha) / runs)

ts_alpha_lb_ci <- ts_alpha -
  qnorm(0.975) * sqrt(ts_alpha * (1 - ts_alpha) / runs)
ts_alpha_ub_ci <- ts_alpha +
  qnorm(0.975) * sqrt(ts_alpha * (1 - ts_alpha) / runs)

# Comparison of two-sided and one-sided empirical levels of
# 1) OSLR with correct reference curve
# 2) OSLR with correct reference curve (Wu's variance)
# 3) OSLR with estimated reference curve (no correction)
# 3.1) OSLR with estimated reference curve from exponential distribution if available (no correction)
# 3.2) OSLR with estimated reference curve from Nelson-Aalen estimation if available (no correction)
# 4) OSLR with estimated reference curve (no correction, Wu's variance)
# 4.1) OSLR with estimated reference curve from exponential distribution if available (no correction, Wu's variance)
# 4.2) OSLR with estimated reference curve from Nelson-Aalen estimation if available (no correction, Wu's variance)
# 5) Corrected test
# 5.1) Corrected test with estimates from exponential distribution if available
# 5.2) Corrected test with estimates from Nelson-Aalen estimation if available
# 6) Corrected test (Wu's variance for OSLR part)
# 6.1) Corrected test with estimates from exponential distribution if available (Wu's variance for OSLR part)
# 6.2) Corrected test with estimates from Nelson-Aalen estimation if available (Wu's variance for OSLR part)
# 7) Log-rank test

ts_rates <- aggregate(
  cbind(
    pts_oslr_pqv,
    pts_oslr_Wu,
    pts_mle_uncor_pqv,
    pts_mle_exp_uncor_pqv,
    pts_na_uncor_pqv,
    pts_mle_uncor_Wu,
    pts_mle_exp_uncor_Wu,
    pts_na_uncor_Wu,
    pts_mle_cor_pqv,
    pts_mle_exp_cor_pqv,
    pts_na_cor_pqv,
    pts_mle_cor_Wu,
    pts_mle_exp_cor_Wu,
    pts_na_cor_Wu,
    pts_lr
  ) ~ hr + alloc_ratio,
  data = results_all,
  FUN = function(x) mean(x <= ts_alpha)
)
# Convert to long format for ggplot
ts_rates_long <- pivot_longer(
  ts_rates,
  cols = c(
    pts_oslr_pqv,
    pts_oslr_Wu,
    pts_mle_uncor_pqv,
    pts_mle_exp_uncor_pqv,
    pts_na_uncor_pqv,
    pts_mle_uncor_Wu,
    pts_mle_exp_uncor_Wu,
    pts_na_uncor_Wu,
    pts_mle_cor_pqv,
    pts_mle_exp_cor_pqv,
    pts_na_cor_pqv,
    pts_mle_cor_Wu,
    pts_mle_exp_cor_Wu,
    pts_na_cor_Wu,
    pts_lr
  ),
  names_to = "Test",
  values_to = "rate"
)

os_left_rates <- aggregate(
  cbind(
    pos_oslr_pqv,
    pos_oslr_Wu,
    pos_mle_uncor_pqv,
    pos_mle_exp_uncor_pqv,
    pos_na_uncor_pqv,
    pos_mle_uncor_Wu,
    pos_mle_exp_uncor_Wu,
    pos_na_uncor_Wu,
    pos_mle_cor_pqv,
    pos_mle_exp_cor_pqv,
    pos_na_cor_pqv,
    pos_mle_cor_Wu,
    pos_mle_exp_cor_Wu,
    pos_na_cor_Wu,
    pos_lr
  ) ~ hr + alloc_ratio,
  data = results_all,
  FUN = function(x) mean(x <= os_alpha)
)
os_left_rates_long <- pivot_longer(
  os_left_rates,
  cols = c(
    pos_oslr_pqv,
    pos_oslr_Wu,
    pos_mle_uncor_pqv,
    pos_mle_exp_uncor_pqv,
    pos_na_uncor_pqv,
    pos_mle_uncor_Wu,
    pos_mle_exp_uncor_Wu,
    pos_na_uncor_Wu,
    pos_mle_cor_pqv,
    pos_mle_exp_cor_pqv,
    pos_na_cor_pqv,
    pos_mle_cor_Wu,
    pos_mle_exp_cor_Wu,
    pos_na_cor_Wu,
    pos_lr
  ),
  names_to = "Test",
  values_to = "rate"
)

os_right_rates <- aggregate(
  cbind(
    pos_oslr_pqv,
    pos_oslr_Wu,
    pos_mle_uncor_pqv,
    pos_mle_exp_uncor_pqv,
    pos_na_uncor_pqv,
    pos_mle_uncor_Wu,
    pos_mle_exp_uncor_Wu,
    pos_na_uncor_Wu,
    pos_mle_cor_pqv,
    pos_mle_exp_cor_pqv,
    pos_na_cor_pqv,
    pos_mle_cor_Wu,
    pos_mle_exp_cor_Wu,
    pos_na_cor_Wu,
    pos_lr
  ) ~ hr + alloc_ratio,
  data = results_all,
  FUN = function(x) mean((1 - x) <= os_alpha)
)
os_right_rates_long <- pivot_longer(
  os_right_rates,
  cols = c(
    pos_oslr_pqv,
    pos_oslr_Wu,
    pos_mle_uncor_pqv,
    pos_mle_exp_uncor_pqv,
    pos_na_uncor_pqv,
    pos_mle_uncor_Wu,
    pos_mle_exp_uncor_Wu,
    pos_na_uncor_Wu,
    pos_mle_cor_pqv,
    pos_mle_exp_cor_pqv,
    pos_na_cor_pqv,
    pos_mle_cor_Wu,
    pos_mle_exp_cor_Wu,
    pos_na_cor_Wu,
    pos_lr
  ),
  names_to = "Test",
  values_to = "rate"
)

rename_procedures <- function(my_df) {
  my_df$Test[my_df$Test %in% c("pts_lr", "pos_lr")] <- "TSLR"
  my_df$Test[
    my_df$Test %in% c("pts_mle_cor_pqv", "pos_mle_cor_pqv")
  ] <- "Corrected OSLR (pqv)"
  my_df$Test[
    my_df$Test %in% c("pts_mle_exp_cor_pqv", "pos_mle_exp_cor_pqv")
  ] <- "Corrected OSLR (exp. dist., pqv)"
  my_df$Test[
    my_df$Test %in% c("pts_na_cor_pqv", "pos_na_cor_pqv")
  ] <- "Corrected OSLR (non-par., pqv)"
  my_df$Test[
    my_df$Test %in% c("pts_mle_cor_Wu", "pos_mle_cor_Wu")
  ] <- "Corrected OSLR"
  my_df$Test[
    my_df$Test %in% c("pts_mle_exp_cor_Wu", "pos_mle_exp_cor_Wu")
  ] <- "Corrected OSLR (exp. dist.)"
  my_df$Test[
    my_df$Test %in% c("pts_na_cor_Wu", "pos_na_cor_Wu")
  ] <- "Corrected OSLR (non-par.)"
  my_df$Test[
    my_df$Test %in% c("pts_mle_uncor_pqv", "pos_mle_uncor_pqv")
  ] <- "Uncorrected OSLR (pqv)"
  my_df$Test[
    my_df$Test %in% c("pts_mle_exp_uncor_pqv", "pos_mle_exp_uncor_pqv")
  ] <- "Uncorrected OSLR (exp. dist., pqv)"
  my_df$Test[
    my_df$Test %in% c("pts_na_uncor_pqv", "pos_na_uncor_pqv")
  ] <- "Uncorrected OSLR (non-par., pqv)"
  my_df$Test[
    my_df$Test %in% c("pts_mle_uncor_Wu", "pos_mle_uncor_Wu")
  ] <- "Uncorrected OSLR"
  my_df$Test[
    my_df$Test %in% c("pts_mle_exp_uncor_Wu", "pos_mle_exp_uncor_Wu")
  ] <- "Uncorrected OSLR (exp. dist.)"
  my_df$Test[
    my_df$Test %in% c("pts_na_uncor_Wu", "pos_na_uncor_Wu")
  ] <- "Uncorrected OSLR (non-par.)"
  my_df$Test[my_df$Test %in% c("pts_oslr_pqv", "pos_oslr_pqv")] <- "OSLR (pqv)"
  my_df$Test[my_df$Test %in% c("pts_oslr_Wu", "pos_oslr_Wu")] <- "OSLR"

  return(my_df)
}

ts_rates_long <- rename_procedures(ts_rates_long)
os_left_rates_long <- rename_procedures(os_left_rates_long)
os_right_rates_long <- rename_procedures(os_right_rates_long)

p_vs_np <- c(
  "TSLR",
  "Corrected OSLR",
  "Corrected OSLR (non-par.)",
  "Uncorrected OSLR",
  "Uncorrected OSLR (non-par.)"
)

# Define tibble to show reference line only for T1E rate
hr_levels <- c("Type I error", "Power")
ref <- tibble(
  hr = factor(hr_levels[1], levels = hr_levels), # nur linkes Panel!
  y = os_alpha,
  ymin = os_alpha_lb_ci,
  ymax = os_alpha_ub_ci
)

# Define colours, linetypes and shapes for all methods
cols <- c(
  "Uncorrected OSLR" = "#C00000",
  "Uncorrected OSLR (non-par.)" = "#C00000",
  "Corrected OSLR" = "#009E73",
  "Corrected OSLR (non-par.)" = "#009E73",
  "TSLR" = "#333333"
)
np <- grepl("non-par", names(cols))
ltys <- setNames(ifelse(np, "dashed", "solid"), names(cols))
shps <- setNames(ifelse(np, 17, 16), names(cols)) # 17 = Dreieck, 16 = Punkt

p_vs_np_plot <-
  os_left_rates_long |>
  filter(Test %in% p_vs_np, hr %in% c(1, 0.8)) |>
  mutate(hr = factor(hr, levels = c(1, 0.8), labels = hr_levels)) |>
  ggplot(aes(x = alloc_ratio, y = rate)) +
  geom_rect(
    data = ref,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    alpha = 0.25
  ) +
  geom_hline(data = ref, aes(yintercept = y)) +
  geom_line(aes(colour = Test, linetype = Test), lwd = 1.5) +
  geom_point(aes(colour = Test, shape = Test), size = 3) +
  scale_colour_manual(values = cols, name = "Method") +
  scale_linetype_manual(values = ltys, name = "Method") +
  scale_shape_manual(values = shps, name = "Method") +
  facet_wrap(~hr, scales = "free_y") + # oder scales = "fixed"
  ylim(0, NA) +
  xlab(bquote(allocation ~ ratio ~ "(" * n[b] ~ "= 100)")) +
  ylab("rejection rate") +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme(
    strip.text = element_text(size = 14), # Facet-Überschriften
    axis.title = element_text(size = 14), # Achsentitel
    axis.text = element_text(size = 12), # Achsenlabels (Ticks)
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"), # fett
    legend.position = "bottom",
    legend.key.width = unit(1.6, "cm"),
    legend.key.height = unit(0.8, "cm")
  )
ggsave(
  p_vs_np_plot,
  filename = "results/plots/w_np/p_vs_np.pdf",
  width = 10,
  height = 5
)
