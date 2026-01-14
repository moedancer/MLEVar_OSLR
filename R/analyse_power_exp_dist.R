### Script to analyse power results from single scenario files with exponential distribution

library(ggplot2)
library(tidyr)
library(ggpubr)

# Load and combine results from all single scenario files
results_all <- NULL

for (tmp_file in list.files(
  path = "results/single_scenarios",
  pattern = "t2e_raw_KAPPA1_"
)) {
  load(paste("results/single_scenarios/", tmp_file, sep = ""))
  results_all <- rbind(results_all, results)
}

# Compute missing sample size
results_all$n_a <- results_all$n_b / results_all$alloc_ratio

### Compute different variance estimates
## Compute Wu's variance with estimated parameters
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
## Compute total variance
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

## Compute raw test statistics with MLE estimate instead of true value
# Version 1: Use Weibull estimate
results_all$raw_test_mle <- results_all$oslr + results_all$hazard_diff
# Version 2: Use exponential estimate
results_all$raw_test_mle_exp <- results_all$oslr + results_all$hazard_diff_exp

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
# 4) OSLR with estimated reference curve (no correction, Wu's variance)
# 4.1) OSLR with estimated reference curve from exponential distribution if available (no correction, Wu's variance)
# 5) Corrected test
# 5.1) Corrected test with estimates from exponential distribution if available
# 6) Corrected test (Wu's variance for OSLR part)
# 6.1) Corrected test with estimates from exponential distribution if available (Wu's variance for OSLR part)
# 7) Log-rank test

ts_rates <- aggregate(
  cbind(
    pts_oslr_pqv,
    pts_oslr_Wu,
    pts_mle_uncor_pqv,
    pts_mle_exp_uncor_pqv,
    pts_mle_uncor_Wu,
    pts_mle_exp_uncor_Wu,
    pts_mle_cor_pqv,
    pts_mle_exp_cor_pqv,
    pts_mle_cor_Wu,
    pts_mle_exp_cor_Wu,
    pts_lr
  ) ~ n_b + alloc_ratio + shape,
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
    pts_mle_uncor_Wu,
    pts_mle_exp_uncor_Wu,
    pts_mle_cor_pqv,
    pts_mle_exp_cor_pqv,
    pts_mle_cor_Wu,
    pts_mle_exp_cor_Wu,
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
    pos_mle_uncor_Wu,
    pos_mle_exp_uncor_Wu,
    pos_mle_cor_pqv,
    pos_mle_exp_cor_pqv,
    pos_mle_cor_Wu,
    pos_mle_exp_cor_Wu,
    pos_lr
  ) ~ n_b + alloc_ratio + shape,
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
    pos_mle_uncor_Wu,
    pos_mle_exp_uncor_Wu,
    pos_mle_cor_pqv,
    pos_mle_exp_cor_pqv,
    pos_mle_cor_Wu,
    pos_mle_exp_cor_Wu,
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
    pos_mle_uncor_Wu,
    pos_mle_exp_uncor_Wu,
    pos_mle_cor_pqv,
    pos_mle_exp_cor_pqv,
    pos_mle_cor_Wu,
    pos_mle_exp_cor_Wu,
    pos_lr
  ) ~ n_b + alloc_ratio + shape,
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
    pos_mle_uncor_Wu,
    pos_mle_exp_uncor_Wu,
    pos_mle_cor_pqv,
    pos_mle_exp_cor_pqv,
    pos_mle_cor_Wu,
    pos_mle_exp_cor_Wu,
    pos_lr
  ),
  names_to = "Test",
  values_to = "rate"
)

### GO ON HERE!

rename_procedures <- function(my_df) {
  my_df$Test[my_df$Test %in% c("pts_lr", "pos_lr")] <- "TSLR"
  my_df$Test[
    my_df$Test %in% c("pts_mle_cor_pqv", "pos_mle_cor_pqv")
  ] <- "Corrected OSLR"
  my_df$Test[
    my_df$Test %in% c("pts_mle_exp_cor_pqv", "pos_mle_exp_cor_pqv")
  ] <- "Corrected OSLR (exp. dist.)"
  my_df$Test[
    my_df$Test %in% c("pts_mle_cor_Wu", "pos_mle_cor_Wu")
  ] <- "Corrected OSLR (Wu)"
  my_df$Test[
    my_df$Test %in% c("pts_mle_exp_cor_Wu", "pos_mle_exp_cor_Wu")
  ] <- "Corrected OSLR (exp. dist., Wu)"
  my_df$Test[
    my_df$Test %in% c("pts_mle_uncor_pqv", "pos_mle_uncor_pqv")
  ] <- "Uncorrected OSLR"
  my_df$Test[
    my_df$Test %in% c("pts_mle_exp_uncor_pqv", "pos_mle_exp_uncor_pqv")
  ] <- "Uncorrected OSLR (exp. dist.)"
  my_df$Test[
    my_df$Test %in% c("pts_mle_uncor_Wu", "pos_mle_uncor_Wu")
  ] <- "Uncorrected OSLR (Wu)"
  my_df$Test[
    my_df$Test %in% c("pts_mle_exp_uncor_Wu", "pos_mle_exp_uncor_Wu")
  ] <- "Uncorrected OSLR (exp. dist., Wu)"
  my_df$Test[my_df$Test %in% c("pts_oslr_pqv", "pos_oslr_pqv")] <- "OSLR"
  my_df$Test[my_df$Test %in% c("pts_oslr_Wu", "pos_oslr_Wu")] <- "OSLR (Wu)"

  return(my_df)
}

ts_rates_long <- rename_procedures(ts_rates_long)
os_left_rates_long <- rename_procedures(os_left_rates_long)
os_right_rates_long <- rename_procedures(os_right_rates_long)

plotdir <- "results/plots"

n_b_vec <- unique(results_all$n_b)
alloc_ratio_vec <- unique(results_all$alloc_ratio)
kappa_vec <- unique(results_all$shape)

for (kappa_temp in kappa_vec) {
  for (alloc_ratio_temp in alloc_ratio_vec) {
    for (n_b_temp in n_b_vec) {
      ts_rates_nfix_kfix_temp <-
        ggplot(
          data = ts_rates_long[
            ts_rates_long$shape == kappa_temp &
              ts_rates_long$n_b == n_b_temp,
          ],
          aes(x = alloc_ratio, y = rate)
        ) +
        geom_line(aes(colour = Test)) +
        geom_point(aes(colour = Test)) +
        geom_abline(intercept = ts_alpha, slope = 0) +
        annotate(
          "rect",
          xmin = -Inf,
          xmax = Inf,
          ymin = ts_alpha_lb_ci,
          ymax = ts_alpha_ub_ci,
          alpha = 0.25
        ) +
        ylim(0, NA) +
        xlab(bquote(
          allocation ~ ratio ~ "(" * n[b] ~ "=" ~ .(n_b_temp) * ")"
        )) +
        ylab("rejection rate")
      ggsave(
        filename = paste(
          "n_kappa_fixed/ts_rates_n",
          n_b_temp,
          "_k",
          sub(x = kappa_temp, pattern = "\\.", replacement = "dec"),
          ".pdf",
          sep = ""
        ),
        path = plotdir,
        plot = ts_rates_nfix_kfix_temp,
        device = "pdf",
        height = 5,
        width = 7
      )

      ts_rates_pifix_kfix_temp <-
        ggplot(
          data = ts_rates_long[
            ts_rates_long$shape == kappa_temp &
              ts_rates_long$alloc_ratio == alloc_ratio_temp,
          ],
          aes(x = n_b, y = rate)
        ) +
        geom_line(aes(colour = Test)) +
        geom_point(aes(colour = Test)) +
        geom_abline(intercept = ts_alpha, slope = 0) +
        annotate(
          "rect",
          xmin = -Inf,
          xmax = Inf,
          ymin = ts_alpha_lb_ci,
          ymax = ts_alpha_ub_ci,
          alpha = 0.25
        ) +
        ylim(0, NA) +
        xlab(bquote(
          n[b] ~ "(" * allocation ~ ratio ~ "=" ~ .(alloc_ratio_temp) * ")"
        )) +
        ylab("rejection rate")
      ggsave(
        filename = paste(
          "alloc_kappa_fixed/ts_rates_pi",
          sub(x = alloc_ratio_temp, pattern = "\\.", replacement = "dec"),
          "_k",
          sub(x = kappa_temp, pattern = "\\.", replacement = "dec"),
          ".pdf",
          sep = ""
        ),
        path = plotdir,
        plot = ts_rates_pifix_kfix_temp,
        device = "pdf",
        height = 5,
        width = 7
      )

      ts_rates_kfix_temp <-
        ggarrange(
          ts_rates_nfix_kfix_temp,
          ts_rates_pifix_kfix_temp,
          ncol = 2,
          nrow = 1,
          common.legend = TRUE,
          legend = "right"
        )
      ggsave(
        filename = paste(
          "combined/ts_rates_k",
          sub(x = kappa_temp, pattern = "\\.", replacement = "dec"),
          "_n",
          n_b_temp,
          "_pi",
          sub(x = alloc_ratio_temp, pattern = "\\.", replacement = "dec"),
          ".pdf",
          sep = ""
        ),
        path = plotdir,
        plot = ts_rates_kfix_temp,
        device = "pdf",
        height = 5,
        width = 12
      )

      os_left_rates_nfix_kfix_temp <-
        ggplot(
          data = os_left_rates_long[
            os_left_rates_long$shape == kappa_temp &
              os_left_rates_long$n_b == n_b_temp,
          ],
          aes(x = alloc_ratio, y = rate)
        ) +
        geom_line(aes(colour = Test)) +
        geom_point(aes(colour = Test)) +
        geom_abline(intercept = os_alpha, slope = 0) +
        annotate(
          "rect",
          xmin = -Inf,
          xmax = Inf,
          ymin = os_alpha_lb_ci,
          ymax = os_alpha_ub_ci,
          alpha = 0.25
        ) +
        ylim(0, NA) +
        xlab(bquote(
          allocation ~ ratio ~ "(" * n[b] ~ "=" ~ .(n_b_temp) * ")"
        )) +
        ylab("rejection rate")
      ggsave(
        filename = paste(
          "n_kappa_fixed/os_left_rates_n",
          n_b_temp,
          "_k",
          sub(x = kappa_temp, pattern = "\\.", replacement = "dec"),
          ".pdf",
          sep = ""
        ),
        path = plotdir,
        plot = os_left_rates_nfix_kfix_temp,
        device = "pdf",
        height = 5,
        width = 7
      )

      os_left_rates_pifix_kfix_temp <-
        ggplot(
          data = os_left_rates_long[
            os_left_rates_long$shape == kappa_temp &
              os_left_rates_long$alloc_ratio == alloc_ratio_temp,
          ],
          aes(x = n_b, y = rate)
        ) +
        geom_line(aes(colour = Test)) +
        geom_point(aes(colour = Test)) +
        geom_abline(intercept = os_alpha, slope = 0) +
        annotate(
          "rect",
          xmin = -Inf,
          xmax = Inf,
          ymin = os_alpha_lb_ci,
          ymax = os_alpha_ub_ci,
          alpha = 0.25
        ) +
        ylim(0, NA) +
        xlab(bquote(
          n[b] ~ "(" * allocation ~ ratio ~ "=" ~ .(alloc_ratio_temp) * ")"
        )) +
        ylab("rejection rate")
      ggsave(
        filename = paste(
          "alloc_kappa_fixed/os_left_rates_pi",
          sub(x = alloc_ratio_temp, pattern = "\\.", replacement = "dec"),
          "_k",
          sub(x = kappa_temp, pattern = "\\.", replacement = "dec"),
          ".pdf",
          sep = ""
        ),
        path = plotdir,
        plot = os_left_rates_pifix_kfix_temp,
        device = "pdf",
        height = 5,
        width = 7
      )

      os_left_rates_kfix_temp <-
        ggarrange(
          os_left_rates_nfix_kfix_temp,
          os_left_rates_pifix_kfix_temp,
          ncol = 2,
          nrow = 1,
          common.legend = TRUE,
          legend = "right"
        )
      ggsave(
        filename = paste(
          "combined/os_left_rates_k",
          sub(x = kappa_temp, pattern = "\\.", replacement = "dec"),
          "_n",
          n_b_temp,
          "_pi",
          sub(x = alloc_ratio_temp, pattern = "\\.", replacement = "dec"),
          ".pdf",
          sep = ""
        ),
        path = plotdir,
        plot = os_left_rates_kfix_temp,
        device = "pdf",
        height = 5,
        width = 12
      )
    }
  }
}
