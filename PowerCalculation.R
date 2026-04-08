library(ggplot2)
library(ggpubr)
library(tidyr)
library(patchwork)
library(cowplot)
library(grid)
library(gridExtra)

# Inputs you supply:
n <- 30

control_ssq <- 50 

# Full Data results
d_full_between <- 1.4208
d_full_long <- 1.2190
d_full_short <- 1.0568
d_full_long_conservative <- 1.3271
d_full_short_conservative <- 1.2752

# Delta-Only Results
d_delta_between <- 1.2796
d_delta_long <- 1.0134
d_delta_short <- 0.897
d_delta_long_conservative <- 1.1909
d_delta_short_conservative <- 1.1521

mean_diff_between_cond <- 0.4 # The experimental condition's average = 50 - 50 * 0.4 = 30

carryover_magnitude_short_seq = seq(from = 0.4, to = 0.9, by = 0.05)
carryover_magnitude_long_seq = seq(from = -0.05, to = -0.5, by = -0.05)
icc_seq = seq(from = 0.3, to = 0.8, by = 0.1)


# Estimates the power of an experiment given two conditions: control and experimental
estimate_power_between <- function(mean_diff_between_cond, outcome_effect_size) {
  muD <- mean_diff_between_cond * control_ssq # mean difference in change: (A_post-A_pre)-(B_post-B_pre)
  # print(muD)
  sd_dA <- (control_ssq) / outcome_effect_size       # SD of SSQ change in control condition
  sd_dB <- (control_ssq - muD) / outcome_effect_size        # SD of SSQ change in experimental condition
  
  sdD = (sd_dA+sd_dB)/2
  dz = muD / sdD
    
  # Between-subjects power
  power <- power.t.test(n = n, delta = muD, sd = (sd_dA+sd_dB)/2, type = "two.sample", alternative = "two.sided", sig.level = 0.05)$power # for non-paired
    
  return (list(sdD = sdD, dz = dz, power = power))
}

estimate_power_within <- function(
    carryover_magnitude, # the amount of carryover effect relative to the main effect (Me/M)
    mean_diff_between_cond, # the mean difference between control and experimental conditions relative to the control condition
    outcome_effect_size, # the effect size of the estimated per condition effect size of the experiment outcome
    r_delta, # the interclass correlation coefficient between the control and experimental conditions; could use different
    within_or_between # within- or between-subjects
    ) {
  
  muD <- (1-carryover_magnitude*0.5) * mean_diff_between_cond * control_ssq # mean difference in change: (A_post-A_pre)-(B_post-B_pre)
  # print(muD)
  sd_dA <- (control_ssq) / outcome_effect_size       # SD of change in A
  sd_dB <- (control_ssq - muD) / outcome_effect_size        # SD of change in B

  
  sdD <- sqrt(sd_dA^2 + sd_dB^2 - 2*r_delta*sd_dA*sd_dB)
  dz <- muD / sdD
  
  if (within_or_between) {
    # Within-subjects power
    power <- power.t.test(n = n, delta = muD, sd = sdD, type = "paired", alternative = "two.sided", sig.level = 0.05)$power # for paired tests
  } else {
    # Between-subjects power
    power <- power.t.test(n = n, delta = muD, sd = (sd_dA+sd_dB)/2, type = "two.sample", alternative = "two.sided", sig.level = 0.05)$power # for non-paired
  }
  
  return (list(sdD = sdD, dz = dz, power = power))
}

# Estimate the power of between subjects design for full and delta results
power_full_between <- estimate_power_between(mean_diff_between_cond, d_full_between)$power
power_full_between
power_delta_between <- estimate_power_between(mean_diff_between_cond, d_delta_between)$power
power_delta_between
  
column_names <- c("OutcomeEffectSize", "StudyType", "CarryoverMagnitude", "ICC", "Power")

all_estimation <- function(d, 
                            seq_carryover_magnitude, 
                            seq_icc_delta) {
  power_results <- data.frame(matrix(nrow = 0, ncol = length(column_names)))
  colnames(power_results) <- column_names
  
  for (carryover_magnitude in seq_carryover_magnitude) {
    for (r_delta in seq_icc_delta) {
      power_estimation = estimate_power_within(carryover_magnitude, mean_diff_between_cond, d, r_delta, TRUE)
      # cat("Short, Carryover_magnitude: ", carryover_magnitude, ",r_delta: ", r_delta, ",Power: ", power_estimation$power, "\n")
      new_row <- data.frame(OutcomeEffectSize = d, StudyType = "Short", CarryoverMagnitude = carryover_magnitude, ICC = r_delta, Power = power_estimation$power)
      
      power_results <- rbind(power_results, new_row)
    }
  }
  
  return (power_results)
}


visualize_estimations <- function(power_results, between_line, carryover_magnitude_seq, figure_title) {
  power_results$ICC <- factor(power_results$ICC)
  
  is_desc <- length(carryover_magnitude_seq) >= 2 &&
    carryover_magnitude_seq[1] > carryover_magnitude_seq[length(carryover_magnitude_seq)]
  
  x_limits <- range(carryover_magnitude_seq, na.rm = TRUE)  # always c(min, max)
  
  p <- ggplot(power_results,
              aes(x = CarryoverMagnitude, y = Power, color = ICC, group = ICC)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) + 
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE
      )
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    labs(
      x = NULL,
      y = NULL,
      color = "ICC",
      title = figure_title
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    geom_hline(yintercept = between_line, linetype = "dashed", color = "grey40")
  
  # Apply exactly one x scale, then return it
  if (is_desc) {
    p + scale_x_reverse(limits = x_limits, breaks = carryover_magnitude_seq)
  } else {
    p + scale_x_continuous(limits = x_limits, breaks = carryover_magnitude_seq)
  }
  
  return(p)
}


power_estimate_and_visualize <- function(d, carryover_seq, power_between, title) {
  power_results <- all_estimation(d, carryover_seq, icc_seq)
  visualize_estimations(power_results, power_between, carryover_seq, title)
}

create_all_plots <- function() {
  # Power estimations
  # Full short
  plot_full_short <- power_estimate_and_visualize(
    d_full_short, 
    carryover_magnitude_short_seq, 
    power_full_between,
    "Full Data, Short"
    )
  
  # Full short conservative
  plot_full_short_conservative <- power_estimate_and_visualize(
    d_full_short_conservative, carryover_magnitude_short_seq, power_full_between,
                                                               "Full Data, Short (Conservative)")
  
  # Full long
  plot_full_long <- power_estimate_and_visualize(
    d_full_long, carryover_magnitude_long_seq, power_full_between,
                                                 "Full Data, Long")
  
  # Full long conservative
  plot_full_long_conservative <- power_estimate_and_visualize(
    d_full_long_conservative, carryover_magnitude_long_seq, power_full_between,
                                                              "Full Data, Long (Conservative)")
  
  
  # Delta short 
  plot_delta_short <- power_estimate_and_visualize(
    d_delta_short, carryover_magnitude_short_seq, power_delta_between,
    "Delta Only, Short")
  
  # Delta short conservative
  plot_delta_short_conservative <- power_estimate_and_visualize(
    d_delta_short_conservative, carryover_magnitude_short_seq, power_delta_between,
    "Delta Only, Short (Conservative)")
  
  # Delta long
  plot_delta_long <- power_estimate_and_visualize(
    d_delta_long, carryover_magnitude_long_seq, power_delta_between,
    "Delta Only, Long")
  
  # Delta long conservative
  plot_delta_long_conservative <- power_estimate_and_visualize(
    d_delta_long_conservative, carryover_magnitude_long_seq, power_delta_between,
    "Delta Only, Long (Conservative)")
  
  # --- ensure NO per-panel axis titles; keep legend only in one plot initially ---
  strip_axes <- theme(axis.title.x = element_blank(),
                      axis.title.y = element_blank())
  

  # ---- remove per-panel axis titles (keep ticks) ----
  strip_axes <- theme(axis.title.x = element_blank(),
                      axis.title.y = element_blank())
  
  plots <- list(
    plot_full_short + strip_axes,
    plot_full_short_conservative + strip_axes,
    plot_full_long + strip_axes,
    plot_full_long_conservative + strip_axes,
    plot_delta_short + strip_axes,
    plot_delta_short_conservative + strip_axes,
    plot_delta_long + strip_axes,
    plot_delta_long_conservative + strip_axes
  )
  
  # ---- Arrange 2 rows × 4 columns with ONE shared legend ----
  fig <- ggarrange(
    plotlist = plots,
    nrow = 2, ncol = 4,
    labels = letters[1:8],        # (a) ... (h)
    common.legend = TRUE,
    legend = "bottom",
    align = "hv"
  )
  
  # ---- Add global axis titles (one copy only) ----
  fig2 <- annotate_figure(
    fig,
    left = text_grob("Statistical Power", rot = 90, size = 15),
    bottom = text_grob("Carryover Magnitude", size = 15)
  )
  
  fig2
}

create_all_plots()
