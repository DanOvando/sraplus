library(ggplot2)
library(tidyr)
library(dplyr)
library(sraplus)


catch_only_driors <- sraplus::format_driors(
  taxa = example_taxa,
  catch = cod$catch,
  years = cod$year,
  use_heuristics = FALSE,
  terminal_state = 0.5,
  terminal_state_cv = 0.25
)


catch_only_fit <- fit_sraplus(driors = catch_only_driors,
                              engine = "sir",
                              draws = 1e5,
                              n_keep = 2000,
                              estimate_proc_error = FALSE, 
                              estimate_shape = TRUE,
                              tune_prior_predictive = TRUE)


plot_sraplus(catch_only_fit)


plot_prior_posterior(catch_only_fit, catch_only_driors)



sar_driors <- format_driors(
  taxa = example_taxa,
  catch = cod$catch,
  years = cod$year,
  initial_state = 1,
  initial_state_cv = 0.25,
  use_heuristics = FALSE,
  sar = 10,
  sar_cv = NA,
  use_b_reg = FALSE,
  b_ref_type = "k")

sar_fit <- fit_sraplus(driors = sar_driors,
                              engine = "sir",
                              draws = 1e5,
                              n_keep = 2000,
                              estimate_proc_error = FALSE, 
                              estimate_shape = TRUE,
                              tune_prior_predictive = TRUE)

plot_prior_posterior(sar_fit, sar_driors)


sar_driors_2 <- format_driors(
  taxa = example_taxa,
  catch = cod$catch,
  years = cod$year,
  initial_state = 1,
  initial_state_cv = 0.25,
  use_heuristics = FALSE,
  sar = 10,
  sar_cv = NA,
  use_b_reg = TRUE,
  b_ref_type = "k")

sar_fit_2 <- fit_sraplus(driors = sar_driors_2,
                        engine = "sir",
                        draws = 1e5,
                        n_keep = 2000,
                        estimate_proc_error = FALSE, 
                        estimate_shape = TRUE,
                        tune_prior_predictive = FALSE)

plot_sraplus(sar_fit, sar_fit_2)

diagnose_sraplus(sar_fit, sar_driors)

plot_prior_posterior(sar_fit, sar_driors)


fmi_sar_driors <- format_driors(
  taxa = example_taxa,
  catch = cod$catch,
  years = cod$year,
  initial_state = 1,
  initial_state_cv = 0.25,
  use_heuristics = FALSE,
  sar = 10,
  fmi = c("research" = 1, "management" = 1, "socioeconomics" = 1, 'enforcement' = 1),
  sar_cv = NA,
  use_b_reg = FALSE,
  b_ref_type = "k")

sraplus::plot_driors(fmi_sar_driors)

fmi_sar_fit <- fit_sraplus(
  driors = fmi_sar_driors,
  engine = "sir",
  draws = 1e6,
  n_keep = 2000,
  estimate_shape = FALSE,
  estimate_proc_error = FALSE
)

plot_sraplus(fmi_sar = fmi_sar_fit,
             sar = sar_fit)



