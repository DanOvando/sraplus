
library(sraplus)
library(ggplot2)
library(dplyr)
library(tmbstan)
example_taxa <- "gadus morhua"

data <-
  readr::read_csv(here::here("data",
                             "Data_Effort_CPUE_forRH.csv")) %>% 
  janitor::clean_names()

head(data) 

effort <- data %>% 
  filter(region == "North America") %>% 
  group_by(year) %>% 
  summarise(effort = sum(nominal_effort_k_w_days)) %>% 
  filter(year %in% cod$year)


driors <- format_driors(
  taxa = example_taxa,
  catch = cod$catch,
  years = cod$year,
  effort = effort$effort,
  initial_b = 1,
  initial_b_sd = 0.2,
  terminal_b = NA,
  terminal_b_sd = 2,
  index_years = effort$year,
  use_heuristics = FALSE,
  q_slope = 0.025)



ml_fit <- fit_sraplus(
  driors = driors,
  engine = "tmb",
  model = "sraplus_tmb",
  cleanup = FALSE,
  chains = 4,
  cores = 4
)

plot_sraplus(ml_fit, years = driors$years)

bayes_fit <- fit_sraplus(
  driors = driors,
  engine = "stan",
  model = "sraplus_tmb",
  cleanup = FALSE,
  chains = 4,
  cores = 4
)

rstanarm::launch_shinystan(bayes_fit$fit)

plot_sraplus( bayes_fit = bayes_fit, years = driors$years)


plot_sraplus(ml_fit = ml_fit, bayes_fit = bayes_fit, years = driors$years)


plot(scale(ml_fit$results$mean[ml_fit$results$variable == "index_hat_t"]))
lines(scale(cod$catch/effort$effort))

# hello
