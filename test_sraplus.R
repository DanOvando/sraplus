library(tidyverse)

library(sraplus)

library(tmbstan)

sim <- sraplus_simulator(sigma_proc = 0.2, sigma_u = 0.1, q_slope = 0.025, r = 0.4, years = 50,q = 1e-3)

sim$pop %>% 
  ggplot(aes(year, depletion)) + 
  geom_point()

pop <- sim$pop

pop %>% 
  ggplot() + 
  geom_point(aes(year, scale(biomass))) + 
  geom_line(aes(year, scale(biomass * q), color = "index")) + 
  geom_line(aes(year, scale(catch / effort), color = "cpue"))

example_taxa <- "gadus morhua"


driors <- format_driors(
  taxa = example_taxa,
  catch = pop$catch,
  years = pop$year,
  initial_b = pop$depletion[1],
  initial_b_sd = 0.05,
  terminal_b = dplyr::last(pop$depletion),
  terminal_b_sd = 0.01,
  growth_rate = 0.4,
  growth_rate_cv = 0.1,
  use_heuristics = FALSE)


plot_driors(driors)

sir_fit <- fit_sraplus(driors = driors,
                       engine = 'sir',
                       draws = 5e4)

ml_driors <- format_driors(taxa = example_taxa,
                        catch = pop$catch,
                        years = pop$year,
                        index = pop$biomass * 1e-3,
                        index_years = pop$year,
                        initial_b = 1,
                        initial_b_sd = 0.05,
                        terminal_b = NA,
                        growth_rate = 0.4,
                        growth_rate_cv = 0.1)

plot_driors(ml_driors)


ml_fit <- fit_sraplus(driors = ml_driors,
                      engine = "tmb",
                      model = "sraplus_tmb", cleanup = FALSE)

plot_sraplus(ml_fit = ml_fit, years = ml_driors$years)


bayes_fit <- fit_sraplus(driors = ml_driors,
                      engine = "stan",
                      model = "sraplus_tmb", cleanup = FALSE,
                      n_keep = 4000)

plot_sraplus(ml_fit = ml_fit, bayes_fit = bayes_fit, sir_fit = sir_fit,years = ml_driors$years)


r_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_r",transformations = "exp") + 
  geom_vline(aes(xintercept = sim$params$r), color = "red")


m_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_m",transformations = "exp") + 
  geom_vline(aes(xintercept = sim$params$m), color = "red")

q_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_q",transformations = "exp") + 
  geom_vline(aes(xintercept = sim$pop$q[1]), color = "red")

sigma_proc_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_sigma_proc",transformations = "exp") + 
  geom_vline(aes(xintercept = sim$params$sigma_proc), color = "red")

plot(ml_fit$results$mean[ml_fit$results$variable == "r"], sim$params$r)
  abline(a = 0, b = 1)
  
ml_m_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_m",transformations = "exp") + 
  geom_vline(aes(xintercept = sim$params$m), color = "red")


bayes_fit$results %>% 
  filter(variable == "depletion") %>% 
  ggplot(aes(year, mean)) + 
  geom_line() + 
  geom_point(data = pop, aes(year, depletion))

ml_fit$results %>% 
  filter(variable == "depletion") %>% 
  ggplot(aes(1:length(mean), mean)) + 
  geom_line() + 
  geom_point(data = pop, aes(year, depletion))


bayes_fit$results %>% 
  filter(variable == "b_div_bmsy") %>% 
  ggplot(aes(year, mean)) + 
  geom_line() + 
  geom_point(data = pop, aes(year, b_bmsy))


# test effort -------------------------------------------------------------


  
effort_ml_driors <- format_driors(taxa = example_taxa,
                           catch = pop$catch,
                           years = pop$year,
                           effort = pop$effort,
                           index_years = pop$year,
                           initial_b = 1,
                           initial_b_sd = 0.05,
                           terminal_b = NA,
                           growth_rate = 0.4,
                           growth_rate_cv = 0.1,
                           q_slope = 0.025)


effort_ml_fit <- fit_sraplus(driors = effort_ml_driors,
                      engine = "tmb",
                      model = "sraplus_tmb", cleanup = FALSE)


effort_bayes_fit <- fit_sraplus(driors = effort_ml_driors,
                             engine = "stan",
                             model = "sraplus_tmb", cleanup = FALSE)

plot_sraplus(index_fit = ml_fit, effort_fit = effort_ml_fit, bayes_effort_fit =  effort_bayes_fit, years = ml_driors$years)


effort_ml_fit$results %>% 
  filter(variable == "depletion") %>% 
  ggplot(aes(1:length(mean), mean)) + 
  geom_line() + 
  geom_point(data = pop, aes(year, depletion))


ml_fit$results %>% 
  filter(variable == "depletion") %>% 
  ggplot(aes(1:length(mean), mean)) + 
  geom_line() + 
  geom_point(data = pop, aes(year, depletion))