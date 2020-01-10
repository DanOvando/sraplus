library(tidyverse)

library(sraplus)

library(tmbstan)

rstan_options(auto_write = TRUE)

set.seed(42)
sim <-
  sraplus_simulator(
    sigma_proc = 0,
    sigma_u = 0,
    q_slope = 0,
    r = 0.4,
    years = 25,
    q = 1e-3,
    m = 1.01,
    init_u_umsy = 1
  )

sim$pop %>% 
  ggplot(aes(year, u)) + 
  geom_point()
sim$pop %>% 
  ggplot(aes(year, biomass / 1000)) + 
  geom_point()

pop <- sim$pop

pop %>% 
  ggplot() + 
  geom_point(aes(year, scale(biomass))) + 
  geom_line(aes(year, scale(biomass * q), color = "index")) + 
  geom_line(aes(year, scale(catch / effort), color = "cpue"))

example_taxa <- "gadus sdfg"


driors <- format_driors(
  taxa =
    example_taxa,
  catch = pop$catch,
  years = pop$year,
  initial_state = pop$depletion[1],
  initial_state_cv = 0.05,
  terminal_state = dplyr::last(pop$depletion),
  terminal_state_cv = 0.01,
  growth_rate_prior = 0.4,
  growth_rate_prior_cv = 0.1,
  use_heuristics = FALSE
)


plot_driors(driors)

sir_fit <- fit_sraplus(driors = driors,
                       engine = 'sir',
                       draws = 1e5,
                       estimate_k = FALSE)

sir_diagnostics <- diagnose_sraplus(sir_fit, driors)

ml_driors <- format_driors(taxa = example_taxa,
                        catch = pop$catch,
                        years = pop$year,
                        index = pop$biomass * 1e-3 * exp(rnorm(length(pop$biomass),0,.4)),
                        index_years = pop$year,
                        initial_state = 1,
                        initial_state_cv = 0.05,
                        terminal_state = NA,
                        growth_rate_prior = 0.4,
                        growth_rate_prior_cv = 0.5,
                        sigma_r_prior = 0.05,
                        sigma_r_prior_cv = 0.1)

plot_driors(ml_driors)


ml_fit <- fit_sraplus(driors = ml_driors,
                      engine = "tmb",
                      model = "sraplus_tmb",
                      estimate_shape = FALSE, 
                      estimate_proc_error = FALSE,
                      estimate_k = FALSE,
                      eps = 1e-6)

diagnose_sraplus(ml_fit, ml_driors)

plot_sraplus(ml_fit = ml_fit, years = ml_driors$years)

plot_prior_posterior(ml_fit, ml_driors)



bayes_fit <- fit_sraplus(driors = ml_driors,
                      engine = "stan",
                      model = "sraplus_tmb",
                      n_keep = 2000,
                      estimate_shape = TRUE,
                      estimate_proc_error = TRUE)

diagnose_sraplus(bayes_fit, ml_driors)

test <- names(bayes_fit$fit)

# pair_vars <- test[!str_detect(test,"uc_proc")]
# 
# pairs(bayes_fit$fit, pars = pair_vars)


plot_sraplus(ml_fit = ml_fit, bayes_fit = bayes_fit,years = ml_driors$years)


# 
# a = tidybayes::gather_draws(bayes_fit$fit,log_f_t[year]) %>% 
#   mutate(u = 1 / (1 + exp(-.value)))
# 
# a %>% 
#   ggplot(aes(year, u)) + 
#   geom_smooth() + 
#   geom_point(data = pop,aes(year,effort * q))

r_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_r",transformations = "exp") + 
  geom_vline(aes(xintercept = sim$params$r), color = "red")


# m_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_m",transformations = "exp") +
#   geom_vline(aes(xintercept = sim$params$m), color = "red")

q_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_q",transformations = "exp") + 
  geom_vline(aes(xintercept = sim$pop$q[1]), color = "red")

sigma_proc_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_sigma_proc",transformations = "exp") +
  geom_vline(aes(xintercept = sim$params$sigma_proc), color = "red")

# plot(ml_fit$results$mean[ml_fit$results$variable == "r"], sim$params$r)
#   abline(a = 0, b = 1)
  
# ml_m_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_m",transformations = "exp") + 
#   geom_vline(aes(xintercept = sim$params$m), color = "red")


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


# plot prior-posterior ----------------------------------------------------



log_sigma_obs_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_sigma_obs", freq = FALSE) + 
  geom_density(data = data_frame(log_sigma_obs = rnorm(1000,-3,.1)), aes(log_sigma_obs), fill = "lightgrey",
               alpha = 0.5)


log_sigma_proc_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_sigma_proc", freq = FALSE) + 
  geom_density(data = data_frame(log_sigma_proc = rnorm(1000,-3,.5)), aes(log_sigma_proc), fill = "lightgrey",
               alpha = 0.5)

log_q_guess <- log(mean(ml_driors$index / ml_driors$catch))


log_q_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_q", freq = FALSE) + 
  geom_density(data = data_frame(log_q = rnorm(1000,log_q_guess,.2)), aes(log_q), fill = "lightgrey",
               alpha = 0.5)


# log_m_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_m", freq = FALSE) + 
#   geom_density(data = data_frame(log_m = rnorm(1000,.1,.5)), aes(log_m), fill = "lightgrey",
#                alpha = 0.5)



# test effort -------------------------------------------------------------

set.seed(42)

sim <-
  sraplus_simulator(
    sigma_proc = 0.05,
    sigma_u = 0.1,
    q_slope = 0.05,
    r = 0.4,
    years = 25,
    q = 1e-3,
    m = 1.01,
    init_u_umsy = 1
  )

pop <- sim$pop  

sim$pop %>% 
  ggplot(aes(year, depletion)) + 
  geom_point()

sim$pop %>% 
  ggplot(aes(year, catch)) + 
  geom_point()

pop %>% 
  ggplot() + 
  geom_point(aes(year, scale(biomass))) + 
  geom_line(aes(year, scale(biomass * 1e-3), color = "index")) + 
  geom_line(aes(year, scale(catch / effort), color = "cpue")) + 
  geom_line(aes(year, scale(u), color = "u"))



effort_ml_driors <- format_driors(taxa = example_taxa,
                           catch = pop$catch,
                           years = pop$year,
                           effort = pop$effort,
                           index_years = pop$year,
                           initial_state = 1,
                           initial_state_cv = 0.05,
                           terminal_state = NA,
                           growth_rate = 0.4,
                           growth_rate_cv = 0.1,
                           q_slope = 0.025,
                           q_slope_cv = 0.01)

index_ml_driors <- format_driors(taxa = example_taxa,
                                  catch = pop$catch,
                                  years = pop$year,
                                  index = pop$biomass * 1e-3,
                                  index_years = pop$year,
                                  initial_state = 1,
                                  initial_state_cv = 0.05,
                                  terminal_state = NA,
                                  growth_rate = 0.4,
                                  growth_rate_cv = 0.1,
                                 sigma_r = 0.05,
                                 sigma_r_cv = 0.05)


index_bayes_fit <- fit_sraplus(driors = index_ml_driors,
                               engine = "stan",
                               model = "sraplus_tmb",
                               n_keep = 4000)


index_ml_fit <- fit_sraplus(driors = index_ml_driors,
                            engine = "tmb",
                            model = "sraplus_tmb")

plot_sraplus(
  index_ml = index_ml_fit,
  index_bayes = index_bayes_fit,
  years = ml_driors$years
)


effort_ml_fit <- fit_sraplus(driors = effort_ml_driors,
                      engine = "tmb",
                      model = "sraplus_tmb",
                      estimate_qslope = TRUE)


effort_bayes_fit <- fit_sraplus(driors = effort_ml_driors,
                             engine = "stan",
                             model = "sraplus_tmb",
                             adapt_delta = 0.9,
                             max_treedepth = 10,
                             n_keep = 6000,
                             chains = 2, 
                             cores = 2,
                             estimate_qslope = FALSE)


effort_q_hat <- bayesplot::mcmc_hist(as.matrix(effort_bayes_fit$fit), "log_q",transformations = "exp") + 
  geom_vline(aes(xintercept = sim$pop$q[1]), color = "red")


# a = tidybayes::gather_draws(effort_bayes_fit$fit,log_f_t[year]) %>% 
#   mutate(u = 1 / (1 + exp(-.value)))
# 
# a %>% 
#   ggplot(aes(year, u)) + 
#   geom_smooth() + 
#   geom_point(data = pop,aes(year,effort * q))


plot_sraplus(
  effort_bayes =  effort_bayes_fit,
  years = effort_ml_driors$years
)


plot_sraplus(
  index_bayes = index_bayes_fit,
  effort_bayes =  effort_bayes_fit,
  effort_ml = effort_ml_fit,
  index_ml = index_ml_fit,
  years = ml_driors$years
)


effort_ml_fit$results %>% 
  filter(variable == "depletion") %>% 
  ggplot(aes(1:length(mean), mean)) + 
  geom_line() + 
  geom_point(data = pop, aes(year, depletion))

# effort_u = tidybayes::gather_draws(effort_bayes_fit$fit,log_f_t[year]) %>% 
#   mutate(u = 1 / (1 + exp(-.value)))

# index_u = tidybayes::gather_draws(index_bayes_fit$fit,log_f_t[year]) %>% 
#   mutate(u = 1 / (1 + exp(-.value)))

# index_catch = tidybayes::gather_draws(index_bayes_fit$fit,catch_hat_t[year])

# ggplot() +
#   geom_violin(data = effort_u, aes(year, u, group = year, fill = "cpue")) +
#   geom_violin(data = index_u, aes(year, u, group = year, fill = "index")) +
#   geom_point(data = pop, aes(year, ((effort * q) / (effort * q + 0.2)) * (1 - exp(-(effort * q + 0.2)))))

index_ml_fit$results %>% 
  filter(variable == "depletion") %>% 
  ggplot(aes(1:length(mean), mean)) + 
  geom_line() + 
  geom_point(data = pop, aes(year, depletion))