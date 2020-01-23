library(tidyverse)

library(sraplus)

library(tmbstan)

rstan_options(auto_write = TRUE)

set.seed(42)

sigma_obs = .1

sigma_proc_ratio = 1

q = 0.0067

sim <-
  sraplus_simulator(
    sigma_proc = sigma_obs * sigma_proc_ratio,
    sigma_u = 0,
    q_slope = 0,
    r = 0.2,
    years = 20,
    q = q,
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


# driors <- format_driors(
#   taxa =
#     example_taxa,
#   catch = pop$catch,
#   years = pop$year,
#   initial_state = pop$depletion[1],
#   initial_state_cv = 0.05,
#   terminal_state = dplyr::last(pop$depletion),
#   terminal_state_cv = 0.01,
#   growth_rate_prior = 0.4,
#   growth_rate_prior_cv = 0.1,
#   use_heuristics = FALSE,
#   shape_prior = 2
# )
# 
# 
# plot_driors(driors)
# 
# sir_fit <- fit_sraplus(driors = driors,
#                        engine = 'sir',
#                        draws = 1e5,
#                        estimate_k = FALSE)
# 
# sir_diagnostics <- diagnose_sraplus(sir_fit, driors)

ml_driors <- format_driors(taxa = example_taxa,
                        catch = pop$catch,
                        years = pop$year,
                        index = pop$biomass * q * exp(rnorm(length(pop$biomass),-sigma_obs^2/2,sigma_obs)),
                        index_years = pop$year,
                        initial_state = 1,
                        initial_state_cv = 0.025,
                        terminal_state = NA,
                        shape_prior = 1.01,
                        growth_rate_prior = 0.4,
                        growth_rate_prior_cv = 0.5,
                        sigma_ratio_prior = 1,
                        sigma_ratio_prior_cv = .1,
                        sigma_obs_prior = 0.1,
                        sigma_obs_prior_cv = 1
                        )

plot_driors(ml_driors)


ml_fit <- fit_sraplus(driors = ml_driors,
                      engine = "tmb",
                      model = "sraplus_tmb",
                      estimate_shape = FALSE, 
                      estimate_proc_error = TRUE,
                      estimate_k = TRUE,
                      learn_rate = 2e-1,
                      n_keep = 2000,
                      eps = 1e-12,
                      adapt_delta = 0.95,
                      max_treedepth = 12)

diagnose_sraplus(ml_fit, ml_driors)

plot_sraplus(ml_fit = ml_fit, years = ml_driors$years)

plot_prior_posterior(ml_fit, ml_driors)

sraplus::summarize_sralpus(ml_fit)


bayes_fit <- fit_sraplus(driors = ml_driors,
                      engine = "stan",
                      model = "sraplus_tmb",
                      estimate_shape = FALSE, 
                      estimate_proc_error = TRUE,
                      estimate_k = TRUE,
                      learn_rate = 2e-1,
                      n_keep = 2000,
                      eps = 1e-3,
                      adapt_delta = 0.95,
                      marginalize_q = FALSE,
                      max_treedepth = 12)

diagnose_sraplus(bayes_fit, ml_driors)

plot_prior_posterior(bayes_fit, ml_driors)


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

sigma_obs_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_sigma_obs",transformations = "exp") +
  geom_vline(aes(xintercept = sigma_obs), color = "red")

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

log_q_guess <- log(mean(ml_driors$index / ml_driors$catch))


log_q_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_q", freq = FALSE) + 
  geom_density(data = data_frame(log_q = rnorm(1000,log_q_guess,.2)), aes(log_q), fill = "lightgrey",
               alpha = 0.5)



# test effort -------------------------------------------------------------

set.seed(42)

sigma_obs <- 0.05

sigma_proc_ratio <- 1

q = 0.006

sim <-
  sraplus_simulator(
    sigma_proc = sigma_obs * sigma_proc_ratio,
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



effort_driors <- format_driors(taxa = example_taxa,
                           catch = pop$catch,
                           years = pop$year,
                           effort = pop$effort * exp(rnorm(length(pop$effort), -sigma_obs^2/2, sigma_obs)),
                           index_years = pop$year,
                           initial_state = 1,
                           initial_state_cv = 0.05,
                           terminal_state = NA,
                           growth_rate_prior = 0.4,
                           growth_rate_prior_cv = 0.1,
                           q_slope_prior = 0.025,
                           q_slope_prior_cv = 0.01)

plot_driors(effort_driors)

index_driors <- format_driors(taxa = example_taxa,
                                  catch = pop$catch,
                                  years = pop$year,
                                  index = pop$biomass * q * exp(rnorm(length(pop$effort), -sigma_obs^2/2, sigma_obs)),
                                  index_years = pop$year,
                                  initial_state = 1,
                                  initial_state_cv = 0.05,
                                  terminal_state = NA,
                                  growth_rate_prior = 0.4,
                                  growth_rate_prior_cv = 0.1,
                                 sigma_ratio_prior = 1,
                                 sigma_ratio_prior_cv = 0.5)

plot_driors(index_driors)


index_bayes_fit <- fit_sraplus(driors = index_driors,
                         engine = "stan",
                         model = "sraplus_tmb",
                         estimate_shape = FALSE, 
                         estimate_proc_error = TRUE,
                         estimate_k = TRUE,
                         learn_rate = 2e-1,
                         n_keep = 2000,
                         eps = 1e-3,
                         adapt_delta = 0.95,
                         marginalize_q = FALSE,
                         max_treedepth = 12)


index_ml_fit <- fit_sraplus(driors = index_driors,
                               engine = "tmb",
                               model = "sraplus_tmb",
                               estimate_shape = FALSE, 
                               estimate_proc_error = TRUE,
                               estimate_k = TRUE,
                               learn_rate = 2e-1,
                               n_keep = 2000,
                               eps = 1e-3,
                               adapt_delta = 0.95,
                               marginalize_q = FALSE,
                               max_treedepth = 12)

plot_sraplus(
  index_ml = index_ml_fit,
  index_bayes = index_bayes_fit,
  years = index_driors$years
)

effort_bayes_fit <- fit_sraplus(driors = effort_driors,
                               engine = "stan",
                               model = "sraplus_tmb",
                               estimate_shape = FALSE, 
                               estimate_proc_error = TRUE,
                               estimate_k = TRUE,
                               learn_rate = 2e-1,
                               n_keep = 2000,
                               eps = 1e-3,
                               adapt_delta = 0.95,
                               marginalize_q = FALSE,
                               max_treedepth = 12)


effort_ml_fit <- fit_sraplus(driors = effort_driors,
                            engine = "tmb",
                            model = "sraplus_tmb",
                            estimate_shape = FALSE, 
                            estimate_proc_error = TRUE,
                            estimate_k = TRUE,
                            learn_rate = 2e-1,
                            n_keep = 2000,
                            eps = 1e-3,
                            adapt_delta = 0.95,
                            marginalize_q = FALSE,
                            max_treedepth = 12)

plot_prior_posterior(effort_ml_fit, effort_driors)

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
  years = effort_driors$years
)


plot_sraplus(
  index_bayes = index_bayes_fit,
  effort_bayes =  effort_bayes_fit,
  effort_ml = effort_ml_fit,
  index_ml = index_ml_fit,
  years = effort_driors$years
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