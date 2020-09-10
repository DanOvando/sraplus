library(tidyverse)

library(sraplus)

# library(tmbstan)

Sys.unsetenv("PKG_CXXFLAGS")

example_taxa <- "gadus morhua"


set.seed(42)

sigma_obs = 0.025

sigma_proc_ratio = 0.1

q = 1

sim <-
  sraplus_simulator(
    sigma_proc = sigma_obs * sigma_proc_ratio,
    sigma_u = 0,
    q_slope = 0,
    r = 0.2,
    years = 25,
    q = q,
    m = 1.1,
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



com_driors <-
  format_driors(
    catch = pop$catch,
    years = pop$year,
    use_heuristics = FALSE,
    initial_state = NA,
    initial_state_cv = .2,
    terminal_state = NA,
    terminal_state_cv = 0.1,
    b_ref_type = "k",
    use_catch_prior = TRUE,
    taxa = example_taxa,
    shape_prior_source = "thorson"
  )

plot_driors(com_driors)

com_fit <-
  fit_sraplus(
    driors = com_driors,
    include_fit = TRUE,
    engine = "sir",
    draws = 1e6,
    tune_prior_predictive = TRUE
  )

plot_prior_posterior(com_fit, com_driors)

ml_driors <- format_driors(taxa = example_taxa,
                           catch = pop$catch,
                           years = pop$year,
                           index = pop$biomass * q * exp(rnorm(length(pop$biomass),-sigma_obs^2/2,sigma_obs)),
                           index_years = pop$year,
                           initial_state = NA,
                           initial_state_cv = .01,
                           terminal_state = NA,
                           terminal_state_cv = 0.1,
                           shape_prior = 1.1,
                           growth_rate_prior = NA,
                           growth_rate_prior_cv = 0.5,
                           sigma_ratio_prior = 1,
                           sigma_ratio_prior_cv = .1,
                           q_prior = 1,
                           q_prior_cv = 0.001
)

plot_driors(ml_driors)


ml_fit <- fit_sraplus(driors = ml_driors,
                      engine = "tmb",
                      model = "sraplus_tmb",
                      estimate_proc_error = TRUE,
                      estimate_initial_state = TRUE,
                      estimate_q = FALSE,
                      learn_rate = 2e-1,
                      n_keep = 2000,
                      eps = 1e-12,
                      adapt_delta = 0.95,
                      analytical_q = FALSE,
                      max_treedepth = 12,
                      ci  = 0.89,
                      tune_prior_predictive = TRUE)

diagnose_sraplus(ml_fit, ml_driors)

plot_sraplus(ml_fit = ml_fit, years = ml_driors$years)

plot_prior_posterior(ml_fit, ml_driors)

sraplus::summarize_sralpus(ml_fit)

a <- Sys.time()
bayes_fit <- fit_sraplus(driors = ml_driors,
                      engine = "stan",
                      model = "sraplus_tmb",
                      estimate_proc_error = TRUE,
                      estimate_f = FALSE,
                      estimate_k = TRUE,
                      learn_rate = 2e-1,
                      n_keep = 2000,
                      eps = 1e-2,
                      adapt_delta = 0.8,
                      estimate_q = FALSE,
                      max_treedepth = 12,
                      refresh = 250,
                      estimate_initial_state = FALSE,
                      thin_draws = TRUE,
                      workers = 8)
Sys.time() - a

diagnose_sraplus(bayes_fit, ml_driors)
  
plot_prior_posterior(bayes_fit, ml_driors)

a = get_prior_posterior(bayes_fit,ml_driors)

sraplus::summarize_sralpus(bayes_fit)

test <- names(bayes_fit$fit)

# pair_vars <- test[!str_detect(test,"uc_proc")]
# 
# pairs(bayes_fit$fit, pars = pair_vars)


plot_sraplus(ml_fit = ml_fit, bayes_fit = bayes_fit,years = ml_driors$years)


# test u priors

u_driors <- format_driors(
  taxa = example_taxa,
  catch = pop$catch,
  years = pop$year,
  index = pop$biomass * q * exp(rnorm(
    length(pop$biomass), -sigma_obs ^ 2 / 2, sigma_obs
  )),
  index_years = pop$year,
  u = pop$u,
  u_years = pop$year,
  u_cv = .1,
  f_ref_type = "f",
  initial_state = 1,
  initial_state_cv = 0.025,
  terminal_state = NA,
  shape_prior = 1.01,
  growth_rate_prior = 0.4,
  growth_rate_prior_cv = 0.5,
  sigma_ratio_prior = 1,
  sigma_ratio_prior_cv = .1,
)

plot_driors(u_driors)


u_fit <- fit_sraplus(driors = u_driors,
                      engine = "tmb",
                      model = "sraplus_tmb",
                      estimate_shape = FALSE, 
                      estimate_proc_error = FALSE,
                      estimate_f = FALSE,
                      estimate_k = TRUE,
                      learn_rate = 2e-1,
                      n_keep = 2000,
                      eps = 1e-12,
                      adapt_delta = 0.95,
                      analytical_q = FALSE,
                      max_treedepth = 12,
                      ci  = 0.89)

diagnose_sraplus(u_fit, u_driors)

plot_sraplus(u_fit = u_fit, years = u_fit$years)

plot_prior_posterior(u_fit, u_driors)

sraplus::summarize_sralpus(ml_fit)


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

r_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_anchor",transformations = "identity") + 
  geom_vline(aes(xintercept = sim$params$r), color = "red")

r_hat <- rstan::extract(bayes_fit$fit, "log_anchor")  


# m_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_m",transformations = "exp") +
#   geom_vline(aes(xintercept = sim$params$m), color = "red")

# q_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_q",transformations = "exp") + 
#   geom_vline(aes(xintercept = sim$pop$q[1]), color = "red")

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


# log_q_hat <- bayesplot::mcmc_hist(as.matrix(bayes_fit$fit), "log_q", freq = FALSE) + 
#   geom_density(data = data_frame(log_q = rnorm(1000,log_q_guess,.2)), aes(log_q), fill = "lightgrey",
#                alpha = 0.5)



# test effort -------------------------------------------------------------

set.seed(24)

sigma_obs <- 0.1

sigma_proc_ratio <- 1

q = 0.006

sim <-
  sraplus_simulator(
    sigma_proc = sigma_obs * sigma_proc_ratio,
    sigma_u = 0.2,
    q_slope = 0,
    r = 0.4,
    years = 25,
    q = 1e-3,
    m = 1.01,
    init_u_umsy = .2
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


# u_driors <-  format_driors(taxa = example_taxa,
#                            catch = pop$catch,
#                            years = pop$year,
#                            index = pop$biomass * q * exp(rnorm(length(pop$effort), -sigma_obs^2/2, sigma_obs)),
#                            index_years = pop$year,
#                            # u = pop$u,
#                            # u_years = pop$year,
#                            u_cv = 0.2,
#                            f_ref_type = "fmsy",
#                            f_prior_form = 1,
#                            initial_state = 1,
#                            initial_state_cv = 0.05,
#                            terminal_state = NA,
#                            growth_rate_prior = 0.4,
#                            growth_rate_prior_cv = 0.1,
#                            q_slope_prior = 0.025,
#                            q_slope_prior_cv = 0.01)
# 
# 
# plot_driors(u_driors)
# 
# 
# u_fit <- fit_sraplus(driors = u_driors,
#                                engine = "stan",
#                                model = "sraplus_tmb",
#                                estimate_shape = FALSE, 
#                                estimate_proc_error = TRUE,
#                               estimate_f = TRUE,
#                                n_keep = 2000,
#                                adapt_delta = 0.95,
#                                analytical_q = FALSE,
#                                max_treedepth = 12)
# 
# u_fit_false <- fit_sraplus(driors = u_driors,
#                      engine = "tmb",
#                      model = "sraplus_tmb",
#                      estimate_shape = FALSE, 
#                      estimate_proc_error = FALSE,
#                      estimate_f = FALSE,
#                      n_keep = 2000,
#                      adapt_delta = 0.95,
#                      analytical_q = FALSE,
#                      max_treedepth = 12)
# 
# 
# plot_sraplus(fit_u = u_fit, dont_fit_u = u_fit_false)
# 
# plot_prior_posterior(u_fit, u_driors)



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
                           q_slope_prior = 0.25,
                           q_slope_prior_cv = 0.01)

plot_driors(effort_driors)



index_driors <- format_driors(taxa = example_taxa,
                                  catch = pop$catch,
                                  years = pop$year,
                                  index = pop$biomass * q * exp(rnorm(length(pop$effort), -sigma_obs^2/2, sigma_obs)),
                                  index_years = pop$year,
                                  initial_state = 1,
                                  initial_state_cv = 0.05,
                                  terminal_state = NA)

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
                         analytical_q = FALSE,
                         max_treedepth = 12,
                         thin_draws = TRUE)

plot_sraplus(index_bayes_fit)

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
                               analytical_q = FALSE,
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
                               analytical_q = FALSE,
                               max_treedepth = 12,
                               thin_draws = TRUE)


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
                            analytical_q = FALSE,
                            max_treedepth = 12)

plot_prior_posterior(effort_ml_fit, effort_driors)

# effort_q_hat <- bayesplot::mcmc_hist(as.matrix(effort_bayes_fit$fit), "log_q",transformations = "exp") + 
#   geom_vline(aes(xintercept = sim$pop$q[1]), color = "red")


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