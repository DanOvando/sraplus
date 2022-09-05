
#first, restart R, then
library(remotes)

# remotes::install_github("danovando/sraplus@better-install")


# run test ----------------------------------------------------------------


library(sraplus)

example_taxa <- "gadus morhua"

set.seed(42)

sim <-
  sraplus_simulator(
    sigma_proc = 0.05,
    sigma_u = 0.05,
    q_slope = 0.025,
    r = 0.2,
    years = 25,
    q = 1e-3,
    m = 1.01,
    init_u_umsy = 0.75)

cpue_driors <- format_driors(taxa = example_taxa,
                             catch = sim$pop$catch,
                             years = sim$pop$year,
                             effort = sim$pop$effort,
                             effort_years = sim$pop$year,
                             shape_prior = 1.01,
                             q_slope_prior = 0,
                             q_slope_prior_cv = 0.2)


cpue_fit <- fit_sraplus(driors = cpue_driors,
                        engine = "tmb",
                        model = "sraplus_tmb",
                        adapt_delta = 0.9,
                        max_treedepth = 10,
                        n_keep = 2000,
                        chains = 1, 
                        cores = 1,
                        estimate_qslope = FALSE)


cpue_qslope_driors <- format_driors(taxa = example_taxa,
                                    catch = sim$pop$catch,
                                    years = sim$pop$year,
                                    effort = sim$pop$effort,
                                    effort_years = sim$pop$year,
                                    shape_prior = 1.01,
                                    q_slope_prior = 0.025,
                                    q_slope_prior_cv = 0.05)

cpue_qslope_fit <- fit_sraplus(driors = cpue_qslope_driors,
                               engine = "stan",
                               model = "sraplus_tmb",
                               adapt_delta = 0.9,
                               max_treedepth = 10,
                               n_keep = 1000,
                               chains = 4, 
                               cores = 4,
                               estimate_qslope = TRUE,
                               estimate_proc_error = FALSE)

plot_prior_posterior(cpue_fit, cpue_driors)


plot_sraplus(`CPUE fit no qslope` = cpue_fit, `CPUE fit with qslope` =  cpue_qslope_fit, years = cpue_driors$years)
