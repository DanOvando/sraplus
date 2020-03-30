## ---- echo = FALSE--------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warnings = FALSE,
  eval = TRUE
  # dev = "svg"
)


## ---- eval = FALSE--------------------------------------------------------------
## install.packages("devtools")


## ---- eval = FALSE--------------------------------------------------------------
## remotes::install_github("danovando/sraplus")


## ---- eval = FALSE--------------------------------------------------------------
## 
## install.packages("devtools")
## 
## remotes::install_github("danovando/sraplus")
## 
## 


## -------------------------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(sraplus)



## ----c-msy-1--------------------------------------------------------------------

example_taxa <- "gadus morhua"

data(cod)

head(cod)



## ----c-msy-2--------------------------------------------------------------------


catch_only_driors <- sraplus::format_driors(
  taxa = example_taxa,
  catch = cod$catch,
  years = cod$year,
  use_heuristics = TRUE
)




## ----c-msy-3--------------------------------------------------------------------

plot_driors(catch_only_driors)



## -------------------------------------------------------------------------------

 catch_only_fit <- fit_sraplus(driors = catch_only_driors,
                       engine = "sir",
                       draws = 1e5,
                       n_keep = 2000,
                       estimate_proc_error = FALSE, 
                       estimate_shape = TRUE)
# 
#  catch_only_fit <- fit_sraplus(driors = catch_only_driors,
#                        engine = "sir",
#                        draws = 1e4,
#                        n_keep = 2000,
#                        estimate_proc_error = FALSE, 
#                        estimate_shape = FALSE,
#                        estimate_k = FALSE,
#                        learn_rate = .05)


## -------------------------------------------------------------------------------

head(catch_only_fit$results)



## -------------------------------------------------------------------------------
head(catch_only_fit$fit)


## -------------------------------------------------------------------------------

sraplus::plot_sraplus(catch_only = catch_only_fit, years = catch_only_driors$years)



## ----fmi-sar-1------------------------------------------------------------------

fmi_sar_driors <- format_driors(
  taxa = example_taxa,
  catch = cod$catch,
  years = cod$year,
  initial_state = 1,
  initial_state_cv = 0.25,
  use_heuristics = FALSE,
  sar = 10,
  fmi = c("research" = 0.5, "management" = 0.5, "socioeconomics" = 0.5, 'enforcement' = 0.5),
  sar_cv = NA,
  use_b_reg = FALSE,
  b_ref_type = "k")

sraplus::plot_driors(fmi_sar_driors)


## ----fmi-sar-2------------------------------------------------------------------
fmi_sar_fit <- fit_sraplus(
  driors = fmi_sar_driors,
  engine = "sir",
  draws = 1e6,
  n_keep = 2000,
  estimate_shape = FALSE,
  estimate_proc_error = TRUE
)

plot_sraplus(fmi_sar = fmi_sar_fit,
             catch_only = catch_only_fit,
             years = fmi_sar_driors$years)





## -------------------------------------------------------------------------------

sraplus::diagnose_sraplus(fit = fmi_sar_fit, driors = fmi_sar_driors )



## ----sim-index-1----------------------------------------------------------------

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
  select(year, depletion,catch, effort,u) %>% 
  gather(metric, value, -year) %>% 
  ggplot(aes(year, value)) + 
  geom_point() + 
  facet_wrap(~metric, scales = "free_y") + 
  labs(y = "Value", x = "Year") + 
  sraplus::theme_sraplus()



## -------------------------------------------------------------------------------
index_driors <- format_driors(
  catch = sim$pop$catch,
  years = sim$pop$year,
  index = sim$pop$biomass * 1e-3,
  index_years = sim$pop$year,
  initial_state = 1,
  initial_state_cv = 0.05,
  growth_rate_prior = 0.4,
  growth_rate_prior_cv = 0.1,
  shape_prior = 1.01,
  shape_prior_cv = 0.1)

plot_driors(index_driors)

index_fit <- fit_sraplus(driors = index_driors,
                      engine = "stan",
                      model = "sraplus_tmb", 
                      estimate_proc_error = FALSE,
                      newtonsteps = 10)

plot_sraplus(index = index_fit,years = index_driors$years)

plot_prior_posterior(index_fit, index_driors)




## ----cpue-fit-1-----------------------------------------------------------------

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
    init_u_umsy = 0.75
  )

sim$pop %>% 
  select(year, depletion,catch, effort,u) %>% 
  gather(metric, value, -year) %>% 
  ggplot(aes(year, value)) + 
  geom_point() + 
  facet_wrap(~metric, scales = "free_y") + 
  labs(y = "Value", x = "Year") + 
  sraplus::theme_sraplus()



## ---- results="hide", message=FALSE, warning=FALSE------------------------------


cpue_driors <- format_driors(taxa = example_taxa,
                           catch = sim$pop$catch,
                           years = sim$pop$year,
                           effort = sim$pop$effort,
                           effort_years = sim$pop$year,
                           growth_rate_prior = 0.4,
                           growth_rate_prior_cv = 0.1,
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
                           growth_rate_prior = 0.4,
                           growth_rate_prior_cv = 0.1,
                           shape_prior = 1.01,
                           q_slope_prior = 0.025,
                           q_slope_prior_cv = 0.05)

cpue_qslope_fit <- fit_sraplus(driors = cpue_qslope_driors,
                             engine = "stan",
                             model = "sraplus_tmb",
                             adapt_delta = 0.9,
                             max_treedepth = 10,
                             n_keep = 1000,
                             chains = 1, 
                             cores = 1,
                             estimate_qslope = TRUE,
                             estimate_proc_error = FALSE)



## -------------------------------------------------------------------------------

plot_sraplus(`CPUE fit no qslope` = cpue_fit, `CPUE fit with qslope` =  cpue_qslope_fit, years = cpue_driors$years)



## ---- results = "hide", message = "FALSE", warning = "FALSE"--------------------
cpue_sar_qslope_driors <- format_driors(taxa = example_taxa,
                           catch = sim$pop$catch,
                           years = sim$pop$year,
                           effort = sim$pop$effort,
                           effort_years = sim$pop$year,
                           growth_rate_prior = 0.4,
                           growth_rate_prior_cv = 0.1,
                           shape_prior = 1.01,
                           q_slope_prior = 0.025,
                           q_slope_prior_cv = 0.25,
                           f_ref_type = "f",
                           sar = 2,
                           sar_cv = 0.1)

cpue_sar_qslope_fit <- fit_sraplus(driors = cpue_sar_qslope_driors,
                             engine = "tmb",
                             model = "sraplus_tmb",
                             adapt_delta = 0.9,
                             max_treedepth = 10,
                             n_keep = 2000,
                             chains = 1, 
                             cores = 1,
                             estimate_qslope = TRUE,
                             estimate_proc_error = FALSE)

plot_sraplus(cpue_sar_qslope_fit, years = cpue_sar_qslope_driors$years)


## ---- results = "hide"----------------------------------------------------------

cpue_sar_proc_driors <- format_driors(taxa = example_taxa,
                           catch = sim$pop$catch,
                           years = sim$pop$year,
                           effort = sim$pop$effort,
                           effort_years = sim$pop$year,
                           growth_rate_prior = 0.4,
                           growth_rate_prior_cv = 0.1,
                           shape_prior = 1.01,
                           q_slope_prior = 0,
                           q_slope_prior_cv = 0.25,
                           sar = 4,
                           sar_cv = .05,
                           f_ref_type = "f")

cpue_sar_proc_fit <- fit_sraplus(driors = cpue_sar_proc_driors,
                             engine = "tmb",
                             model = "sraplus_tmb",
                             adapt_delta = 0.9,
                             max_treedepth = 10,
                             n_keep = 2000,
                             chains = 1, 
                             cores = 1,
                             estimate_qslope = FALSE,
                             estimate_proc_error = TRUE)


## -------------------------------------------------------------------------------

plot_sraplus(`no rocess error and no qslope ` = cpue_fit, 
             `no process error with qslope` =  cpue_qslope_fit, 
             `no process error with qslope and sar` = cpue_sar_qslope_fit,
             `process error and sar` = cpue_sar_proc_fit,
             years = cpue_driors$years)




## -------------------------------------------------------------------------------

plot_prior_posterior(cpue_sar_proc_fit, cpue_sar_proc_driors)



## -------------------------------------------------------------------------------
summarize_sralpus(cpue_sar_proc_fit, output = "plot")


## -------------------------------------------------------------------------------
summarize_sralpus(cpue_sar_proc_fit, output = "table") %>% 
  knitr::kable()

