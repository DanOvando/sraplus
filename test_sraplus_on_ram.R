
library(tidyverse)
library(here)

if (!file.exists("ram.zip")){
download.file("https://www.dropbox.com/s/jpgz0a5s5of3qev/RAM%20v4.491%20Files%20(1-14-20).zip?dl=1", destfile = "ram.zip")

unzip("ram.zip")

}

min_years_catch <- 20

crazy_b = 5

crazy_u = 15

q = 1e-3

ram_dirs <- list.files()

ram_dirs <- ram_dirs[str_detect(ram_dirs,"RAM")]

ram_files <- list.files(ram_dirs, recursive = TRUE)

ram_files <- ram_files[str_detect(ram_files,".RData")]

load(here(ram_dirs,ram_files[2]))

stock <- stock %>%
  left_join(area, by = "areaid")
# catches
ram_catches <- tcbest.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_data_frame() %>%
  gather(stockid, catch,-year)

# B/Bmsy
ram_b_v_bmsy <- divbpref.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_data_frame() %>%
  gather(stockid, b_v_bmsy,-year)

# U/Umsy
ram_u_v_umsy <- divupref.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_data_frame() %>%
  gather(stockid, u_v_umsy,-year)

# Effort
ram_effort <- effort.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_data_frame() %>%
  gather(stockid, effort,-year)

# biomass


ram_total_biomass <- tbbest.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_data_frame() %>%
  gather(stockid, total_biomass,-year)

# ssb

ram_ss_biomass <- ssb.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_data_frame() %>%
  gather(stockid, ss_biomass,-year)


ram_exp_rate <- ram_catches %>%
  left_join(ram_total_biomass, by = c("stockid", "year")) %>%
  mutate(exploitation_rate = catch / total_biomass) %>%
  select(-catch, -total_biomass)

# ram_exp_rate <- erbest.data %>%
#   mutate(year = rownames(.) %>% as.integer()) %>%
#   as_data_frame() %>%
#   gather(stockid, exploitation_rate, -year)

# put it together

ram_data <- ram_catches %>%
  left_join(bioparams_values_views, by = "stockid") %>%
  left_join(ram_b_v_bmsy, by = c("stockid", "year")) %>%
  left_join(ram_u_v_umsy, by = c("stockid", "year")) %>%
  left_join(ram_exp_rate, by = c("stockid", "year")) %>%
  left_join(ram_effort, by = c("stockid", "year")) %>%
  left_join(ram_total_biomass, by = c("stockid", "year")) %>%
  left_join(ram_ss_biomass, by = c("stockid", "year")) %>%
  left_join(stock, by = "stockid") %>%
  select(stockid, scientificname, commonname, everything())


# create new variables

ram_data <- ram_data %>%
  mutate(tb_v_tb0 = total_biomass / TB0,
         ssb_v_ssb0 = ss_biomass / SSB0)

# filter data -------------------------------------------------------------


# for now, only include continuous catch series

ram_data <- ram_data %>%
  filter(is.na(catch) == FALSE) %>%
  group_by(stockid) %>%
  mutate(delta_year = year - lag(year)) %>%
  mutate(delta_year = case_when(year == min(year) ~ as.integer(1),
                                TRUE ~ delta_year)) %>%
  mutate(missing_gaps = any(delta_year > 1)) %>%
  filter(missing_gaps == FALSE) %>%
  mutate(n_years = n_distinct(year)) %>%
  filter(n_years >= min_years_catch) %>%
  filter(all(b_v_bmsy < crazy_b, na.rm = TRUE),
         all(u_v_umsy < crazy_u, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(stockid) %>%
  mutate(
    has_tb0 = !all(is.na(TB0)),
    has_tb = !all(is.na(total_biomass)),
    first_catch_year = year[which(catch > 0)[1]]
  ) %>%
  filter(year >= first_catch_year) %>%
  mutate(
    pchange_effort = lead(u_v_umsy) / (u_v_umsy + 1e-6),
    cs_effort = (u_v_umsy - mean(u_v_umsy)) / sd(u_v_umsy),
    index = total_biomass * q,
    approx_cpue = catch / (u_v_umsy / q + 1e-3),
    b_rel = dplyr::case_when(
      has_tb0 ~ total_biomass / max(TB0),
      has_tb ~ total_biomass / max(total_biomass),
      TRUE ~ b_v_bmsy / 2
    )
  ) %>%
  mutate(approx_cpue = pmin(quantile(approx_cpue, 0.9, na.rm = TRUE), approx_cpue)) %>%
  ungroup()

ram_b_plot <- ram_data %>%
  ggplot(aes(x = year, y = b_v_bmsy)) +
  geom_bin2d() +
  scale_fill_viridis_c()

kobe_panels <- ram_data %>%
  filter(year >= 1950) %>%
  mutate(year_block = plyr::round_any(year, 10, floor)) %>%
  ggplot(aes(x = b_v_bmsy, y = u_v_umsy)) +
  geom_bin2d(binwidth = c(0.5, 0.5)) +
  facet_wrap(~ year_block) +
  scale_fill_viridis_c()


# try fitting ram stocks


ram_fits <- ram_data %>% 
  group_by(stockid) %>% 
  mutate(has_index = any(!is.na(index))) %>% 
  filter(has_index) %>% 
  nest()


dat <- ram_fits$data[[100]]

dat <- dat %>% filter(!is.na(catch))

sigma_obs <- 0.01

index_years <- dat$year[!is.na(dat$index)]

driors <- format_driors(taxa = dat$scientificname[1],
                           catch = dat$catch,
                           years = dat$year,
                           index =  dat$index[dat$year %in% index_years]* exp(rnorm(length(index_years),-sigma_obs^2/2,sigma_obs)),
                           index_years =index_years,
                           initial_state = 1,
                           initial_state_cv = .1,
                           terminal_state = NA,
                           shape_prior = 1.01,
                           shape_prior_cv = 1,
                           growth_rate_prior = 0.4,
                           growth_rate_prior_cv = 0.5,
                           sigma_ratio_prior = 1,
                           sigma_ratio_prior_cv = .1,
                           sigma_obs_prior = 0.1,
                           sigma_obs_prior_cv = 1
)

plot_driors(driors)


fit <- fit_sraplus(driors = driors,
                      engine = "stan",
                      model = "sraplus_tmb",
                      estimate_shape = FALSE, 
                      estimate_proc_error = TRUE,
                      estimate_k = TRUE,
                      learn_rate = 2e-1,
                      n_keep = 2000,
                      eps = 1e-12,
                      adapt_delta = 0.95,
                      marginalize_q = FALSE,
                      max_treedepth = 12)

diagnose_sraplus(fit, driors)

plot_sraplus(ram_fit = fit, years = driors$years)

plot_prior_posterior(fit, driors)

fit_b <- fit$results %>% 
  filter(variable == "b_div_bmsy") %>% 
  left_join(dat %>% select(year, b_v_bmsy), by = "year" )

fit_plot <- fit_b %>% 
  ggplot() + 
  geom_point(aes(year, b_v_bmsy, color = "RAM"), size = 2) + 
  geom_line(aes(year, mean, color = "sraplus"))

fit_plot