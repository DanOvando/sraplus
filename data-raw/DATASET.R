## code to prepare prior regressions for sraplus

# fit prior relationships for sraplus


# setup -------------------------------------------------------------------

library(tidyverse)
library(rstan)
library(rstanarm)
library(patchwork)
library(scales)
library(hrbrthemes)
library(here)
library(recipes)
options(mc.cores = parallel::detectCores() / 2)
rstan_options(auto_write = TRUE)
theme_set(theme_classic() + theme(strip.background = element_rect(color = "transparent")))

# options -----------------------------------------------------------------

min_years_catch <- 20

crazy_b <- 4 # threshold for suspect B/Bmsy value

crazy_u <- 5 # threshold for suspect U/Umsy value

draws <- 3000

min_draws <- 2000 # minimum number of unique SIR draws

n_cores <- 6 # number of cores for parallel processing

lookup_fmi_names <- FALSE

future::plan(future::multiprocess, workers = n_cores)

# data(Return)

# return is from here, WARNING, changes rapidly, things break check and make sure this isn't why
# https://drive.google.com/drive/u/0/folders/1J46tM6PYDdPwhx5zGrlHMdxUyGRrky7X?ogsrc=32

functions <- list.files(here::here("R"))

functions <- functions[!functions %in% c("zzz.R", "sysdata.rda")]

purrr::walk(functions, ~ source(here::here("R", .x)))

# load data ---------------------------------------------------------------

load(here("data-raw","Return.Rdata"))

FishLifeData<- Return[c("ParentChild_gz","beta_gv","Cov_gvv")]

FishLifeData$metadata <- "emailed from Thorson, beta version with newparameters"

usethis::use_data(FishLifeData, overwrite = TRUE, internal = TRUE)


if (file.exists(here("data-raw","ram.RData")) == FALSE) {
  # for now storing in my google drive... will need to put this in a better and public location
  
  ram <-
    googledrive::drive_get(path = "~/Databases/RAM/RAM-v4.41-8_20_18/DB-Files-With-Model-Fit-Data/DBdata-model_fits_included.RData") ## if this is your first time running this will be prompted to sign into your google account
  
  googledrive::drive_download(file = googledrive::as_id(ram$id), path = "data/ram.RData")
}

load(here::here("data-raw", "ram.RData"))



min_years_catch = 20
crazy_b = 4
crazy_u = 5
lookup_fmi_names = FALSE


min_years_catch = 20
crazy_b = 4
crazy_u = 5
lookup_fmi_names = FALSE


# process data ------------------------------------------------------------

# catches
ram_catches <- tcbest.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_data_frame() %>%
  gather(stockid, catch, -year)

# B/Bmsy
ram_b_v_bmsy <- divbpref.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_data_frame() %>%
  gather(stockid, b_v_bmsy, -year)

# U/Umsy
ram_u_v_umsy <- divupref.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_data_frame() %>%
  gather(stockid, u_v_umsy, -year)

# Effort
ram_effort <- effort.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_data_frame() %>%
  gather(stockid, effort, -year)

# biomass


ram_total_biomass <- tbbest.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_data_frame() %>%
  gather(stockid, total_biomass, -year)

# ssb

ram_ss_biomass <- ssb.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_data_frame() %>%
  gather(stockid, ss_biomass, -year)


ram_exp_rate <- ram_catches %>%
  left_join(ram_total_biomass, by = c("stockid","year")) %>%
  mutate(exploitation_rate = catch / total_biomass) %>%
  select(-catch,-total_biomass)

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
  mutate(n_years = n_distinct(year)) %>%
  # filter(all(b_v_bmsy < crazy_b, na.rm = TRUE),
  #        all(u_v_umsy < crazy_u, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(stockid) %>%
  mutate(
    has_tb0 = !all(is.na(TB0)),
    has_tb = !all(is.na(total_biomass)),
    first_catch_year = year[which(catch > 0)[1]]
  ) %>%
  # filter(year >= first_catch_year) %>%
  ungroup()


# load other data ---------------------------------------------------------



# load prices

prices <<- read_csv(here("data-raw", "Exvessel Price Database.csv")) %>%
  janitor::clean_names() %>%
  rename(scientificname = scientific_name) %>%
  mutate(log_exvessel = log(exvessel)) %>%
  group_by(asfis_species, pooled_commodity, group_for_pairing) %>%
  mutate(lag_exvessel = lag(exvessel)) %>%
  ungroup()

ram_data <- ram_data %>%
  left_join(prices, by = c("scientificname", "year"))

# fao and effort data

fao_to_effort <-
  read_csv(here("data-raw", "fao-to-bell-region.csv")) %>%
  rename(bell_region = region)

country_to_fao <-
  read_csv(here("data-raw", "country-to-fao-area.csv")) %>%
  unique() %>%
  janitor::clean_names() %>%
  rename(fao_fishing_area = fishing_area_fao_major_fishing_area_1) %>%
  left_join(fao_to_effort, by = "fao_fishing_area")


fao <-
  read_csv(here::here("data-raw","tidy_fao_capture_1950-2016.csv")) %>%
  mutate(id = paste(country, fao_area, common_name, sep = "-")) %>%
  group_by(id) %>%
  mutate(first_year = year[capture > 0 & !is.na(capture)][1]) %>%
  filter(
    year >= first_year,
    capture_units == "Tonnes"
  ) %>%
  ungroup() %>%
  group_by(id) %>%
  mutate(missing_catch = sum(is.na(capture))) %>%
  # filter(missing_catch == 0) %>%
  ungroup() %>%
  mutate(fao_area_code = as.numeric(fao_area_code)) %>%
  filter(
    !str_detect(country, "Totals"),
    isscaap_number < 60
  )

fao_species <- fao %>%
  select(scientific_name, common_name,isscaap_group, isscaap_number) %>%
  unique()

fao_genus <-
  str_split(fao_species$scientific_name, ' ', simplify = TRUE)[,1]

fao_genus = fao_species %>%
  mutate(genus = fao_genus) %>%
  group_by(genus, isscaap_group) %>%
  count() %>%
  group_by(genus) %>%
  filter(n == max(n)) %>%
  select(-n) %>%
  ungroup()


fao$fao_country_name <-
  countrycode::countrycode(fao$country, "country.name", "fao.name")

fao <- fao %>%
  mutate(country = case_when(is.na(fao_country_name) ~ country, TRUE ~ fao_country_name))


fao$continent <-
  countrycode::countrycode(fao$country, "country.name", "continent")

fao_stock_lookup <- fao %>%
  select(scientific_name, common_name, country, fao_area, fao_area_code) %>%
  unique()

plot_prior_fit <- function(fit, split) {
  fit_r2 <- bayes_R2(fit)
  
  
  ppc_plot <- bayesplot::ppc_scatter_avg(y = split$log_value, yrep = posterior_predict(fit)) +
    labs(
      x = "Mean Posterior Predicted log(value)",
      y = "Observed log(value)"
    ) +
    geom_smooth(method = "lm", se = TRUE)
  
  br2_plot <- ggplot() +
    geom_density(data = data_frame(r2 = fit_r2), aes(r2), fill = "lightgrey") +
    labs(x = bquote(R^2), y = "Count")
  
  ppc_plot + br2_plot + plot_layout(ncol = 2, widths = c(3, 1)) +
    plot_annotation(tag_levels = "A")
  
  
}


# fmi data ----------------------------------------------------------------


ram_fmi_sheets <-
  readxl::excel_sheets(here("data-raw", "RAM FMI stocks and taxonomic groupings.xlsx"))


ram_fmi <-
  map(ram_fmi_sheets, ~ readxl::read_xlsx(
    here("data-raw", "RAM FMI stocks and taxonomic groupings.xlsx"),
    sheet = .x
  )) %>%
  set_names(str_replace(tolower(ram_fmi_sheets), " ", "_"))

ram_species <- ram_fmi$ram_species %>%
  janitor::clean_names() %>%
  rename(scientificname = scientificname_ram)

ram_data <<- ram_data %>%
  left_join(ram_species, by = "scientificname")

# stock$country <- ifelse(stock$country == "USA", "United States of America", stock$country)


ram_vs_fmi <- read_csv(here::here("data-raw", "RAM vs FMI.csv"))

ram_fmi_sheets <-
  readxl::excel_sheets(here("data-raw", "RAM FMI stocks and taxonomic groupings.xlsx"))

ram_fmi <-
  map(ram_fmi_sheets, ~ readxl::read_xlsx(
    here("data-raw", "RAM FMI stocks and taxonomic groupings.xlsx"),
    sheet = .x
  )) %>%
  set_names(str_replace(tolower(ram_fmi_sheets), " ", "_"))



fmi <-
  readxl::read_xlsx(here::here("data-raw", "FMI data extract by stock for Ray 2018-11-05.xlsx"),
                    sheet = "summary"
  ) %>%
  janitor::clean_names()


fmi$fao_country_name <-
  countrycode::countrycode(fmi$country_rfmo, "country.name", "fao.name")

fmi <- fmi %>%
  mutate(country_rfmo = case_when(is.na(fao_country_name) ~ country_rfmo, TRUE ~ fao_country_name))


if (lookup_fmi_names == TRUE) {
  fmi_names <- fmi %>%
    select(lookup_code, country_rfmo)
  
  temp <-
    fmi_names$lookup_code %>%
    str_split("\\|", simplify = TRUE) %>%
    as_data_frame() %>%
    set_names("code", "species", "region") %>%
    map_df(str_trim) %>%
    select(species) %>%
    unique() %>%
    mutate(sciname = map(species, ~ taxize::comm2sci(commnames = .x, db = "worms")))
  
  get_names <- function(name) {
    name <- name[[1]]
    if (length(name) == 0) {
      out <- NA
    } else {
      out <- name[1]
    }
    return(out)
  }
  
  
  fmi_scinames <- map_chr(temp$sciname, get_names)
  
  temp <- temp %>%
    mutate(scientificname = fmi_scinames) %>%
    select(-sciname)
  
  saveRDS(temp, file = here::here("data-raw", "fmi_scinames.rds"))
  fmi_scinames <- temp
} else {
  fmi_scinames <-
    readRDS(file = here::here("data-raw", "fmi_scinames.rds"))
}

fmi <- fmi %>%
  left_join(fmi_scinames, by = "species")

stock <- stock %>%
  left_join(area, by = "areaid")

stock <- stock %>%
  left_join(
    fmi %>% select(lookup_code, country_rfmo, scientificname),
    by = c("scientificname", "country" = "country_rfmo")
  )

fmi <- fmi %>%
  left_join(stock %>% select(stockid, lookup_code), by = "lookup_code")


manual_fmi_to_fao <-
  tribble(
    ~country_rfmo,
    ~basin,
    ~fao_area_code,
    "the United States of America",
    "Pac",
    67,
    "the United States of America",
    "Pac",
    77,
    "the United States of America",
    "Atl",
    21,
    "the United States of America",
    "Atl",
    31,
    "Canada",
    "Atl",
    21,
    "Canada",
    "Pac",
    67,
    "the Russian Federation",
    "Atl",
    27,
    "the Russian Federation",
    "Pac",
    61,
    "Spain",
    "Atl",
    27,
    "Spain",
    "Atl",
    34
  )

likely_country_to_fao <- fao_stock_lookup %>%
  group_by(country, fao_area_code) %>%
  count() %>%
  group_by(country) %>%
  mutate(rank_n = n / sum(n)) %>%
  filter(rank_n > 0.15) # serious hack to try and keep it nearshore, knocks out China in 87 but unfortunately also USA in 67.. but those are mostly identified to basin in FMI so we'll take it

fmi <-
  fmi %>%
  mutate(genus = map_chr(scientificname,~str_split(.x," ", simplify = TRUE)[,1])) %>%
  left_join(fao_species, by = c("scientificname" = "scientific_name")) %>%
  left_join(fao_genus, by = "genus") %>%
  mutate(isscaap_group = ifelse(is.na(isscaap_group.x), isscaap_group.y, isscaap_group.x)) %>%
  select(-isscaap_group.x,-isscaap_group.y) %>%
  left_join(fao_stock_lookup,
            by = c("species" = "common_name", "country_rfmo" = "country")
  ) %>%
  rename(best_fao_code = fao_area_code) %>%
  left_join(manual_fmi_to_fao, by = c("country_rfmo", "basin")) %>%
  rename(second_best_fao_code = fao_area_code) %>%
  left_join(likely_country_to_fao, by = c("country_rfmo" = "country")) %>%
  rename(third_best_fao_code = fao_area_code) %>%
  nest(-lookup_code)

find_best_fao <- function(data) {
  # data <- temp$data[[55]]
  
  best_match <- data %>% filter(!is.na(best_fao_code))
  
  second_best_match <- data %>% filter(!is.na(second_best_fao_code))
  
  third_best_match <- data %>% filter(!is.na(third_best_fao_code))
  
  if (nrow(best_match) > 0) {
    out <- best_match %>%
      mutate(fao_area_code = best_fao_code)
  } else if (nrow(second_best_match) > 0) {
    out <- second_best_match %>%
      mutate(fao_area_code = second_best_fao_code)
  } else if (nrow(third_best_match) > 0) {
    out <- third_best_match %>%
      mutate(fao_area_code = third_best_fao_code)
  } else {
    out <- best_match
  }
  
  out <- out %>%
    select(-(fao_area:rank_n))
  return(out)
}

fmi <- fmi %>%
  mutate(best_matches = map(data, find_best_fao)) %>%
  select(-data) %>%
  unnest()


usethis::use_data(fmi, overwrite = TRUE)


# load sar data -----------------------------------------------------------

recent_ram <- ram_data %>%
  group_by(stockid) %>%
  mutate(c_maxc = catch / max(catch, na.rm = TRUE),
         c_meanc = catch / mean(catch, na.rm = TRUE)) %>%
  filter(year > (max(year[!is.na(u_v_umsy)]) - 5)) %>%
  summarise(mean_bbmsy = mean(b_v_bmsy, na.rm = TRUE),
            mean_uumsy = mean(u_v_umsy, na.rm = TRUE),
            mean_f = mean(-log(1 - pmin(0.95,exploitation_rate)), na.rm = TRUE),
            c_div_max_c = mean(c_maxc),
            c_div_mean_c = mean(c_meanc)) %>%
  na.omit()



sar_coverage <- readr::read_csv(here("data-raw","OverlapTable2.csv")) %>%
  janitor::clean_names() %>% 
  group_by(stockid) %>% 
  summarise(mean_tbp_in_stock = mean(tbp_in_stock),
            mean_stock_in_tbp = mean(stock_in_tbp)) %>% 
  ungroup()
  

sar_to_ram <-
  readr::read_csv(here("data-raw", "RamStocksWithID2.csv")) %>%
  janitor::clean_names() %>%
  map_df(stringr::str_trim) %>% # sigh, white spaces in the numerics work on mac but not linux
  modify_at(4:7, as.numeric) %>% # fun fun
  mutate(log_f = log(fstatus + 1e-3)) %>%
  mutate(genus = map_chr(latin_binomial,~str_split(.x," ", simplify = TRUE)[,1])) %>%
  left_join(fao_species, by = c("latin_binomial" = "scientific_name")) %>%
  left_join(fao_genus, by = "genus") %>%
  mutate(isscaap_group = ifelse(is.na(isscaap_group.x), isscaap_group.y, isscaap_group.x)) %>%
  select(-isscaap_group.x, -isscaap_group.y) %>%
  left_join(sar_coverage, by = "stockid") %>%
  left_join(recent_ram, by = "stockid")

sar = sar_to_ram

usethis::use_data(sar, overwrite = TRUE)



# fit fmi models ----------------------------------------------------------


ram_v_fmi <- ram_data %>%
  group_by(stockid) %>%
  mutate(
    c_maxc = catch / max(catch, na.rm = TRUE),
    c_meanc = catch / mean(catch, na.rm = TRUE)
  ) %>%
  filter(year > (max(year) - 5)) %>%
  summarise(
    mean_bbmsy = mean(b_v_bmsy, na.rm = TRUE),
    mean_uumsy = mean(u_v_umsy, na.rm = TRUE),
    mean_f = mean(-log(1 - pmin(0.95,exploitation_rate)), na.rm = TRUE),
    c_div_max_c = mean(c_maxc),
    c_div_mean_c = mean(c_meanc)
  ) %>%
  na.omit() %>%
  gather(metric, value,-stockid, -c_div_max_c, -c_div_mean_c) %>%
  ungroup() %>%
  left_join(fmi, by = "stockid") %>%
  filter(!is.na(lookup_code)) %>%
  select(-basin, -stock, -scientificname, -scientific_name) %>%
  mutate(log_value = log(value)) %>%
  unique() #%>%

ram_v_fmi <- recipe(log_value ~ ., data = ram_v_fmi) %>%
  step_other(isscaap_group) %>%
  prep(data = ram_v_fmi, retain = TRUE) %>%
  juice()

# filter(!(metric == "mean_u" & value < 0.13 & enforcement > 0.6)) # per discussions with coauthors? don't remember why

# huh <- unique(fmi$stockid[!fmi$stockid %in% ram_v_fmi$stockid])
#
# wtf <- ram_data %>%
#   filter(stockid %in% huh) %>%
#   select(year, stockid, b_v_bmsy, u_v_umsy, catch)


random_fmi_tests <- ram_v_fmi %>%
  nest(-metric) %>%
  mutate(splits = map(data, ~ rsample::vfold_cv(.x, v = 2, repeats = 1))) %>%
  select(-data) %>%
  unnest() %>%
  mutate(sampid  = 1:nrow(.))

model_structures <-
  purrr::cross_df(list(
    sampid = random_fmi_tests$sampid,
    model_structure = c(
      "log_value ~ research + management + enforcement + socioeconomics + c_div_max_c + c_div_mean_c" ,
      "log_value ~ (research + management + enforcement + socioeconomics + c_div_max_c  + c_div_mean_c - 1|isscaap_group)",
      "log_value ~ c_div_max_c + c_div_mean_c  + (research + management + enforcement + socioeconomics - 1|isscaap_group)",
      "log_value ~ research + management + enforcement + socioeconomics",
      "log_value ~ (research + management + enforcement + socioeconomics - 1|isscaap_group)",
      "log_value ~ + management + enforcement + socioeconomics + c_div_max_c + c_div_mean_c"
      
    )
  ))

random_fmi_tests <- model_structures %>%
  left_join(random_fmi_tests, by = "sampid")

random_fmi_tests <- random_fmi_tests %>%
  mutate(
    fmi_fit = map2(
      splits,
      model_structure,
      fit_prior_regressions,
      produce = "summary",
      refresh = 0
    )
  )

random_fmi_tests <- random_fmi_tests %>%
  mutate(
    training_performance = map(fmi_fit, "training_summary"),
    testing_performance = map(fmi_fit, "testing_summary")
  )

random_fmi_tests <- random_fmi_tests %>%
  mutate(
    training_rmse = map_dbl(
      training_performance,
      ~ yardstick::rmse_vec(truth = .x$observed,
                            estimate = .x$pp_pred)
    ),
    testing_rmse = map_dbl(
      testing_performance,
      ~ yardstick::rmse_vec(truth = .x$observed,
                            estimate = .x$pp_pred)
    )
  )


best_fmi_models <- random_fmi_tests %>%
  group_by(metric, model_structure) %>%
  summarise(
    mean_testing_rmse = mean(testing_rmse),
    mean_training_rmse = mean(training_rmse)
  ) %>%
  group_by(metric) %>%
  filter(mean_testing_rmse == min(mean_testing_rmse))

best_fmi_models <- best_fmi_models %>%
  mutate(splits = map(metric, ~ filter(ram_v_fmi, metric == .x))) %>%
  mutate(
    best_fmi_fit = map2(
      splits,
      model_structure,
      fit_prior_regressions,
      produce = "results",
      refresh = 100,
      use_splits = FALSE
    )
  )




best_fmi_models <- best_fmi_models %>%
  mutate(fit = map(best_fmi_fit, "fit")) %>% 
  mutate(prior_plot = map2(fit, splits, plot_prior_fit))


usethis::use_data(best_fmi_models,overwrite = TRUE)


# fit sar models --------------------------------------------------------------

ram_v_sar <- sar_to_ram %>%
  gather(metric, value, contains("mean_"), -c_div_mean_c,-mean_stock_in_tbp,-mean_tbp_in_stock) %>%
  mutate(log_value = log(value + 1e-3)) %>%
  mutate(sar_2 = sar ^ 2) %>%
  select(stockid, sar, sar_2, isscaap_group, metric, value, log_value, c_div_mean_c, c_div_max_c,mean_stock_in_tbp) %>% 
  na.omit() %>%
  filter(mean_stock_in_tbp > 50)

ram_v_sar %>% 
  ggplot(aes(sar, log_value, color = mean_stock_in_tbp)) + 
  geom_point() + 
  facet_wrap(~metric)

ram_v_sar <- recipe(log_value ~ ., data = ram_v_sar) %>%
  step_other(isscaap_group) %>%
  prep(data = sar_data, retain = TRUE) %>%
  juice()


random_sar_tests <- ram_v_sar %>%
  nest(-metric) %>%
  mutate(splits = map(data, ~ rsample::vfold_cv(.x, v = 2, repeats = 1))) %>%
  select(-data) %>%
  unnest() %>%
  mutate(sampid  = 1:nrow(.))

model_structures <-
  purrr::cross_df(list(
    sampid = random_sar_tests$sampid,
    model_structure = c(
      "log_value ~ poly(sar,2) + c_div_max_c + c_div_mean_c",
      "log_value ~ sar + c_div_max_c + c_div_mean_c",
      "log_value ~ c_div_max_c + c_div_mean_c + (sar - 1|isscaap_group)",
      "log_value ~ c_div_max_c + c_div_mean_c + (sar + sar_2 - 1|isscaap_group)",
      "log_value ~  (sar + sar_2 - 1|isscaap_group)"
      
      
    )
  ))

random_sar_tests <- model_structures %>%
  left_join(random_sar_tests, by = "sampid")

random_sar_tests <- random_sar_tests %>%
  mutate(fit = map2(
    splits,
    model_structure,
    fit_prior_regressions,
    produce = "summary",
    refresh = 0
  ))

random_sar_tests <- random_sar_tests %>%
  mutate(
    training_performance = map(fit, "training_summary"),
    testing_performance = map(fit, "testing_summary")
  )

random_sar_tests <- random_sar_tests %>%
  mutate(
    training_rmse = map_dbl(
      training_performance,
      ~ yardstick::rmse_vec(truth = .x$observed,
                            estimate = .x$pp_pred)
    ),
    testing_rmse = map_dbl(
      testing_performance,
      ~ yardstick::rmse_vec(truth = .x$observed,
                            estimate = .x$pp_pred)
    )
  )

# random_sar_tests %>%
#   ggplot(aes(model_structure, testing_rmse, color = metric)) +
#   geom_point() +
#   coord_flip()

best_sar_models <- random_sar_tests %>%
  group_by(metric, model_structure) %>%
  summarise(mean_rmse = mean(testing_rmse)) %>%
  group_by(metric) %>%
  filter(mean_rmse == min(mean_rmse))

best_sar_models <- best_sar_models %>%
  mutate(splits = map(metric, ~ filter(ram_v_sar, metric == .x))) %>%
  mutate(
    best_fit = map2(
      splits,
      model_structure,
      fit_prior_regressions,
      produce = "results",
      refresh = 100,
      use_splits = FALSE
    )
  )


best_sar_models <- best_sar_models %>%
  mutate(fit = map(best_fit, "fit")) %>% 
  mutate(prior_plot = map2(fit, splits, plot_prior_fit))


sar_v_f_plot <- bayesplot::ppc_intervals(
  x = best_sar_models$splits[[3]]$sar,
  y = best_sar_models$splits[[3]]$log_value,
  yrep = posterior_predict(best_sar_models$fit[[3]])
) +
  labs(
    x = "SAR",
    y = "log(U/Umsy)"
  )


usethis::use_data(best_sar_models,overwrite = TRUE)

# catch priors ------------------------------------------------------------



# catch_data <- ram_data %>%
#   group_by(stockid) %>%
#   mutate(c_div_maxc = catch / max(catch, na.rm = TRUE),
#          c_div_meanc = catch / mean(catch, na.rm = TRUE),
#          c_length = as.numeric(1:length(catch))) %>%
#   gather(metric, value, b_v_bmsy, u_v_umsy,exploitation_rate) %>%
#   select(stockid, year, contains('c_'), isscaap_group, metric, value) %>%
#   mutate(log_value = log(value + 1e-3)) %>%
#   unique() %>%
#   na.omit() %>%
#   ungroup() %>%
#   mutate(c_length = as.numeric(scale(c_length)))
#
#   # filter(!(metric == "mean_u" & value < 0.13 & enforcement > 0.6)) # per discussions with coauthors? don't remember why
# random_catch_tests <- catch_data %>%
#   nest(-metric) %>%
#   mutate(splits = map(data, ~rsample::vfold_cv(.x, v = 2, repeats = 1))) %>%
#   select(-data) %>%
#   unnest() %>%
#   mutate(sampid  = 1:nrow(.))
#
# model_structures <-
#   purrr::cross_df(list(
#     sampid = random_catch_tests$sampid,
#     model_structure = c(
#       "log_value ~ c_div_maxc + c_div_meanc + c_length")
#   ))
#
# random_catch_tests <- model_structures %>%
#   left_join(random_catch_tests, by = "sampid")
#
# random_catch_tests <- random_catch_tests %>%
#   mutate(
#     fit = map2(
#       splits,
#       model_structure,
#       fit_prior_regressions,
#       produce = "summary",
#       refresh = 100
#     )
#   )
#
# random_catch_tests <- random_catch_tests %>%
#   mutate(
#     training_performance = map(fit, "training_summary"),
#     testing_performance = map(fit, "testing_summary")
#   )
#
# random_catch_tests <- random_catch_tests %>%
#   mutate(
#     training_rmse = map_dbl(
#       training_performance,
#       ~ yardstick::rmse_vec(truth = .x$observed,
#                             estimate = .x$pp_pred)
#     ),
#     testing_rmse = map_dbl(
#       testing_performance,
#       ~ yardstick::rmse_vec(truth = .x$observed,
#                             estimate = .x$pp_pred)
#     )
#   )
#
# random_catch_tests %>%
#   ggplot(aes(model_structure, testing_rmse, color = metric)) +
#   geom_point() +
#   coord_flip()
#
# best_catch_tests <- random_catch_tests %>%
#   group_by(metric, model_structure) %>%
#   summarise(mean_rmse = mean(testing_rmse)) %>%
#   group_by(metric) %>%
#   filter(mean_rmse == min(mean_rmse))
#
# best_catch_tests <- best_catch_tests %>%
#   mutate(splits = map(metric, ~ filter(catch_data, metric == .x))) %>%
#   mutate(
#     best_fit = map2(
#       splits,
#       model_structure,
#       fit_prior_regressions,
#       produce = "results",
#       refresh = 100,
#       use_splits = FALSE
#     )
#   )
#
#
# best_catch_tests <- best_catch_tests %>%
#   mutate(fit = map(best_fit, "fit")) %>%
#   mutate(prior_plot = map2(fit, splits, plot_prior_fit))

