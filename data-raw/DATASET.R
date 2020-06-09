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
library(kernlab)
library(tidymodels)

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


# load(here("data-raw","Return.Rdata"))

# readr::write_rds(Return,here::here("data-raw","Return.rds"))
# 
# rm(Return)

# Return <- readr::read_rds(here::here("data-raw","Return.rds"))


# FishLifeData<- Return[c("ParentChild_gz","beta_gv","Cov_gvv")]
# 
# FishLifeData$metadata <- "emailed from Thorson, beta version with newparameters"


if (!file.exists(here("data-raw","ram.zip"))) {
  
  download.file("https://www.dropbox.com/s/jpgz0a5s5of3qev/RAM%20v4.491%20Files%20(1-14-20).zip?dl=1", destfile = here::here("data-raw","ram.zip"))
  
  unzip(here::here("data-raw","ram.zip"), exdir = "data-raw")
  
}



# process RAM data --------------------------------------------------------

ram_dirs <- list.files("data-raw")

ram_dirs <- ram_dirs[str_detect(ram_dirs,"RAM v\\d")]

ram_files <- list.files(ram_dirs, recursive = TRUE)

ram_files <- ram_files[str_detect(ram_files,".RData")]

ram_files <- ram_files[str_detect(ram_files,"Model Fit")]

load(here(ram_dirs,ram_files[1]))

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

# filter data 

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
    b_rel = dplyr::case_when(
      has_tb0 ~ total_biomass / max(TB0),
      has_tb ~ total_biomass / max(total_biomass),
      TRUE ~ b_v_bmsy / 2
    )
  ) %>%
  ungroup()




# process data ------------------------------------------------------------


cod <- ram_data %>% 
  filter(stocklong.x == "Atlantic cod ICES 3a(west)-4-7d")
usethis::use_data(cod, overwrite = TRUE, internal = FALSE)


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

fao_taxa <- list(fao_species = fao_species,
                 fao_genus = fao_genus)

# usethis::use_data(fao_taxa,FishLifeData, overwrite = TRUE, internal = TRUE)

# usethis::use_data(fao_taxa,FishLifeData, overwrite = TRUE, internal = TRUE)

usethis::use_data(fao_taxa, overwrite = TRUE)

# usethis::use_data(FishLifeData, overwrite = TRUE, internal = TRUE)



fao$fao_country_name <-
  countrycode::countrycode(fao$country, "country.name", "un.name.en")

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


# load fmi data ----------------------------------------------------------------


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


ram_fmi_linkages <-
  readxl::read_xlsx(
    here::here("data-raw", "RAM-FMI linkages for DO 2019-09-16.xlsx"),
    sheet = "RAM-FMI linkages",
    skip = 1
  ) %>%
  janitor::clean_names() %>%
  select(-contains("primary")) %>%
  select(stockid,
         areaid,
         iso3_pref,
         scientificname,
         fao_scientific_name) %>%
  mutate(iso3_pref = stringr::str_trim(stringr::str_replace_all(iso3_pref, ' \\| ', "\\|")))


ram_fmi_fao_regions_linkages <-
  readxl::read_xlsx(
    here::here("data-raw", "RAM-FMI linkages for DO 2019-09-16.xlsx"),
    sheet = "csv",
    skip = 0
  ) %>%
  janitor::clean_names() %>%
  mutate(lookup_code = stringr::str_trim(stringr::str_replace_all(cspref, ' \\| ', "\\|"))) %>%
  select(-contains("_recent"), stockid, region, primary_country, faor)


fmi <-
  readxl::read_xlsx(here::here("data-raw", "FMI data extract by stock for Ray 2018-11-05.xlsx"),
                    sheet = "summary"
  ) %>%
  janitor::clean_names()


fmi$fao_country_name <-
  countrycode::countrycode(fmi$country_rfmo, "country.name", "un.name.en")

fmi$region <-   countrycode::countrycode(fmi$country_rfmo, "country.name", "region")

# fmi <- fmi %>%
#   mutate(country_rfmo = case_when(is.na(fao_country_name) ~ country_rfmo, TRUE ~ fao_country_name))
# 
# fmi$country <- fmi$country_rfmo
# 
# fmi$country_rfmo <-
#   case_when(
#     (fmi$region == "Southern Europe"  &
#        (fmi$basin == "Med" |
#           is.na(fmi$basin))) |
#       (fmi$basin == "Med" | is.na(fmi$basin))  ~ "GFCM",
#     fmi$region %in%  c("Northern Europe", "Western Europe") |
#       (fmi$region == "Southern Europe" &
#          !(fmi$basin == "Med" | is.na(fmi$basin))) ~ "ICES",
#     TRUE ~ fmi$country_rfmo
#   )

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
  left_join(fmi_scinames %>% rename(scientificname2 = scientificname), by = "species")

fmi$lookup_code <- stringr::str_trim(stringr::str_replace_all(fmi$lookup_code,' \\| ', "\\|"))

fmi <- fmi %>% 
  left_join(ram_fmi_linkages, by = c("lookup_code" = "iso3_pref")) %>% 
  left_join(ram_fmi_fao_regions_linkages, by = c("stockid", "lookup_code")) %>% 
  mutate(scientificname = ifelse(is.na(scientificname), scientificname2, scientificname)) %>% 
  select(-scientificname2)


fmi <-
  fmi %>%
  mutate(genus = map_chr(scientificname,~str_split(.x," ", simplify = TRUE)[,1])) %>%
  left_join(fao_species, by = c("scientificname" = "scientific_name")) %>%
  left_join(fao_genus, by = "genus") %>%
  mutate(isscaap_group = ifelse(is.na(isscaap_group.x), isscaap_group.y, isscaap_group.x)) %>%
  select(-isscaap_group.x,-isscaap_group.y)

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
            c_div_max_c = mean(c_maxc)) %>%
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
    c_div_max_c = mean(c_maxc)) %>%
  gather(metric, value,-stockid, -c_div_max_c) %>%
  ungroup() %>%
  left_join(fmi, by = "stockid") %>%
  filter(!is.na(lookup_code)) %>%
  select(-basin, -stock, -scientificname,-species,-contains("fao"),-contains(".x")) %>%
  mutate(log_value = log(value)) %>%
  unique() %>% 
  na.omit() %>% 
  mutate_at(c("research", "management", "enforcement", "socioeconomics"), ~ .x + 1e-6) %>% 
  filter(isscaap_group != "Tunas, bonitos, billfishes")



ram_v_fmi %>% 
  filter(metric == "mean_uumsy") %>% 
  ggplot(aes(management, socioeconomics, color = value)) + 
  geom_point(size = 4) + 
  theme_classic() + 
  scale_color_viridis_c()

ram_v_fmi %>% 
  filter(metric == "mean_bbmsy") %>% 
  ggplot(aes(management, enforcement, color = value)) + 
  geom_point(size = 4) + 
  theme_classic() + 
  scale_color_viridis_c()

ram_v_fmi %>% 
  filter(metric == "mean_uumsy") %>% 
  ggplot(aes(socioeconomics, value, color = value)) + 
  geom_point(size = 4) + 
  geom_smooth(method = "lm") +
  theme_classic() + 
  scale_color_viridis_c()


ram_v_fmi %>%
  filter(metric == "mean_uumsy") %>%
  pivot_longer(
    c(mean_fmi, research, socioeconomics, management, enforcement),
    names_to = "fmi_metric",
    values_to = "fmi_score"
  ) %>%
  ggplot(aes(fmi_score, value, color = value)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm") +
  theme_ipsum() +
  scale_color_viridis_c() +
  facet_wrap( ~ fmi_metric)

ram_v_fmi %>%
  filter(metric == "mean_uumsy") %>%
  pivot_longer(
    c(mean_fmi, research, socioeconomics, management, enforcement),
    names_to = "fmi_metric",
    values_to = "fmi_score"
  ) %>%
  ggplot(aes(fmi_score, value, color = value)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm") +
  theme_ipsum() +
  scale_color_viridis_c() +
  facet_wrap( ~ fmi_metric)




ram_fmi <- ram_v_fmi

usethis::use_data(ram_fmi, overwrite = TRUE)

a = fmi %>% 
  filter(!is.na(lookup_code))

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

ram_v_fmi %>%
  ggplot(aes(mean_fmi, log_value)) +
  geom_point() +
  facet_wrap(~metric)  +
  geom_smooth()

random_fmi_tests <- ram_v_fmi %>%
  nest(-metric) %>%
  mutate(splits = map(data, ~ rsample:: vfold_cv(.x, v = 3, repeats = 5))) %>%
  select(-data) %>%
  unnest() %>%
  mutate(sampid  = 1:nrow(.))


# fmi_splits <- initial_split(ram_v_fmi %>% filter(metric == "mean_bbmsy") %>% 
#                        select(value,research, management,enforcement, socioeconomics, c_div_max_c,country_rfmo), prop = 0.75)
# 
# training_data <- training(fmi_splits)
# 
# testing_data <- testing(fmi_splits)
# 
# aa_splits <- rsample::group_vfold_cv(training_data, group = country_rfmo)
# 
# tune_grid <-
#   parameters(
#     min_n(range(2, 10)),
#     tree_depth(range(4, 15)),
#     learn_rate(range = c(-2, 0)),
#     mtry(),
#     loss_reduction(),
#     sample_size(range = c(1, 1)),
#     trees(range = c(500, 2000))
#   ) %>%
#   dials::finalize(mtry(), x = training_data %>% select(-(1:2)))
# 
# xgboost_grid <- grid_max_entropy(tune_grid, size = 30) 
# 
# xgboost_model <-
#   parsnip::boost_tree(
#     mode = "regression",
#     mtry = tune(),
#     min_n = tune(),
#     loss_reduction = tune(),
#     sample_size =tune(), 
#     learn_rate = tune(),
#     tree_depth = tune(),
#     trees =tune()
#   ) %>%
#   parsnip::set_engine("xgboost") 
# 
# 
# xgboost_workflow <- workflows::workflow() %>% 
#   add_formula(value ~ research + management + enforcement + socioeconomics + c_div_max_c) %>% 
#   add_model(xgboost_model)
# 
# 
# set.seed(234)
# doParallel::registerDoParallel(cores = parallel::detectCores() - 2)
# 
# xgboost_tuning <- tune_grid(
#   xgboost_workflow,
#   resamples = aa_splits,
#   grid = xgboost_grid,
#   control = control_grid(save_pred = TRUE)
# )
# 
# # xgboost_tuning
# 
# collect_metrics(xgboost_tuning) %>%
#   select(mean, mtry:sample_size, .metric) %>%
#   pivot_longer(mtry:sample_size, names_to = "dial", values_to = "level") %>%
#   ggplot(aes(level, mean)) +
#   geom_point() +
#   facet_grid(.metric ~ dial, scales = "free")
# 
# # show_best(xgboost_tuning, "rmse")
# # 
# 
# best_rmse <- tune::select_best(xgboost_tuning, metric = "rmse")
# 
# final_workflow <- finalize_workflow(
#   xgboost_workflow,
#   best_rmse
# )
# 
# fmi_fit <- last_fit(
#   final_workflow,
#   fmi_splits
# )
# 
# fmi_fit$.predictions[[1]] %>% 
#   ggplot(aes(value, .pred)) + 
#   geom_point()
# 
# fmi_fit$.metrics
# 
# 
# huh <- stan_glm(value ~ research + management + enforcement + socioeconomics + c_div_max_c, 
#                 data = training_data, family = Gamma)
# 
# pp_huh <- posterior_predict(huh, newdata = testing_data, type = "response") %>% 
#   colMeans()
# 
# fmi_fit$.predictions[[1]]$stan_fit <- pp_huh
# 
# fmi_fit$.predictions[[1]] %>% 
#   pivot_longer(c(.pred,stan_fit),names_to = "model", values_to = "prediction") %>% 
#   ggplot(aes(value, prediction, color = model)) + 
#   geom_point()
# 
# 
# fmi_fit$.predictions[[1]] %>% 
#   pivot_longer(c(.pred,stan_fit),names_to = "model", values_to = "prediction") %>% 
#   group_by(model) %>% 
#   yardstick::rsq(truth = value, estimate = prediction)

model_structures <-
  purrr::cross_df(list(
    sampid = random_fmi_tests$sampid,
    model_structure = c(
      "log_value ~ research + management + enforcement + socioeconomics + c_div_max_c" ,
      # "log_value ~ (research + management + enforcement + socioeconomics + c_div_max_c - 1|isscaap_group)",
      # "log_value ~ c_div_max_c  + (research + management + enforcement + socioeconomics - 1|isscaap_group)",
      "log_value ~ research + management + enforcement + socioeconomics",
      # "log_value ~ (research + management + enforcement + socioeconomics - 1|isscaap_group)",
      "log_value ~ + management + enforcement + socioeconomics + c_div_max_c"
      # "log_value ~ (log(research) + log(management) + log(enforcement) + log(socioeconomics) - 1|isscaap_group)"
      
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
      refresh = 500,
      iter = 2000
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
      use_splits = FALSE,
      iter = 5000
    )
  )

# random_fmi_tests %>% 
#   ggplot(aes(model_structure, testing_rmse)) + 
#   geom_violin() + 
#   coord_flip() + 
#   facet_wrap(~metric) + 
#   theme_minimal()




fmi_models <- best_fmi_models %>%
  mutate(fit = map(best_fmi_fit, "fit")) %>% 
  # mutate(prior_plot = map2(fit, splits, plot_prior_fit)) %>% 
  select(-best_fmi_fit)

usethis::use_data(fmi_models,overwrite = TRUE)

# sraplus::fmi_models$fit[[1]] %>% rstanarm::launch_shinystan()
# fit sar models --------------------------------------------------------------

ram_v_sar <- sar_to_ram %>%
  gather(metric, value, contains("mean_"),-mean_stock_in_tbp,-mean_tbp_in_stock) %>%
  mutate(log_value = log(value + 1e-3)) %>%
  mutate(sar_2 = sar ^ 2) %>%
  select(stockid, sar, sar_2, isscaap_group, metric, value, log_value, c_div_max_c,mean_stock_in_tbp) %>% 
  filter(mean_stock_in_tbp > 25 | is.na(mean_stock_in_tbp)) %>% 
  select(-mean_stock_in_tbp) %>% 
  na.omit()

ram_v_sar %>% 
  ggplot(aes(log(sar), log_value)) + 
  geom_point() + 
  facet_wrap(~metric)

ram_v_sar <- recipe(log_value ~ ., data = ram_v_sar) %>%
  step_other(isscaap_group) %>%
  prep(data = sar_data, retain = TRUE) %>%
  juice()


huh <- stan_gamm4(log_value ~ s(sar) + s(c_div_max_c), data = ram_v_sar)


random_sar_tests <- ram_v_sar %>%
  nest(-metric) %>%
  mutate(splits = map(data, ~ rsample::vfold_cv(.x, v = 3, repeats = 5))) %>%
  select(-data) %>%
  unnest() %>%
  mutate(sampid  = 1:nrow(.))

model_structures <-
  purrr::cross_df(list(
    sampid = random_sar_tests$sampid,
    model_structure = c(
      # "log_value ~  (log(sar) + log(sar_2) - 1|isscaap_group)",
      "log_value ~ poly(sar,2) + c_div_max_c",
      "log_value ~ log(sar) + c_div_max_c",
      "log_value ~ sar + c_div_max_c"
      # "log_value ~ c_div_max_c + (sar - 1|isscaap_group)",
      # "log_value ~ c_div_max_c + (sar + sar_2 - 1|isscaap_group)",
      # "log_value ~  (sar + sar_2 - 1|isscaap_group)",
      # "log_value ~ c_div_max_c + (log(sar) - 1|isscaap_group)"
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
    refresh = 500,
    iter = 2000
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
      use_splits = FALSE,
      iter = 5000
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

sar_models <- best_sar_models %>% 
  select(-best_fit, -prior_plot)

usethis::use_data(sar_models,overwrite = TRUE)



# catch based initial state prior ----------------------------------------------

# random idea: can you classify catch series into groups and use that as a predictor?

ram_catches <- ram_data %>% 
  ungroup() %>% 
  select(stockid, year, catch) %>% 
  group_by(stockid) %>% 
  mutate(stock_year = 1:length(catch),
         n_years = length(catch)) %>% 
  mutate(scaled_catch = scale(catch)) %>% 
  ungroup() %>% 
  filter(n_years > 25,
         stock_year <= 25)

ram_catches %>% 
  ggplot(aes(stock_year, scaled_catch, color = stockid)) + 
  geom_line(show.legend = FALSE)

ram_catches <- ram_catches %>% 
  select(stockid, stock_year, scaled_catch) %>% 
  pivot_wider(names_from = stock_year, values_from = scaled_catch) %>% 
  ungroup() 

nstocks <- nrow(ram_catches)

map_dbl(ram_catches, ~sum(is.na(.x)))

a = ram_catches %>% select(-stockid) %>% as.matrix()

catch_pca <- specc(a,centers = 8)

centers(catch_pca)
size(catch_pca)
withinss(catch_pca)

cluster <- as.vector(catch_pca)

ram_catches$cluster <- cluster


ram_catches <- ram_catches  %>%
  pivot_longer(
    c(-stockid, -cluster),
    names_to = "stock_year",
    values_to = "catch",
  ) %>% 
  mutate(stock_year = as.integer(stock_year))

ram_catches %>% 
  ggplot(aes(stock_year, catch, group = stockid)) + 
  geom_line() + 
  facet_wrap(~cluster)



class_data <- ram_catches %>%
  pivot_wider(names_from = stock_year, values_from = catch) %>%
  ungroup() %>%
  mutate(cluster = as.factor(cluster))

class_splits <- rsample::initial_split(class_data,strata = cluster)

class_model <- parsnip::nearest_neighbor(neighbors = 15,
                                         mode = "classification") %>% 
  set_engine(engine = 'kknn')

class_fit <- class_model %>% 
  parsnip::fit(cluster ~ . , data = training(class_splits) %>% select(-stockid))

class_fit$fit %>% summary()

training_data <- training(class_splits) %>% 
  bind_cols(predict(class_fit, new_data = .)) %>% 
  mutate(split = "training")


testing_data <- testing(class_splits) %>% 
  bind_cols(predict(class_fit, new_data = .)) %>% 
  mutate(split = "testing")

cluster_predictions <- training_data %>% 
  bind_rows(testing_data) %>% 
  rename(predicted_cluster = .pred_class)

cluster_model_performance <- cluster_predictions %>%
  group_by(split, cluster) %>% 
  summarise(accuracy = mean(cluster == predicted_cluster))

cluster_model_performance %>% 
  ggplot(aes(cluster, accuracy, fill = split)) + 
  geom_col(position = "dodge")
 
cluster_predictions %>% 
  group_by(split) %>% 
  summarise(accuracy = mean(cluster == predicted_cluster)) %>% 
  pivot_wider(names_from = "split", values_from = "accuracy") %>% 
  mutate(testing_loss = testing / training - 1)



# fit a model now

status_model_data <- ram_data %>% 
  filter(stockid %in% unique(cluster_predictions$stockid)) %>% 
  left_join(cluster_predictions %>% select(stockid, predicted_cluster), by = "stockid") %>% 
  left_join(sraplus::fao_taxa$fao_species %>% select(scientific_name, isscaap_group), by = c("scientificname" = "scientific_name")) %>% 
  group_by(stockid) %>%
  mutate(c_div_maxc = catch / max(catch, na.rm = TRUE),
         c_div_meanc = catch / mean(catch, na.rm = TRUE),
         c_length = log(length(catch)),
         fishery_year = 1:length(catch)) %>%
  mutate(c_roll_meanc = RcppRoll::roll_meanr(c_div_meanc,5)) %>% 
  gather(metric, value, b_v_bmsy, u_v_umsy,exploitation_rate) %>%
  select(stockid, year, contains('c_'), metric, value, predicted_cluster,fishery_year) %>%
  mutate(log_value = log(value + 1e-3)) %>%
  unique() %>%
  na.omit() %>%
  ungroup() 

# try and predict initial status


initial_state_splits <-
  initial_split(
    status_model_data %>%  group_by(stockid) %>% filter(metric == "u_v_umsy") %>% ungroup()
  )

a <- status_model_data %>%  group_by(stockid) %>% filter(metric == "u_v_umsy") %>% ungroup()

training_stocks <- sample(unique(a$stockid), round(n_distinct(a$stockid) * .75), replace = FALSE)

a$training <- !a$stockid %in% training_stocks

a <- a %>% 
  arrange(training, stockid, year)

training <- a %>% 
  filter(stockid %in% training_stocks)

testing <- a %>% 
  filter(!stockid %in% training_stocks)


catch_split <- initial_time_split(a, prop = last(which(a$training == FALSE)) / nrow(a))

training_data <- training(catch_split)

testing_data <- testing(catch_split)

aa_splits <- rsample::group_vfold_cv(training_data, group = stockid)

tune_grid <-
  parameters(
    min_n(range(2, 10)),
    tree_depth(range(4, 15)),
    learn_rate(range = c(-2, 0)),
    mtry(),
    loss_reduction(),
    sample_size(range = c(1, 1)),
    trees(range = c(500, 2000))
  ) %>%
  dials::finalize(mtry(), x = training_data %>% select(-(1:2)))

xgboost_grid <- grid_max_entropy(tune_grid, size = 10)

xgboost_model <-
  parsnip::boost_tree(
    mode = "regression",
    mtry = tune(),
    min_n = tune(),
    loss_reduction = tune(),
    sample_size =tune(),
    learn_rate = tune(),
    tree_depth = tune(),
    trees =tune()
  ) %>%
  parsnip::set_engine("xgboost")


xgboost_workflow <- workflows::workflow() %>%
  add_formula(value ~ c_div_maxc + fishery_year + predicted_cluster + c_div_meanc + c_roll_meanc) %>%
  add_model(xgboost_model)


set.seed(234)
doParallel::registerDoParallel(cores = parallel::detectCores() - 2)

xgboost_tuning <- tune_grid(
  xgboost_workflow,
  resamples = aa_splits,
  grid = xgboost_grid,
  control = control_grid(save_pred = TRUE)
)

# xgboost_tuning

collect_metrics(xgboost_tuning) %>%
  select(mean, mtry:sample_size, .metric) %>%
  pivot_longer(mtry:sample_size, names_to = "dial", values_to = "level") %>%
  ggplot(aes(level, mean)) +
  geom_point() +
  facet_grid(.metric ~ dial, scales = "free")

# show_best(xgboost_tuning, "rmse")
#

best_rmse <- tune::select_best(xgboost_tuning, metric = "rmse")

final_workflow <- finalize_workflow(
  xgboost_workflow,
  best_rmse
)



catch_model <-
  parsnip::boost_tree(
    mode = "regression",
    mtry = best_rmse$mtry,
    min_n = best_rmse$min_n,
    loss_reduction = best_rmse$loss_reduction,
    sample_size = best_rmse$sample_size, 
    learn_rate = best_rmse$learn_rate,
    tree_depth = best_rmse$tree_depth,
    trees = best_rmse$trees
  ) %>%
  parsnip::set_engine("xgboost") %>%
  parsnip::fit(value ~ c_div_maxc + fishery_year + predicted_cluster + c_div_meanc + c_roll_meanc, data = training_data)


catch_model %>%
  vip::vi() %>% 
  vip::vip(geom = "point")


training_data$boost_fit <- predict(catch_model, new_data = training_data)$.pred

testing_data$boost_fit <- predict(catch_model, new_data = testing_data)$.pred

catch_model_predictions <- training_data %>% 
  mutate(set = "training") %>% 
  bind_rows(testing_data %>% mutate(set = "testing")) 

catch_model_predictions %>% 
  ggplot(aes(value, boost_fit)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~set)


huh <- stan_glm(value ~ research + management + enforcement + socioeconomics + c_div_max_c,
                data = training_data, family = Gamma)

pp_huh <- posterior_predict(huh, newdata = testing_data, type = "response") %>%
  colMeans()

fmi_fit$.predictions[[1]]$stan_fit <- pp_huh

fmi_fit$.predictions[[1]] %>%
  pivot_longer(c(.pred,stan_fit),names_to = "model", values_to = "prediction") %>%
  ggplot(aes(value, prediction, color = model)) +
  geom_point()


fmi_fit$.predictions[[1]] %>%
  pivot_longer(c(.pred,stan_fit),names_to = "model", values_to = "prediction") %>%
  group_by(model) %>%
  yardstick::rsq(truth = value, estimate = prediction)




status_model_data %>%  
  group_by(stockid) %>% 
  filter(fishery_year == max(fishery_year), metric == "b_v_bmsy") %>% 
  ungroup() %>% 
  ggplot(aes(log(c_div_meanc), log(value))) + 
  geom_point() + 
  geom_smooth(method = "lm")

status_model_data %>%  
  group_by(stockid) %>% 
  filter(fishery_year == max(fishery_year), metric == "b_v_bmsy") %>% 
  ungroup() %>% 
  ggplot(aes(c_div_meanc, c_div_maxc, color = value)) + 
  geom_point() + 
  scale_color_viridis_c()


with_cluster <-
  stan_glmer(
    value ~  c_div_meanc + c_div_maxc + (1 |predicted_cluster),
    data =  training,
    family = gaussian(link = "log")
  )




witho_cluster <-
  stan_glm(
    value ~ predicted_cluster + c_div_meanc + c_div_maxc + c_length + log(fishery_year),
    data =  training(initial_state_splits),
    family = gaussian(link = "log")
  )

with_cluster$loo <- loo(with_cluster)

witho_cluster$loo <- loo(witho_cluster)

model_list <- stanreg_list(with_cluster, witho_cluster, model_names = c("With Clusters", "Without Clusters"))
loo_compare(model_list)


# 
# 
with_cluster <-
  rand_forest(trees = 1000) %>%
  set_engine("ranger") %>%
  set_mode("regression")

test_with_cluster <- with_cluster %>%
  fit(
    value ~ .,
    data = training(initial_state_splits) %>% select(value, c_div_maxc, c_div_meanc, c_length, predicted_cluster, fishery_year, isscaap_group)
  )

# test_rf$fit
plot(with_cluster)

training_with_cluster <- training(initial_state_splits) %>% 
  mutate(predicted_value = predict(with_cluster, type = "response")) %>% 
  mutate(split = "training")

pp <- posterior_predict(with_cluster, newdata =  testing(initial_state_splits), type = "response")

testing_with_cluster <- testing(initial_state_splits) %>% 
  mutate(predicted_value = colMeans(pp)) %>% 
  mutate(split = "testing")



with_cluster_performance <- training_with_cluster %>% 
  bind_rows(testing_with_cluster)

with_cluster_performance %>% 
  ggplot(aes(value, predicted_value)) + 
  geom_point() +
  geom_abline(aes(slope = 1, intercept = 0)) + 
  geom_smooth(method = "lm") +
  facet_wrap(~split)

with_cluster_performance %>% 
  group_by(split) %>% 
  yardstick::rmse(truth = value, estimate = predicted_value)

without_cluster <-
  stan_glm(
    value ~ c_div_maxc + c_div_meanc +  c_length ,
    data =  training(initial_state_splits),
    family = Gamma
  )

plot(without_cluster)

training_without_cluster <- training(initial_state_splits) %>% 
  mutate(predicted_value = predict(without_cluster, type = "response")) %>% 
  mutate(split = "training")

pp <- posterior_predict(without_cluster, newdata =  testing(initial_state_splits), type = "response")

testing_without_cluster <- testing(initial_state_splits) %>% 
  mutate(predicted_value = colMeans(pp)) %>% 
  mutate(split = "testing")

without_cluster_performance <- training_without_cluster %>% 
  bind_rows(testing_without_cluster)

without_cluster_performance %>% 
  ggplot(aes(value, predicted_value)) + 
  geom_point() +
  geom_abline(aes(slope = 1, intercept = 0)) + 
  geom_smooth(method = "lm") +
  facet_wrap(~split)

without_cluster_performance %>% 
  group_by(split) %>% 
  yardstick::rmse(truth = value, estimate = predicted_value)




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

