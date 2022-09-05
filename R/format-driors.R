#' Format data and priors for sraplus
#'
#' @param taxa the genus and species of the species (case insensitive)
#' @param initial_state reference point in the initial year, units of depletion or B/Bmsy (set ref_type accordingly).
#' when \code{initial_state} is 1, implies either B/K = 1 or B/Bmsy = 1
#' @param initial_state_cv CV associated with initial state reference point
#' @param terminal_state reference point in the terminal year, units of depletion or B/Bmsy (set ref_type accordingly).
#' when \code{initial_state} is 1, implies either B/K = 1 or B/Bmsy = 1
#' @param terminal_state_cv CV associated with terminal state reference point
#' @param u u/umsy data over time
#' @param u_years years in which u/umsy data are available
#' @param u_cv cv associated with u/umsy data
#' @param terminal_u vector of priors on u/umsy in the terminal years
#' @param terminal_u_cv vector of cvs on u/umsy in the terminal years
#' @param catch vector of catches over lifetime of fishery
#' @param years vector of years that the catch data correspond to
#' @param index vector of an abundance index
#' @param effort vector of an effort series
#' @param index_years the years in which abundance index data are available
#' @param effort_years years in which effort data are available
#' @param use_heuristics logical,TRUE uses catch-msy hueristics for priors, FALSE requires user to pass them
#' @param fmi named vector of fisheries management index scrores
#' @param fmi_cv overwrite posterior predictive cv of FMI derived U/Umsy (NA to use empirical CV)
#' @param sar swept area ratio
#' @param sar_cv overwrite posterior predictive cv of SAR derived U/Umsy (NA to use empirical CV)
#' @param growth_rate_prior manually pass prior on the growth rate r in the Pella-Tomlinson model
#' @param growth_rate_prior_cv manually pass cv  for prior on the growth rate r in the Pella-Tomlinson model
#' @param q_slope_prior prior on the slope of catchability q
#' @param m natural mortality
#' @param k_prior prior on carrying capacity
#' @param k_prior_cv CV of prior on carrying capacity
#' @param shape_prior prior on shape parameter of Pella-Tomlinson model
#' @param shape_prior_cv CV of prior on shape parameter of Pella-Tomlinson model
#' @param q_prior_cv CV of prior on q itself (prior on q set in \code{fit_sraplus})
#' @param sigma_obs_prior prior on observation error
#' @param sigma_obs_prior_cv cv of prior on observation error
#' @param b_ref_type the units of the initial and final depletion. k means that units are relative to carrying capacity, b bmsy
#' @param f_ref_type the units of the fishing mortality rates. fmsy means F/Fmsy, f means fishing mortality rate f
#' @param use_b_reg logical marking whether to use the regression to put priors on biomass as well as fishing mortality
#' @param q_slope_prior_cv the cv on the prior on q_slope
#' @param isscaap_group the isscaap group of the stock. If NA will try and find
#' @param prob the probability range of the posterior predictive to use in construction of the priors
#' @param sigma_ratio_prior
#' @param sigma_ratio_prior_cv
#'
#' @return a list of data and priors
#' @export
#'
format_driors <-
  function(taxa = "lutjanus griseus",
           initial_state = 1,
           initial_state_cv = 0.1,
           terminal_state = NA,
           terminal_state_cv = 0.2,
           u = NA,
           u_years = NA,
           u_cv = 0.2,
           terminal_u = NA,
           terminal_u_cv = NA,
           catch = NA,
           years = NA,
           index = NA,
           effort = NA,
           b_ref_type = "k",
           f_ref_type = "fmsy",
           index_years = 1,
           effort_years = NA,
           use_heuristics = FALSE,
           use_b_reg = FALSE,
           growth_rate_prior = NA,
           growth_rate_prior_cv = 0.25,
           fmi = c(
             "research" = NA,
             "management" = NA,
             "enforcement" = NA,
             "socioeconomics" = NA
           ),
           fmi_cv = NA,
           sar = NA,
           sar_cv = NA,
           q_slope_prior = 0,
           q_slope_prior_cv = 0.1,
           m = NA,
           k_prior = NA,
           k_prior_cv = 2,
           sigma_ratio_prior = 1,
           sigma_ratio_prior_cv = 0.05,
           shape_prior = NA,
           shape_prior_cv = 0.25,
           q_prior = 1e-3,
           q_prior_cv = 1,
           sigma_obs_prior = 0.05,
           sigma_obs_prior_cv = 1,
           isscaap_group = NA,
           f_prior_form = 0,
           shape_prior_source = "thorson",
           prob = 0.9,
           use_fmsy_based_r = FALSE,
           use_catch_priors = FALSE) {

    # make sure vectors are in ascending year order
    
    year_order <- sort(years)
    
    year_index <- match(years, year_order)
    
    year <- sort(years)
    
    catch <- catch[year_index]
    
    if (use_b_reg == TRUE |
        (is.na(initial_state) &
         use_heuristics == FALSE) |
        (is.na(terminal_state) &
         use_heuristics == FALSE & use_catch_priors == TRUE)) {
      if (!is.na(initial_state) & b_ref_type == "k") {
        initial_state <-
          initial_state * (1 / shape_prior ^ (-1 / (shape_prior - 1)))
        
      }
      
      b_ref_type <-  "b"
      
    }
    
    if (use_heuristics == TRUE) {
      warning(
        "You are using catch heursitics as your stock assessment. Consider manually setting priors on terminal depletion based on expert opinion, or using a proxy such as swept area ratio"
      )
      
    }
    
    delta_years <- years - lag(years)
    
    if (any(delta_years[!is.na(delta_years)] > 1)) {
      stop(
        "You have missing years in your catch data. sraplus needs continuous catch data at this time: either interpolate missing years or locate missing data"
      )
    }
    
    genus <- tolower(strsplit(taxa, split = ' ')[[1]][1])
    
    if (is.na(isscaap_group)) {
      if (any(grepl(
        tolower(taxa),
        tolower(sraplus::fao_taxa$fao_species$scientific_name)
      ))) {
        isscaap_group = sraplus::fao_taxa$fao_species$isscaap_group[tolower(sraplus::fao_taxa$fao_species$scientific_name) == tolower(taxa)][1]
        
      } else if (any(grepl(
        tolower(genus),
        tolower(sraplus::fao_taxa$fao_genus$genus)
      ))) {
        isscaap_group = sraplus::fao_taxa$fao_genus$isscaap_group[tolower(sraplus::fao_taxa$fao_genus$genus) == genus][1]
        
        
      } else{
        isscaap_group = "unknown"
      }
      
      
    }
    
    
    # predict catch history cluster
    
    # create initial state prior
    
    if (is.na(initial_state) &
        use_heuristics == FALSE & length(catch) >= 25) {
      # catch <- ram_data$catch[ram_data$stockid == ram_data$stockid[1]]
      
      tmp <- data.frame(year = 1:length(catch), catch = catch) %>%
        dplyr::mutate(
          c_div_meanc = catch / mean(catch, na.rm = TRUE),
          log_fishery_length = log(length(catch))
        ) %>%
        dplyr::mutate(scaled_catch = as.numeric(scale(catch)))
      
      wide_tmp <- tmp %>%
        dplyr::select(year, scaled_catch) %>%
        tidyr::pivot_wider(names_from = year, values_from = scaled_catch) %>%
        janitor::clean_names()
      
      tmp$predicted_cluster <-
        parsnip::predict.model_fit(sraplus::cluster_fit, new_data = wide_tmp)[[1]]
      
      tmp <- tmp[tmp$year == 1, ]
      
      pred_init_state <-
        rstanarm::posterior_predict(sraplus::init_state_model,
                                    newdata = tmp,
                                    type = "response")
      
      initial_state <- exp(mean(pred_init_state))
      
      initial_state_cv <- sd(pred_init_state)
      
    }
    
    
    if (is.na(terminal_state) &
        use_heuristics == FALSE &
        length(catch) >= 25 & use_catch_priors == TRUE) {
      # catch <- ram_data$catch[ram_data$stockid == ram_data$stockid[1]]
      
      tmp <- data.frame(year = 1:length(catch), catch = catch) %>%
        dplyr::mutate(
          c_div_maxc = catch / max(catch, na.rm = TRUE),
          c_div_meanc = catch / mean(catch, na.rm = TRUE),
          c_length = log(length(catch)),
          fishery_year = 1:length(catch)
        ) %>%
        mutate(c_roll_meanc = RcppRoll::roll_meanr(c_div_meanc, 5))  %>%
        dplyr::mutate(scaled_catch = as.numeric(scale(catch)))
      
      wide_tmp <- tmp %>%
        dplyr::select(year, scaled_catch) %>%
        tidyr::pivot_wider(names_from = year, values_from = scaled_catch) %>%
        janitor::clean_names()
      
      tmp$predicted_cluster <-
        parsnip::predict.model_fit(sraplus::cluster_fit, new_data = wide_tmp)[[1]]
      
      tmp <- tmp[tmp$year == max(tmp$year), ]
      
      pred_init_state <-
        rstanarm::posterior_predict(sraplus::catch_b_model,
                                    newdata = tmp,
                                    type = "response")
      
      terminal_state <- exp(mean(pred_init_state))
      
      terminal_state_cv <- sd(pred_init_state)
      
    }
    
    
    if (!is.na(sar)) {
      if (f_ref_type == "fmsy") {
        tempmod <- sar_models$fit[sar_models$metric == "mean_uumsy"][[1]]
      } else if (f_ref_type == "f") {
        tempmod <- sar_models$fit[sar_models$metric == "mean_f"][[1]]
      }
      
      
      temp <- dplyr::tibble(
        sar = sar,
        c_div_max_c = dplyr::last(catch / max(catch)),
        c_div_mean_c = dplyr::last(catch / mean(catch)),
        isscaap_group = isscaap_group,
        sar_2 = sar ^ 2
      )
      
      pp <-
        rstanarm::posterior_predict(tempmod, newdata = temp)
      
      pp <-
        pp[pp > quantile(pp, (1 - prob) / 2) &
             pp < quantile(pp, 1 - (1 - prob) / 2)]
      
      terminal_u <- c(terminal_u, exp(mean(pp)))
      
      usd <- ifelse(is.na(sar_cv), sd(pp), sar_cv)
      
      terminal_u_cv <- c(terminal_u_cv, usd)
      
      if (use_b_reg == TRUE) {
        tempmod <- sar_models$fit[sar_models$metric == "mean_bbmsy"][[1]]
        
        pp <-
          rstanarm::posterior_predict(tempmod, newdata = temp)
        
        pp <-
          pp[pp > quantile(pp, (1 - prob) / 2) &
               pp < quantile(pp, 1 - (1 - prob) / 2)]
        
        terminal_state <- exp(mean(pp))
        
        terminal_state_cv <- ifelse(is.na(sar_cv), sd(pp), sar_cv)
      }
      
    }
    
    
    if (any(!is.na(fmi))) {
      temp <- purrr::map_df(fmi,  ~ . + 1e-6)
      
      temp$c_div_max_c = dplyr::last(catch / max(catch))
      
      temp$c_div_mean_c = dplyr::last(catch / mean(catch))
      
      temp$isscaap_group = isscaap_group
      if (f_ref_type == "fmsy") {
        tempmod <- fmi_models$fit[fmi_models$metric == "mean_uumsy"][[1]]
      } else if (f_ref_type == "f") {
        tempmod <- fmi_models$fit[fmi_models$metric == "mean_f"][[1]]
      }
      
      
      pp <-
        rstanarm::posterior_predict(tempmod, newdata = temp)
      
      pp <-
        pp[pp > quantile(pp, (1 - prob) / 2) &
             pp < quantile(pp, 1 - (1 - prob) / 2)]
      
      terminal_u <- c(terminal_u, exp(mean(pp)))
      
      terminal_u_cv <-
        c(terminal_u_cv, ifelse(is.na(fmi_cv), sd(pp), fmi_cv))
      
      if (use_b_reg == TRUE) {
        tempmod <- fmi_models$fit[fmi_models$metric == "mean_bbmsy"][[1]]
        
        pp <-
          rstanarm::posterior_predict(tempmod, newdata = temp)
        
        pp <-
          pp[pp > quantile(pp, (1 - prob) / 2) &
               pp < quantile(pp, 1 - (1 - prob) / 2)]
        
        terminal_state <- exp(mean(pp))
        
        terminal_state_cv <- ifelse(is.na(fmi_cv), sd(pp), fmi_cv)
      }
      
    }
    
    
    genus_species <-
      taxa %>% stringr::str_split(" ", simplify = TRUE)
    
    safe <- purrr::safely(sraplus::search_species)
    
    shh <- purrr::quietly(safe)
    
    fish_search <-
      shh(Genus = genus_species[1], Species = genus_species[2])$result
    
    if (is.null(fish_search$error)) {
      taxon <- fish_search$result$match_taxonomy[1] %>%
        stringr::str_split("_") %>%
        unlist()
      
    } else{
      taxon <-
        shh(Genus = genus_species[1])$result$result$match_taxonomy[1] %>%
        stringr::str_split("_") %>%
        unlist()
      
    }
    
    fishlife_taxa <-
      ifelse(is.null(taxon),
             "no_fishlife_match",
             paste(taxon, collapse = "_"))
    
    params_mvn <-
      c("r", "ln_var", "M", "ln_Fmsy")
    if (is.null(fish_search$error)) {
      taxon <-
        dplyr::tibble(
          Class = taxon[[1]],
          Order = taxon[[2]],
          Family = taxon[[3]],
          Genus = taxon[[4]],
          Species = taxon[[5]]
        )
      
      sp <- shh(
        Class = taxon["Class"],
        Order = taxon["Order"],
        Family = taxon["Family"],
        Genus = taxon["Genus"],
        Species = taxon["Species"],
      )$result$result$match_taxonomy[1]
      
      
      taxa_location <-
        grep(sp, sraplus::FishBase_and_RAM$ParentChild_gz[, "ChildName"])[1]
      
      # taxa_location <-
      #   grep(sp, FishLifeData$ParentChild_gz[, "ChildName"])[1]
      #
      
      mean_lh <- sraplus::FishBase_and_RAM$beta_gv[taxa_location,]
      
      # mean_lh <- FishLifeData$beta_gv[taxa_location,]
      
      cov_lh <-
        sraplus::FishBase_and_RAM$Cov_gvv[taxa_location, ,]
      
      # cov_lh <- FishLifeData$Cov_gvv[taxa_location, ,]
      
      mean_lh <- mean_lh[which(names(mean_lh) %in% params_mvn)]
      
      cov_lh <-
        cov_lh[which(rownames(cov_lh) %in% params_mvn), which(colnames(cov_lh) %in% params_mvn)] %>% diag()
      
      f_msy <- exp(mean_lh["ln_Fmsy"] + 0.5 * cov_lh["ln_Fmsy"])
      
      r_implied <-
        (f_msy / (1 - 1 / shape_prior)) * (shape_prior - 1)
      
      if (shape_prior_source == "thorson" & is.na(shape_prior)) {
        if (tolower(taxon$Order) %in% tolower(sraplus::ssb0msy_to_ssb0$tax_group)) {
          # pull MSY
          msy_k <-
            sraplus::ssb0msy_to_ssb0$ssbmsy_to_ssb0[tolower(sraplus::ssb0msy_to_ssb0$tax_group) == tolower(taxon$Order)]
          
          shape_prior_cv <-
            sraplus::ssb0msy_to_ssb0$ssbmsy_to_ssb0_sd[tolower(sraplus::ssb0msy_to_ssb0$tax_group) == tolower(taxon$Order)]
        } else {
          msy_k <-
            sraplus::ssb0msy_to_ssb0$ssbmsy_to_ssb0[tolower(sraplus::ssb0msy_to_ssb0$tax_group) == "other"]
          
          shape_prior_cv <-
            sraplus::ssb0msy_to_ssb0$ssbmsy_to_ssb0_sd[tolower(sraplus::ssb0msy_to_ssb0$tax_group) == "other"]
          
        }
        
        shapefoo <- function(log_m, msy_k) {
          m <- exp(log_m)
          
          msy_k_hat <- m ^ (-1 / (m - 1))
          
          obj <- (log(msy_k_hat) - log(msy_k)) ^ 2
          
        }
        
        implied_shape <-
          nlminb(log(1.001), shapefoo, msy_k = msy_k)
        
        shape_prior <- exp(implied_shape$par)
        
        
      } # close thorson shape prior
      
      if (shape_prior_source == "fishlife" & is.na(shape_prior)) {
        if (f_msy ==  mean_lh["r"]) {
          fudge <- 1e-3
        } else {
          fudge <- 0
        }
        
        shape_prior <-   as.numeric(mean_lh["r"]  / (f_msy + fudge))
        
      }
      
      
      if (use_fmsy_based_r & shape_prior_source != "fishlife") {
        mean_lh["r"] <- r_implied
        
      }
    } else {
      mean_lh <-
        c("r" = median(sraplus::FishBase_and_RAM$beta_gv[, "r"]),
          m = exp(median(
            sraplus::FishBase_and_RAM$beta_gv[, "M"]
          )))
      
      f_msy <-
        exp(median(sraplus::FishBase_and_RAM$beta_gv[, "ln_Fmsy"]))
      
      
      r_implied <-
        (f_msy / (1 - 1 / shape_prior)) * (shape_prior - 1)
      
      cov_lh <-
        c(
          "r" = sd(sraplus::FishBase_and_RAM$beta_gv[, "r"]),
          m = sd(sraplus::FishBase_and_RAM$beta_gv[, "M"])
        )
      
      if (shape_prior_source == "thorson" & is.na(shape_prior)) {
        msy_k <-
          sraplus::ssb0msy_to_ssb0$ssbmsy_to_ssb0[tolower(sraplus::ssb0msy_to_ssb0$tax_group) == "other"]
        
        shape_prior_cv <-
          sraplus::ssb0msy_to_ssb0$ssbmsy_to_ssb0_sd[tolower(sraplus::ssb0msy_to_ssb0$tax_group) == "other"]
        
        shapefoo <- function(log_m, msy_k) {
          m <- exp(log_m)
          
          msy_k_hat <- m ^ (-1 / (m - 1))
          
          obj <- (log(msy_k_hat) - log(msy_k)) ^ 2
          
        }
        
        implied_shape <-
          nlminb(log(1.001), shapefoo, msy_k = msy_k)
        
        shape_prior <- exp(implied_shape$par)
        
        
      } # close thorson shape prior
      
      if (shape_prior_source == "fishlife" & is.na(shape_prior)) {
        if (f_msy ==  mean_lh["r"]) {
          fudge <- 1e-3
        } else {
          fudge <- 0
        }
        
        shape_prior <-   mean_lh["r"]  / (f_msy + fudge)
        
      }
      
      if (use_fmsy_based_r & shape_prior_source != "fishlife") {
        mean_lh["r"] <- r_implied
        
      }
      
    } # close if/else match in fishlife
    
    
    
    
    if (is.na(initial_state)) {
      initial_state <-
        ifelse(b_ref_type == "k", 1, (1 / shape_prior ^ (-1 / (shape_prior - 1))))
    }
    
    if (is.na(initial_state)) {
      initial_state <- if (catch[1] / max(catch) < 0.2)
        c(0.7)
      else
        0.4
      
      initial_state_cv <- 0.1
    }
    
    
    if (any(!is.na(terminal_u))) {
      log_terminal_u <- log(terminal_u[!is.na(terminal_u)])
      
      log_terminal_u_cv <- terminal_u_cv[!is.na(terminal_u)]
      #
      # if (is.na(terminal_state)){ # if no terminal b prior use approx from U/Umsy
      #
      #   terminal_state <- pmax(.05,2.5 - mean(terminal_u[!is.na(terminal_u)]))
      #
      #   terminal_state_cv <- 0.2
      #
      #   if (ref_type == "k"){
      #
      #     ref_type = "b"
      #
      #     initial_state <- initial_state * 2.5
      #   }
      #
      # }
      
    } else {
      log_terminal_u <- NA
      
      log_terminal_u_cv <- NA
    }
    
    if (use_heuristics == TRUE) {
      temp <-
        if (catch[1] / max(catch, na.rm = TRUE) < 0.2)
          0.7
      else
        0.4
      
      initial_state <-  dplyr::case_when(b_ref_type == "k" ~ temp,
                                         TRUE ~ temp * (1 / shape_prior ^
                                                          (-1 / (shape_prior - 1))))
      
      initial_state_cv <-
        ifelse(is.na(initial_state_cv), 0.2, initial_state_cv)
      
      temp_terminal <-
        ifelse((dplyr::last(catch) / max(catch)) > 0.5, 0.6, 0.2)
      
      terminal_state <-
        dplyr::case_when(b_ref_type == "k" ~ temp_terminal,
                         TRUE ~ temp_terminal * (1 / shape_prior ^ (-1 / (shape_prior - 1))))
      
      terminal_state_cv <-
        ifelse(is.na(terminal_state_cv), 0.2, terminal_state_cv)
      
    }
    
    if (!all(is.na(effort_years))) {
      index_years = effort_years
      
    }
    
    # winkler life table translates one to the other
    # state space model one or the other variances goes to zero
    # can stabilize in some ways
    
    # see gausian prior on the sum of the variances and the raio of the variances
    # informative prior on the observation or process variance
    # observation error of the index
    # normally distributed on the log variance ratio, centered at zero
    # thorson munch ono 2014
    
    if (use_b_reg == TRUE) {
      log_terminal_u = NA
      
    }
    driors <-
      list(
        catch = catch,
        years = years,
        k_prior = ifelse(is.na(k_prior), 10 * max(catch), k_prior),
        k_prior_cv = sqrt(log(k_prior_cv ^ 2 + 1)),
        terminal_state = terminal_state,
        terminal_state_cv = sqrt(log(terminal_state_cv ^ 2 + 1)),
        initial_state = initial_state,
        initial_state_cv = sqrt(log(initial_state_cv ^ 2 + 1)),
        index = index,
        effort = effort,
        u = u,
        u_years = u_years,
        u_cv = u_cv,
        index_years = index_years,
        effort_years = index_years,
        growth_rate_prior = ifelse(is.na(growth_rate_prior), mean_lh["r"], growth_rate_prior),
        growth_rate_prior_cv = ifelse(
          is.na(growth_rate_prior_cv),
          sqrt(cov_lh["r"]),
          growth_rate_prior_cv
        ),
        sigma_ratio_prior = sigma_ratio_prior ,
        sigma_ratio_prior_cv = sqrt(log(sigma_ratio_prior_cv ^ 2 + 1)),
        # sigma_r_prior = ifelse(is.na(sigma_r_prior), exp(mean_lh["ln_var"]) / 2, sigma_r_prior),
        # sigma_r_prior_cv = ifelse(is.na(sigma_r_prior_cv), sqrt(cov_lh["ln_var"]), sigma_r_prior_cv),
        m =  ifelse(is.na(m), exp(mean_lh["M"]), m),
        log_terminal_u = log_terminal_u,
        log_terminal_u_cv =  sqrt(log(log_terminal_u_cv ^ 2 + 1)),
        terminal_u = exp(log_terminal_u),
        terminal_u_cv =  sqrt(log(log_terminal_u_cv ^ 2 + 1)),
        q_slope_prior = q_slope_prior,
        q_slope_prior_cv = sqrt(log(q_slope_prior_cv ^ 2 + 1)),
        shape_prior = shape_prior,
        shape_prior_cv =  sqrt(log(shape_prior_cv ^ 2 + 1)),
        q_prior = q_prior,
        q_prior_cv = sqrt(log(q_prior_cv ^ 2 + 1)),
        sigma_obs_prior = sigma_obs_prior,
        sigma_obs_prior_cv = sqrt(log(sigma_obs_prior_cv ^ 2 + 1)),
        fishlife_taxa = fishlife_taxa,
        input_taxa = taxa,
        f_prior_form = f_prior_form
      )
    
    
    driors$b_ref_type <- b_ref_type
    
    driors$f_ref_type <- f_ref_type
    
    return(driors)
    
  } # close function