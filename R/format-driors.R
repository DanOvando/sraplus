#' Format data and priors for sraplus
#'
#' @param taxa the genus and species of the species (case insensitive)
#' @param initial_state reference point in the initial year, units of depletion or B/Bmsy (set ref_type accordingly). 
#' when \code{initial_state} is 1, implies either B/K = 1 or B/Bmsy = 1
#' @param initial_state_cv CV associated with initial state reference point
#' @param terminal_state reference point in the terminal year, units of depletion or B/Bmsy (set ref_type accordingly). 
#' when \code{initial_state} is 1, implies either B/K = 1 or B/Bmsy = 1
#' @param terminal_state_cv CV associated with terminal state reference point
#' @param u_v_umsy u/umsy data over time
#' @param u_years years in which u/umsy data are available
#' @param u_cv cv associated with u/umsy data
#' @param final_u vector of priors on u/umsy in the terminal years
#' @param final_u_cv vector of cvs on u/umsy in the terminal years
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
#' @param sigma_r_prior prior on process error
#' @param sigma_r_prior_cv CV of prior on process error
#' @param shape_prior prior on shape parameter of Pella-Tomlinson model
#' @param shape_prior_cv CV of prior on shape parameter of Pella-Tomlinson model
#' @param q_prior_cv CV of prior on q itself (prior on q set in \code{fit_sraplus})
#' @param sigma_obs_prior prior on observation error
#' @param sigma_obs_prior_cv cv of prior on observation error
#' @param b_ref_type 
#' @param f_ref_type 
#' @param use_b_reg 
#' @param q_slope_prior_cv 
#' @param isscaap_group 
#' @param prob 
#'
#' @return a list of data and priors
#' @export
#'
format_driors <-
  function(taxa = "lutjanus griseus",
           initial_state = 1,
           initial_state_cv = 0.1,
           terminal_state = NA,
           terminal_state_cv = 0.1,
           u_v_umsy = NA,
           u_years = NA,
           u_cv = 0.2,
           final_u = NA,
           final_u_cv = NA,
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
           growth_rate_prior_cv = 0.5,
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
           sigma_r_prior = 0.05,
           sigma_r_prior_cv = 0.5,
           shape_prior = 1.01,
           shape_prior_cv = 0.25,
           q_prior_cv = 1,
           sigma_obs_prior = 0.05,
           sigma_obs_prior_cv = .5,
           isscaap_group = NA,
           prob = 0.9) {
    
    if (use_heuristics == TRUE){
      
      warning("WARNING: You are using catch heursitics as your stock assessment")
      
    }
    
    genus <- tolower(strsplit(taxa, split = ' ')[[1]][1])
    
    if (is.na(isscaap_group)){
      
      
      if (any(grepl(tolower(taxa), tolower(fao_taxa$fao_species$scientific_name)))) {
        
        isscaap_group = fao_taxa$fao_species$isscaap_group[tolower(fao_taxa$fao_species$scientific_name) == tolower(taxa)][1]
        
      } else if (any(grepl(tolower(genus),tolower(fao_taxa$fao_genus$genus)))) {
        
        isscaap_group = fao_taxa$fao_genus$isscaap_group[tolower(fao_taxa$fao_genus$genus) == genus][1]
        
        
      } else{
        isscaap_group = "unknown"
      }
      
      
    }

    if (!is.na(sar)) {
      

      # isscaap_group = "Flounders, halibuts, soles"
      
      # tempmod <- best_sar_models$fit[best_sar_models$metric == "mean_uumsy"][[1]]
      
      if (f_ref_type == "fmsy"){
        
        tempmod <- best_sar_models$fit[best_sar_models$metric == "mean_uumsy"][[1]]
      } else if (f_ref_type == "f"){
        
        tempmod <- best_sar_models$fit[best_sar_models$metric == "mean_f"][[1]]
      }
      
      
      temp <- dplyr::tibble(sar = sar,
                            c_div_max_c = last(catch / max(catch)),
                            c_div_mean_c = last(catch / mean(catch)),
                            isscaap_group = isscaap_group,
                            sar_2 = sar^2)
      
      # factor_levels <- levels(tempmod$data$isscaap_group)
      
      # levels(temp$isscaap_group) <- factor_levels
      
      
      pp <-
        rstanarm::posterior_predict(tempmod, newdata = temp)
      
      pp <- pp[pp > quantile(pp, (1 - prob)/2) & pp < quantile(pp, 1 - (1 - prob)/2)]
      
      final_u <- c(final_u, exp(mean(pp)))
      
      usd <- ifelse(is.na(sar_cv), sd(pp), sar_cv)
      
      final_u_cv <- c(final_u_cv, usd)
      
      if (use_b_reg == TRUE){
        
        tempmod <- best_sar_models$fit[best_sar_models$metric == "mean_bbmsy"][[1]]
        
        pp <-
          rstanarm::posterior_predict(tempmod, newdata = temp)
        
        pp <- pp[pp > quantile(pp, (1 - prob)/2) & pp < quantile(pp, 1 - (1 - prob)/2)]
        
        terminal_state <- exp(mean(pp))
        
        terminal_state_cv <- ifelse(is.na(sar_cv), sd(pp), sar_cv)
      }
      
    }
    
    
    if (any(!is.na(fmi))) {
      
      
      temp <- purrr::map_df(fmi,  ~ . + 1e-6)
      
      temp$c_div_max_c = last(catch / max(catch))
      
      temp$c_div_mean_c = last(catch / mean(catch))
      
      temp$isscaap_group = isscaap_group
      if (f_ref_type == "fmsy"){
      
      tempmod <- best_fmi_models$fit[best_fmi_models$metric == "mean_uumsy"][[1]]
      } else if (f_ref_type == "f"){
        tempmod <- best_fmi_models$fit[best_fmi_models$metric == "mean_f"][[1]]
      }
    
    
      pp <-
        rstanarm::posterior_predict(tempmod, newdata = temp)
      
      pp <- pp[pp > quantile(pp, (1 - prob)/2) & pp < quantile(pp, 1 - (1 - prob)/2)]
      
      final_u <- c(final_u, exp(mean(pp)))
      
      final_u_cv <-
        c(final_u_cv, ifelse(is.na(fmi_cv), sd(pp), fmi_cv))
      
      if (use_b_reg == TRUE){
      
      tempmod <- best_fmi_models$fit[best_fmi_models$metric == "mean_bbmsy"][[1]]
      
      pp <-
        rstanarm::posterior_predict(tempmod, newdata = temp)
      
      pp <- pp[pp > quantile(pp, (1 - prob)/2) & pp < quantile(pp, 1 - (1 - prob)/2)]
      
      terminal_state <- exp(mean(pp))
      
      terminal_state_cv <- ifelse(is.na(fmi_cv), sd(pp), fmi_cv)
      }
      
    }
    
    
    genus_species <-
      taxa %>% stringr::str_split(" ", simplify = TRUE)
    
    shh <- purrr::safely(FishLife::Search_species)
    
    fish_search <-
      shh(Genus = genus_species[1], Species = genus_species[2])
    
    if (is.null(fish_search$error)) {
      taxon <- fish_search$result$match_taxonomy[1] %>%
        stringr::str_split("_") %>%
        unlist()
      
    } else{
      shh <- purrr::safely(FishLife::Search_species)
      
      taxon <-
        shh(Genus = genus_species[1])$result$match_taxonomy[1] %>%
        stringr::str_split("_") %>%
        unlist()
      
    }
    
    params_mvn <-
      c("r", "ln_var","M")
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
        ParentChild_gz = FishLifeData$ParentChild_gz
      )$result$match_taxonomy[1]
      
      taxa_location <-
        grep(sp, FishLifeData$ParentChild_gz[, "ChildName"])[1]
      
      mean_lh <- FishLifeData$beta_gv[taxa_location, ]
      
      cov_lh <- FishLifeData$Cov_gvv[taxa_location, , ]
      
      mean_lh <- mean_lh[which(names(mean_lh) %in% params_mvn)]
      
      cov_lh <-
        cov_lh[which(rownames(cov_lh) %in% params_mvn), which(colnames(cov_lh) %in% params_mvn)] %>% diag()
      
    } else {
      mean_lh <-
        c("r" = mean(FishLifeData$beta_gv[, "r"]),
          m = exp(mean(FishLifeData$beta_gv[, "M"])))
      
      cov_lh <-
        c("r" = sd(FishLifeData$beta_gv[, "r"]),
          m = sd(FishLifeData$beta_gv[, "M"]))
      
      
    }
    if (is.na(initial_state)) {
      initial_state <- ifelse(b_ref_type == "k", 1, 2.5)
    }
    
    if (is.na(initial_state)) {
      initial_state <- if (catch[1] / max(catch) < 0.2)
        c(0.7)
      else
        0.4
      
      initial_state_cv <- 0.1
    }
    
    
    if (any(!is.na(final_u))) {
      log_final_u <- log(final_u[!is.na(final_u)])
      
      log_final_u_cv <- final_u_cv[!is.na(final_u)]
      # 
      # if (is.na(terminal_state)){ # if no terminal b prior use approx from U/Umsy
      #   
      #   terminal_state <- pmax(.05,2.5 - mean(final_u[!is.na(final_u)]))
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
      log_final_u <- NA
      
      log_final_u_cv <- NA
    }
    
    if (use_heuristics == TRUE) {
      
      temp <-
        if (catch[1] / max(catch, na.rm = TRUE) < 0.2)
          0.7
      else
        0.4
      
      initial_state <-  dplyr::case_when(b_ref_type == "k" ~ temp,
                                     TRUE ~ temp * 2.5)
      
      initial_state_cv <- 0.2
      
      temp_terminal <-
        ifelse((dplyr::last(catch) / max(catch)) > 0.5, 0.6, 0.2)
      
      terminal_state <-
        dplyr::case_when(b_ref_type == "k" ~ temp_terminal,
                         TRUE ~ temp_terminal * 2.5)
      
      terminal_state_cv <- 0.2
      
    }
    
    if (!all(is.na(effort_years))){
      
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
    
    driors <-
      list(
        catch = catch,
        years = years,
        k_prior = ifelse(is.na(k_prior), 10 * max(catch), k_prior),
        k_prior_cv = k_prior_cv,
        terminal_state = terminal_state,
        terminal_state_cv = terminal_state_cv,
        initial_state = initial_state,
        initial_state_cv = initial_state_cv,
        index = index,
        effort = effort,
        u_v_umsy = u_v_umsy,
        u_years = u_years,
        u_cv = u_cv,
        index_years = index_years,
        effort_years = index_years,
        growth_rate_prior = ifelse(is.na(growth_rate_prior), mean_lh["r"],growth_rate_prior),
        growth_rate_prior_cv = ifelse(is.na(growth_rate_prior_cv),sqrt(cov_lh["r"]),growth_rate_prior_cv),
        sigma_r_prior = ifelse(is.na(sigma_r_prior),exp(mean_lh["ln_var"]) / 2, sigma_r_prior),
        sigma_r_prior_cv = ifelse(is.na(sigma_r_prior_cv), sqrt(cov_lh["ln_var"]), sigma_r_prior_cv),
        m =  ifelse(is.na(m), exp(mean_lh["M"]),m),
        log_final_u = log_final_u,
        log_final_u_cv = log_final_u_cv,
        q_slope_prior = q_slope_prior + 1e-6,
        q_slope_prior_cv = q_slope_prior_cv,
        shape_prior = shape_prior,
        shape_prior_cv = shape_prior_cv,
        q_prior_cv = q_prior_cv,
        sigma_obs_prior = sigma_obs_prior,
        sigma_obs_prior_cv = sigma_obs_prior_cv
      )
    
    driors$b_ref_type <- b_ref_type
    
    driors$f_ref_type <- f_ref_type
    
    return(driors)
    
  } # close function