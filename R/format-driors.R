#' Format data and priors for sraplus
#'
#' @param taxa the genus and species of the species (case insensitive)
#' @param initial_b b reference point in the initial year, units of depletion of B/Bmsy (set ref_type accordingly)
#' @param initial_b_sd sigma associated with initial b reference point
#' @param terminal_b b reference point in the terminal year, units of depletion of B/Bmsy (set ref_type accordingly)
#' @param terminal_b_sd sigma associated with terminal b reference point
#' @param carry prior on carrying capacity, deprecated
#' @param carry_sd cv associated with prior on carrying capacity, deprecated
#' @param u_v_umsy u/umsy data over time
#' @param u_years years in which u/umsy data are available
#' @param u_sd cv associated with u/umsy data
#' @param final_u vector of priors on u/umsy in the terminal years
#' @param final_u_sd vector of cvs on u/umsy in the terminal years
#' @param catch vector of catches over lifetime of fishery
#' @param years vector of years that the catch data correspond to
#' @param index vector of an abundance index
#' @param effort vector of an effort series
#' @param ref_type k if initial and final depletions are in depletion units, b if in b/bmsy units
#' @param index_years the years in which abundance index data are available
#' @param effort_years years in which effort data are available
#' @param use_heuristics logical,TRUE uses catch-msy hueristics for priors, FALSE requires user to pass them
#' @param fmi named vector of fisheries management index scrores
#' @param fmi_sd overwrite fmi prediction sd
#' @param sar swept area ratio
#' @param sar_sd overwrite sar prediction sd
#' @param f_sd deprecated
#'
#' @return a list of data and priors
#' @export
#'
format_driors <-
  function(taxa = "gadus morhua",
           initial_b = 1,
           initial_b_sd = 0.1,
           terminal_b = NA,
           terminal_b_sd = 0.1,
           carry = NA,
           carry_sd = 0.1,
           u_v_umsy = NA,
           u_years = NA,
           u_sd = 0.1,
           final_u = NA,
           final_u_sd = NA,
           catch = NA,
           years = NA,
           index = NA,
           effort = NA,
           ref_type = "k",
           index_years = 1,
           effort_years = 1,
           use_heuristics = FALSE,
           fmi = c(
             "research" = NA,
             "management" = NA,
             "enforcement" = NA,
             "socioeconomics" = NA
           ),
           fmi_sd = NA,
           sar = NA,
           sar_sd = NA,
           f_sd = 0.1) {
    if (!is.na(sar)) {
      temp <- dplyr::tibble(sar = sar)
      
      pp <-
        rstanarm::posterior_predict(regs$sar_f_reg, newdata = temp)
      
      pp <- pp[pp > quantile(pp, .05) & pp < quantile(pp, 0.95)]
      
      final_u <- c(final_u, exp(mean(pp)))
      
      usd <- ifelse(is.na(sar_sd), sd(pp), sar_sd)
      
      final_u_sd <- c(final_u_sd, usd)
      
    }
    
    
    if (any(!is.na(fmi))) {
      temp <- purrr::map_df(fmi,  ~ .)
      
      pp <-
        rstanarm::posterior_predict(regs$fmi_f_reg, newdata = temp)
      
      pp <- pp[pp > quantile(pp, .05) & pp < quantile(pp, 0.95)]
      
      final_u <- c(final_u, exp(mean(pp)))
      
      final_u_sd <-
        c(final_u_sd, ifelse(is.na(fmi_sd), sd(pp), fmi_sd))
      
      pp <-
        rstanarm::posterior_predict(regs$fmi_b_reg, newdata = temp)
      
      pp <- pp[pp > quantile(pp, .05) & pp < quantile(pp, 0.95)]
      
      terminal_b <- exp(mean(pp))
      
      terminal_b_sd <- ifelse(is.na(fmi_sd), sd(pp), fmi_sd)
      
      ref_type <- "b"
      
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
      c("r", "ln_var")
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
        c("r" = median(FishLifeData$beta_gv[, "r"]),
          "ln_var" = log(0.05))
      
      cov_lh <-
        c("r" = sd(FishLifeData$beta_gv[, "r"]),
          "ln_var" = 0.05)
      
      
    }
    if (is.na(initial_b)) {
      initial_b <- ifelse(ref_type == "k", 1, 2.5)
    }
    
    if (is.na(initial_b)) {
      initial_b <- if (catch[1] / max(catch) < 0.2)
        c(0.7)
      else
        0.4
      
      initial_b_sd <- 0.1
    }
    
    
    if (any(!is.na(final_u))) {
      log_final_u <- log(final_u[!is.na(final_u)])
      
      log_final_u_sd <- final_u_sd[!is.na(final_u)]
      
      if (is.na(terminal_b)){ # if no terminal b prior use approx from U/Umsy
        
        terminal_b <- pmax(.05,2.5 - mean(final_u[!is.na(final_u)]))
        
        terminal_b_sd <- 0.2
        
        if (ref_type == "k"){
          
          ref_type = "b"
          
          initial_b <- initial_b * 2.5
        }
        
      }
      
    } else {
      log_final_u <- NA
      
      log_final_u_sd <- NA
    }
    
    if (use_heuristics == TRUE) {
      ref_type <-  "k"
      
      temp <-
        if (catch[1] / max(catch, na.rm = TRUE) < 0.2)
          0.7
      else
        0.4
      
      initial_b <-  dplyr::case_when(ref_type == "k" ~ temp,
                                     TRUE ~ temp * 2.5)
      
      initial_b_sd <- 0.2
      
      temp_terminal <-
        ifelse((dplyr::last(catch) / max(catch)) > 0.5, 0.6, 0.2)
      
      terminal_b <-
        dplyr::case_when(ref_type == "k" ~ temp_terminal,
                         TRUE ~ temp_terminal * 2.5)
      
      terminal_b_sd <- 0.2
      
    }
    
    driors <-
      list(
        catch = catch,
        years = years,
        carry = carry,
        carry_cv = sqrt(log(carry_sd ^ 2 + 1)),
        terminal_b = terminal_b,
        terminal_b_cv = sqrt(log(terminal_b_sd ^ 2 + 1)),
        initial_b = initial_b,
        initial_b_cv = sqrt(log(initial_b_sd ^ 2 + 1)),
        index = index,
        effort = effort,
        u_v_umsy = u_v_umsy,
        u_years = u_years,
        u_cv = u_sd,
        index_years = index_years,
        effort_years = effort_years,
        growth_rate = mean_lh["r"],
        growth_rate_cv = sqrt(cov_lh["r"]),
        sigma_r = exp(mean_lh["ln_var"]) / 2,
        sigma_r_cv = exp(cov_lh["ln_var"]),
        f_cv = f_sd,
        log_final_u = log_final_u,
        log_final_u_cv = log_final_u_sd
      )
    
    driors$ref_type <- ref_type
    
    return(driors)
    
  } # close function