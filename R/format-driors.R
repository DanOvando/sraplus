#' Format data and priors for sraplus
#'
#' @param taxa the genus and species of the species
#' @param initial_b b reference point in the initial year
#' @param initial_b_cv cv associated with initial b reference point
#' @param terminal_b b reference point in the terminal year
#' @param terminal_b_cv cv associated with terminal b reference point
#' @param carry prior on carrying capacity
#' @param carry_cv cv associated with prior on carrying capacity
#' @param u_v_umsy u/umsy data over time
#' @param u_years years in which u/umsy data are available
#' @param u_cv cv associated with u/umsy data
#' @param final_u vector of priors on u/umsy in the terminal years
#' @param final_u_cv vector of cvs on u/umsy in the terminal years
#' @param catch vector of catches over lifetime of fishery
#' @param years vector of years that the catch data correspond to
#' @param index vector of an abundance index
#' @param effort vector of an effort series
#' @param ref_type k if initial and final depletions are in depletion units, b if in b/bmsy units
#' @param index_years the years in which abundance index data are available
#' @param effort_years years in which effort data are available
#' @param use_heuristics logical,TRUE uses catch-msy hueristics for priors, FALSE requires user to pass them
#' @param f_cv no idea, at all
#'
#' @return a list of data and priors
#' @export
#'
#' @examples
#' #' \dontrun{
#' sum("a")
#' }
format_driors <-
  function(taxa,
           initial_b = 1,
           initial_b_cv = 0.1,
           terminal_b = 0.25,
           terminal_b_cv = 0.1,
           carry = NA,
           carry_cv = 0.1,
           u_v_umsy = NA,
           u_years = NA,
           u_cv = 0.1,
           final_u = NA,
           final_u_cv = NA,
           catch = NA,
           years = NA,
           index = NA,
           effort = NA,
           ref_type = "k",
           index_years = 1,
           effort_years = 1,
           use_heuristics = FALSE,
           f_cv = 0.1) {


    if (use_heuristics == TRUE){


      temp <-
        if (catch[1] / max(catch, na.rm = TRUE) < 0.2)
          0.7
      else
        0.4

      initial_b <-  dplyr::case_when(ref_type == "k" ~ temp,
                                   TRUE ~ temp * 2)

      initial_b_cv <- 0.1

      temp_terminal <-
        ifelse((last(catch) / max(catch)) > 0.5, 0.6, 0.2)

      terminal_b <-  dplyr::case_when(ref_type == "k" ~ temp_terminal,
                                   TRUE ~ temp_terminal * 2)

      terminal_b_cv <- 0.1

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

    taxa_location <- grep(sp, FishLifeData$ParentChild_gz[, "ChildName"])[1]

    mean_lh <- FishLifeData$beta_gv[taxa_location,]

    cov_lh <- FishLifeData$Cov_gvv[taxa_location, ,]


    mean_lh <- mean_lh[which(names(mean_lh) %in% params_mvn)]

    cov_lh <-
      cov_lh[which(rownames(cov_lh) %in% params_mvn), which(colnames(cov_lh) %in% params_mvn)] %>% diag()

    } else {


      mean_lh <-  c("r" = median(FishLifeData$beta_gv[,"r"]),"ln_var" = log(0.05))

      cov_lh <- c("r" = sd(FishLifeData$beta_gv[,"r"]),"ln_var" = 0.05)


    }
    if (is.na(initial_b)){

      initial_b <- 1

      ref_type <- "k"

    }


    # if (is.na(initial_b) | is.na(terminal_b)){
    #
    #   ref_type <- "k"
    #
    # }
    #
    if (is.na(initial_b)){

      initial_b <-if (catch[1] / max(catch) < 0.2) c(0.7) else 0.4

      initial_b_cv <- 0.5
    }
    #
    # if (is.na(terminal_b)){
    #
    #   terminal_b = if(last(catch)/max(catch) > 0.5) c(0.5) else c(0.205)
    #
    # }




    driors <-
      list(
        catch = catch,
        years = years,
        carry = carry,
        carry_cv = sqrt(log(carry_cv ^ 2 + 1)),
        terminal_b = terminal_b,
        terminal_b_cv = sqrt(log(terminal_b_cv ^ 2 + 1)),
        initial_b = initial_b,
        initial_b_cv = sqrt(log(initial_b_cv ^ 2 + 1)),
        index = index,
        effort = effort,
        u_v_umsy = u_v_umsy,
        u_years = u_years,
        u_cv = u_cv,
        index_years = index_years,
        effort_years = effort_years,
        growth_rate = mean_lh["r"],
        growth_rate_cv = (cov_lh["r"]),
        sigma_r =exp(mean_lh["ln_var"])/2,
        sigma_r_cv = exp(cov_lh["ln_var"]),
        f_cv = f_cv,
        log_final_u = ifelse(any(!is.na(final_u)),log(final_u[!is.na(final_u)]), NA),
        log_final_u_cv = ifelse(any(!is.na(final_u)),final_u_cv[!is.na(final_u)], NA)
      )

    driors$ref_type <- ref_type

    return(driors)

  } # close function