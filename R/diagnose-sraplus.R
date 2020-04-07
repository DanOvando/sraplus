#' Run minimal diagnostics on sraplus fits
#'
#' @param fit a fitted sraplus model
#' @param driors the driors used to fit the model
#'
#' @return a list of diagnostics
#' @export


diagnose_sraplus <- function(fit, driors) {
  # check quality of fishlife match
  
  # fit <- bayes_fit
  #
  # driors <- ml_driors
  
  diagnostics <- list()
  
  supplied_taxa <-  driors$input_taxa
  
  if (driors$fishlife_taxa != "no_fishlife_match") {
    fishlife_taxa <-
      tolower(paste(
        stringr::str_split(driors$fishlife_taxa, "_", simplify = TRUE)[, 4:5],
        collapse = " "
      ))
    
  } else {
    fishlife_taxa <- driors$fishlife_taxa
  }
  
  diagnostics$fishlife_match <-
    dplyr::case_when(
      supplied_taxa == fishlife_taxa ~ "fishlife matched supplied species",
      fishlife_taxa == "no_fishlife_match" ~ "fishlife couldn't match anything",
      TRUE ~ "fishlife matched supplied genus"
    )
  
  
  
  
  # check model fit
  
  
  if (fit$engine == "sir") {
    diagnostics$distinct_sir_draws <-
      fit$fit$draw_id %>% n_distinct()
    
    
    foo <- function(n_draws, fit) {
      draws <- sample(unique(fit$draw_id), n_draws, replace = FALSE)
      
      fit %>%
        filter(draw_id %in% draws) %>%
        dplyr::group_by(variable, year) %>%
        summarise(mean_value = mean(value))
      
    }
    
    sub_samps <-
      dplyr::tibble(draws = round(seq(1, diagnostics$distinct_sir_draws, length.out = min(20, diagnostics$distinct_sir_draws)))) %>%
      dplyr::mutate(samps = purrr::map(draws, foo, fit = fit$fit)) %>%
      tidyr::unnest(cols = samps)
    
    # sub_samps %>%
    #   ggplot(aes(year, mean_value, color = draws, group = draws)) +
    #   geom_line() +
    #   facet_wrap(~variable, scales = "free_y")
    #
    diagnostics$sir_convergence_plot <- sub_samps %>%
      filter(year == max(year)) %>%
      ggplot(aes(draws, mean_value)) +
      geom_line() +
      facet_wrap(~ variable, scales = "free_y") +
      scale_x_continuous(name = "# of unique SIR samples Used") +
      scale_y_continuous(name = "Mean Terminal Value") +
      theme_sraplus()
    
    
  } else if (fit$engine == "tmb") {
    diagnostics$fit_gradients <- fit$fit$diagnostics
    
    diagnostics$fit_diagnostic_message <-
      fit$fit$fit_diagnostic_message
    
    
  } else {
    diagnostics$divergences <- rstan::get_num_divergent(fit$fit)
    
    diagnostics$max_treedepth <-
      rstan::get_num_max_treedepth(fit$fit)
    
    diagnostics$bfmi <- rstan::get_bfmi(fit$fit)
    
    diagnostics$low_bfmi_chains <-
      rstan::get_low_bfmi_chains(fit$fit)
    
  }
  
  
  return(diagnostics)
  
  
}