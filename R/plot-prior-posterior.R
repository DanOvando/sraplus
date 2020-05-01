#' plot_prior_posterior
#'
#' @param fit 
#' @param driors 
#' @param draws 
#' @param prob 
#'
#' @return
#' @export
#'
plot_prior_posterior <- function(fit, driors,
                                 draws = 1000,
                                 prob  = 0.9){
  
  
#   library(tidyverse)
#   library(tmbstan)
#   library(sraplus)
#   example_taxa <- "gadus morhua"
#
#    draws = 1000
#
#    prob = 0.9
#    set.seed(42)
#
#    sim <-
#      sraplus_simulator(
#        sigma_proc = 0.05,
#        sigma_u = 0.05,
#        q_slope = 0.05,
#        r = 0.2,
#        years = 25,
#        q = 1e-3,
#        m = 1.01,
#        init_u_umsy = 0.75
#      )
#
#    sim$pop %>%
#      select(year, depletion,catch, effort,u) %>%
#      gather(metric, value, -year) %>%
#      ggplot(aes(year, value)) +
#      geom_point() +
#      facet_wrap(~metric, scales = "free_y") +
#      labs(y = "Value", x = "Year") +
#      sraplus::theme_sraplus()
#
#    effort_years <- seq(5,25, by = 2)
#
#    effort_years <- 1:nrow(sim$pop)
#
#
#   cpue_driors <- format_driors(taxa = example_taxa,
#                                         catch = sim$pop$catch,
#                                         years = sim$pop$year,
#                                         effort = sim$pop$effort[effort_years],
#                                         effort_years = effort_years,
#                                u = sim$pop$u_umsy[effort_years],
#                                u_years = effort_years,
#                                         growth_rate_prior = 0.4,
#                                         growth_rate_prior_cv = 0.1,
#                                         shape_prior = 1.01,
#                                         q_slope_prior = 0.025,
#                                         q_slope_prior_cv = 0.25,
#                                sar = 0,
#                                sar_cv = 0.1,
#                                fmi = c("research" = 0, "management" = 0, "socioeconomics" = 0, 'enforcement' = 0),
#                                         f_ref_type = "f")
#   #
# #
#   cpue_driors  <- format_driors(taxa = example_taxa,
#                                            catch = sim$pop$catch,
#                                            years = sim$pop$year,
#                                            growth_rate_prior = 0.4,
#                                            growth_rate_prior_cv = 0.1,
#                                            shape_prior = 1.01,
#                                            q_slope_prior = 0,
#                                            q_slope_prior_cv = 0.25,
#                                            sar = 2,
#                                            sar_cv = 0.1,
#                                            fmi = c("research" = 0.5, "management" = 0.5, "socioeconomics" = 0.5, 'enforcement' = 0.5),
#                                            f_ref_type = "f",
#                                 terminal_state = 0.5)
#
#   plot_driors(cpue_driors)
#   #
#   cpue_fit  <- fit_sraplus(driors = cpue_driors,
#                                                engine = "sir",
#                                                model = "sraplus_tmb",
#                                                adapt_delta = 0.9,
#                                                max_treedepth = 10,
#                                                n_keep = 2000,
#                                                chains = 1,
#                                                cores = 1,
#                                                estimate_qslope = FALSE,
#                                                estimate_proc_error = FALSE)
#   #
#   plot_sraplus(cpue_fit, years = sim$pop$year)
#
#   #
#   fit <- cpue_fit
#
#   driors <- cpue_driors


 

  timeseries <- dplyr::tibble(year = driors$years,
                              catch = driors$catch)
  
  fitseries <- as.data.frame(matrix(NA, nrow = 0, ncol = 5))
  
  colnames(fitseries) <- c("metric","year", "mean", 'lower','upper')
  any_fits <- TRUE
  if (any(!is.na(driors$index))) {
    index_frame <- dplyr::tibble(year = driors$index_years,
                                 index = driors$index)
    timeseries <- timeseries %>%
      dplyr::left_join(index_frame, by = "year")
    
    index_fit <- fit$results %>%
      dplyr::filter(variable == "index_hat_t") %>%
      dplyr::select(mean, lower, upper) %>%
      dplyr::mutate(metric = "index",
                    year = driors$years) %>%
      filter(year %in% driors$index_years) %>%
      mutate(meanmean = mean(mean),
             sdmean = sd(mean)) %>%
      dplyr::mutate(
        mean = (mean - meanmean) / sdmean,
        lower = (lower - meanmean) / sdmean,
        upper = (upper - meanmean) / sdmean
      ) %>%
      dplyr::select(metric, year, mean, lower, upper)
    
    fitseries <- rbind(fitseries, index_fit)
    
  }
  
  if (any(!is.na(driors$u))) {
    u_frame <- dplyr::tibble(year = driors$u_years,
                             u = driors$u)
    timeseries <- timeseries %>%
      dplyr::left_join(u_frame, by = "year")
    
    u_fit <- fit$results %>% 
      dplyr::filter(variable == "u_div_umsy") %>% 
      dplyr::select(mean,lower, upper) %>% 
      dplyr::mutate(metric = "u",
                    year = driors$years) %>% 
      # dplyr::filter(year %in% driors$u_years) %>% 
      dplyr::mutate(
                    meanmean = mean(mean),
                    sdmean = sd(mean)) %>% 
      dplyr::mutate(mean = (mean - meanmean) / sdmean,
                    lower = (lower - meanmean) / sdmean,
                    upper = (upper - meanmean) / sdmean) %>% 
      dplyr::select(metric, year, mean, lower, upper)
    
    fitseries <- rbind(fitseries, u_fit)
    
  }
  
  if (any(!is.na(driors$effort))) {
    e_frame <- dplyr::tibble(year = driors$effort_years,
                             effort = driors$effort)
    timeseries <- timeseries %>%
      dplyr::left_join(e_frame, by = "year") %>% 
      mutate(cpue = catch / effort)
    
    cpue_fit <- fit$results %>% 
      dplyr::filter(variable == "index_hat_t") %>% 
      dplyr::select(mean,lower, upper) %>% 
      dplyr::mutate(metric = "cpue",
                    year = driors$years) %>% 
      filter(year %in% driors$effort_years) %>% 
      mutate(
                    meanmean = mean(mean),
                    sdmean = sd(mean)) %>% 
      dplyr::mutate(mean = (mean - meanmean) / sdmean,
                    lower = (lower - meanmean) / sdmean,
                    upper = (upper - meanmean) / sdmean) %>% 
      dplyr::select(metric, year, mean, lower, upper)
    
    fitseries <- rbind(fitseries, cpue_fit)
    
    
  }
  
  if (nrow(fitseries) == 0){
    
    any_fits <- FALSE
    
    fitseries <- as.data.frame(matrix(NA, nrow = nrow(timeseries), ncol = 3))
    
    colnames(fitseries) <- c("mean", 'lower','upper')
    
  }
  
  timeseries_plot <- timeseries %>%
    tidyr::gather(metric, value,-year) %>% 
    dplyr::group_by(metric) %>% 
    dplyr::mutate(value = scale(value)) %>% 
    dplyr:: ungroup() %>% {
      if (any_fits){
        dplyr::left_join(.,fitseries, by = c("year","metric"))
      } else {
        dplyr::bind_cols(.,fitseries)
      }
    } 
  
  
  if (any_fits){
    timeseries_plot <- timeseries_plot %>% 
    ggplot2::ggplot() +
    ggplot2::geom_line(aes(year, value, color = metric),size = 1, show.legend = FALSE) +
    ggplot2::geom_point(aes(year, value, color = metric),size = 2, show.legend = FALSE) +
    ggplot2::geom_pointrange(aes(year, mean, ymin = lower, ymax = upper, color = "Fit"),alpha = 0.5) +
    ggplot2::facet_grid(metric ~ ., scales = "free_y") +
    ggplot2::scale_y_continuous(name = "") +
    ggplot2::labs(x = "Year")  +
    theme_sraplus(base_size = 12)  +
    theme(legend.position = "top")
  } else {
    
    timeseries_plot <- timeseries_plot %>% 
      ggplot2::ggplot() +
      ggplot2::geom_line(aes(year, value, color = metric),size = 1, show.legend = FALSE) +
      ggplot2::geom_point(aes(year, value, color = metric),size = 2, show.legend = FALSE) +
      ggplot2::facet_grid(metric ~ ., scales = "free_y") +
      ggplot2::scale_y_continuous(name = "") +
      ggplot2::labs(x = "Year")  +
      theme_sraplus(base_size = 12)  +
      theme(legend.position = "top")
    
  }
  
  
  if (length(driors$log_terminal_u) > 1) {
    driors$log_terminal_u1 = driors$log_terminal_u[1]
    
    driors$log_terminal_u1_cv = driors$log_terminal_u_cv[1]
    
    driors$log_terminal_u2 = driors$log_terminal_u[2]
    
    driors$log_terminal_u2_cv = driors$log_terminal_u_cv[2]
    
    driors <-
      purrr::list_modify(driors,
                         "log_terminal_u" = NULL,
                         "log_terminal_u_cv" = NULL)
    
  }
  
  
  vars <- names(driors)
  
  var_count <- table(stringr::str_replace_all(vars, "_cv", ""))
  
  plot_vars <- names(var_count)[var_count == 2]
  
  has_values <-
    purrr::map_lgl(plot_vars, ~ !all(is.na(driors[[.x]])))
  
  plot_vars <- plot_vars[has_values]
  
  
  
  foo <- function(var, driors, n = 1000) {
    if (stringr::str_detect(var, "log")) {
      sims <- rnorm(n, driors[[var]], driors[[paste0(var, "_cv")]])
      
      simquants <-
        quantile(sims, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2))
      
      sims <- sims[sims >= simquants[1] & sims <= simquants[2]]
      
    } else {
      sims <-
        exp(rnorm(n, log(pmax(driors[[var]], 1e-6)), driors[[paste0(var, "_cv")]]))
      
      simquants <-
        quantile(sims, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2))
      
      sims <- sims[sims >= simquants[1] & sims <= simquants[2]]
      
      
    }
  }
  
  priors <- dplyr::tibble(variable = plot_vars) %>%
    dplyr::mutate(draws = purrr::map(variable, foo, driors = driors)) %>%
    tidyr::unnest(cols = draws)
  
 # this is going to be very unelegant given differnces in naming
  
  fits <- fit$results %>% 
    mutate(variable = dplyr::case_when(variable == "r" ~ "growth_rate",
                            variable == "m" ~ "shape",
                            TRUE ~ variable ))
  
  priors$variable <- stringr::str_replace_all(priors$variable,"_prior","") 
  
  priors <- priors %>% 
    dplyr::group_by(variable) %>% 
    dplyr::summarise(mean = mean(draws),
              lower = quantile(draws,probs = c((1 - prob) / 2)),
              upper = quantile(draws,probs = c(1 - (1 - prob) / 2))) %>% 
    dplyr::mutate(source = "Prior")
  
  if ("initial_state" %in% priors$variable) {
    
    if (driors$b_ref_type == "b") {
      
      initial_state <- fits[fits$variable == "b_div_bmsy",]
      
      initial_state <- as.data.frame(initial_state[1,])
      
      initial_state$variable <-  "initial_state"
      
      fits <- rbind(fits, (initial_state))
      
      
    } else if (driors$b_ref_type == "k"){
      
      initial_state <- fits[fits$variable == "depletion",]
      
      initial_state <- as.data.frame(initial_state[1,])
      
      initial_state$variable <-  "initial_state"
      
      fits <- rbind(fits, (initial_state))
      
      
    }
    
    
  }
  
  if ("terminal_state" %in% priors$variable){
    
    if (driors$b_ref_type == "b") {
      
      terminal_state <- fits[fits$variable == "b_div_bmsy",]
    
      terminal_state <- as.data.frame(terminal_state[nrow(terminal_state),])
        
      terminal_state$variable <-  "terminal_state"
        
    
   fits <- rbind(fits, (terminal_state))
    
    
  } else if (driors$b_ref_type == "k"){
    
    terminal_state <- fits[fits$variable == "depletion",]
    
    terminal_state <- as.data.frame(terminal_state[nrow(terminal_state),])
    
    terminal_state$variable <-  "terminal_state"
    
    fits <- rbind(fits, terminal_state)
    
    
  }
    
    
  } # close if terminal state
    
    
  
  posteriors <- fits %>% 
    dplyr::filter(variable %in% priors$variable) %>% 
    dplyr::mutate(source = "Posterior")
  
  prior_posterior <- posteriors[,colnames(priors)] %>% 
    rbind(priors %>% filter(variable %in% posteriors$variable))
  
  prior_posterior_plot <- prior_posterior %>%
    ggplot2::ggplot() +
    ggplot2::geom_pointrange(aes(
      variable,
      y = mean,
      ymin = lower,
      ymax = upper,
      color = source
    ),
    position = position_dodge(width = 0.1)) +
    ggplot2::facet_wrap( ~ variable, scales = "free") +
    # ggplot2::coord_flip() +
    ggplot2::theme_classic() + 
    theme(legend.position = "top")
  
  
  patchwork::wrap_plots(timeseries_plot,
                        prior_posterior_plot, widths = c(1, 2))
  
  
  
  
  
}