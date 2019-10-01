summarize_sralpus <- function(fit, output = "table"){
  
  # library(tidyverse)
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
  #                                u_v_umsy = sim$pop$u_umsy[effort_years],
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
  #   
  #
    # cpue_driors  <- format_driors(taxa = example_taxa,
    #                                          catch = sim$pop$catch,
    #                                          years = sim$pop$year,
    #                                          growth_rate_prior = 0.4,
    #                                          growth_rate_prior_cv = 0.1,
    #                                          shape_prior = 1.01,
    #                                          q_slope_prior = 0,
    #                                          q_slope_prior_cv = 0.25,
    #                                          sar = 2,
    #                                          sar_cv = 0.1,
    #                                          fmi = c("research" = 0.5, "management" = 0.5, "socioeconomics" = 0.5, 'enforcement' = 0.5),
    #                                          f_ref_type = "f",
    #                               terminal_state = 0.5)

    # plot_driors(cpue_driors)
    # #
    # cpue_fit  <- fit_sraplus(driors = cpue_driors,
    #                                              engine = "tmb",
    #                                              model = "sraplus_tmb",
    #                                              adapt_delta = 0.9,
    #                                              max_treedepth = 10,
    #                                              n_keep = 2000,
    #                                              chains = 1,
    #                                              cores = 1,
    #                                              estimate_qslope = FALSE,
    #                                              estimate_proc_error = FALSE)
    # 

    result_summary <- cpue_fit$results %>% 
      dplyr::group_by(variable) %>% 
      filter(year == max(year)) 
    
    result_summary_plot <- result_summary %>% 
      ggplot2::ggplot(aes(variable, mean, ymin = lower, ymax = upper)) + 
      ggplot2::geom_pointrange() + 
      ggplot2::facet_wrap(~variable, scales = "free") + 
      ggplot2::labs(x = '', y = "Estimate", caption = "Point is mean estimated value, bars are 95% interval. Values are in the most recent year") + 
      theme_sraplus() + 
      ggplot2::theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
      )  
    
    if (output == "table"){
      
      out <- result_summary
      
    } else {
      out <- result_summary_plot
      
    }
    
  
  return(out)
}