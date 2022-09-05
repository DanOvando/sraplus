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
  
  
  tmp <- sraplus::get_prior_posterior(fit = fit, driors = driors, draws = draws, prob = prob)
  
  
  if (tmp$any_fits & any(!is.na(tmp$fits$lower))){
    
    timeseries_plot <-  tmp$fits %>% 
    ggplot2::ggplot() +
    ggplot2::geom_line(aes(year, value, color = metric),size = 1, show.legend = FALSE) +
    ggplot2::geom_point(aes(year, value, color = metric),size = 2, show.legend = FALSE) +
    ggplot2::geom_pointrange(data = tmp$fits %>% dplyr::filter(!is.na(lower)),aes(year, mean, ymin = lower, ymax = upper, color = "Fit"),alpha = 0.5) +
    ggplot2::facet_grid(metric ~ ., scales = "free_y") +
    ggplot2::scale_y_continuous(name = "") +
    ggplot2::labs(x = "Year")  +
    theme_sraplus(base_size = 12)  +
    ggplot2::theme(legend.position = "top")
  } else {
    
   
   timeseries_plot <- tmp$fits %>% 
      ggplot2::ggplot() +
      ggplot2::geom_line(aes(year, value, color = metric),size = 1, show.legend = FALSE) +
      ggplot2::geom_point(aes(year, value, color = metric),size = 2, show.legend = FALSE) +
      ggplot2::facet_grid(metric ~ ., scales = "free_y") +
      ggplot2::scale_y_continuous(name = "") +
      ggplot2::labs(x = "Year")  +
      theme_sraplus(base_size = 12)  +
      ggplot2::theme(legend.position = "top")
    
  }
  
  
    prior_posterior_plot <- tmp$prior_posterior %>%
    ggplot2::ggplot() +
    ggplot2::geom_pointrange(aes(
      x = variable,
      y = mean,
      ymin = lower,
      ymax = upper,
      color = source
    ),
    position = ggplot2::position_dodge(width = 0.1)) +
    ggplot2::facet_wrap( ~ variable, scales = "free") +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = "top")

  patchwork::wrap_plots(timeseries_plot,
                        prior_posterior_plot, widths = c(1, 2))
  
  
}