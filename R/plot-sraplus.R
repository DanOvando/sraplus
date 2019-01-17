#' Plot timeseries of estimates from multiple sraplus models
#'
#' @param ... any number of sraplus models
#'
#' @return a ggplot object
#' @export
#'
plot_sraplus <- function(...){
  
  fit_names<- names(list(...))
  
  fits <- list(...)  %>% 
    purrr::set_names(fit_names) 
  
  fitframe <- dplyr::tibble(fit = fit_names,
                            temp = fits) %>% 
    dplyr::mutate(results = purrr::map(temp,"results")) %>% 
    dplyr::select(-temp) %>% 
    tidyr::unnest()
  
  plotvars <- c("b_div_bmsy","c_div_msy","depletion","u_div_umsy")
  
 fitframe %>% 
   dplyr::filter(variable %in% plotvars) %>% 
   dplyr::group_by(variable,fit) %>% 
   dplyr::mutate(year = seq_along(mean)) %>% 
    dplyr::ungroup() %>% 
    ggplot2::ggplot() + 
   ggplot2::geom_ribbon(aes(year, ymin = lower, ymax = upper, fill = fit),
                size = 0.5, alpha = 0.5) +
   ggplot2::geom_line(aes(year, mean, color = fit),
              size = 1) +
   ggplot2::facet_wrap(~variable, scales = "free_y") + 
    sraplus::theme_sraplus() + 
   ggplot2::scale_y_continuous( name = "") +
    labs(x = "Time") +
   ggplot2::scale_fill_discrete(name = "Fit") + 
   ggplot2::scale_color_discrete(name = "Fit")
  
}